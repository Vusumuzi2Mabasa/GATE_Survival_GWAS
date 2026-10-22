#!/usr/bin/env python3
"""Audited PLINK2 preprocessing branch for the synthetic survival cohorts."""
import argparse, csv, gzip, json, math, os, shlex, subprocess
from pathlib import Path


def read(path):
    with open(path) as f:
        return list(csv.DictReader((line for line in f if not line.startswith('##')), delimiter='\t'))


def write(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cohort', choices=['discovery', 'replication'])
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--plink', default=os.environ.get('PLINK2', 'plink2'))
    a = ap.parse_args(); c = a.cohort
    if a.threads < 1: ap.error('--threads must be positive')
    out = Path('results/plink') / c; out.mkdir(parents=True, exist_ok=True)
    # A failed rerun must not retain a success marker.
    (out / 'validation.json').unlink(missing_ok=True)
    ph = read(f'data/raw/{c}_phenotypes.tsv'); audit = read(f'results/qc/{c}_samples.tsv')
    variants = read(f'results/qc/{c}_variants.tsv')
    assert len({r['id'] for r in ph}) == len(ph)
    assert {r['id'] for r in ph} == {r['id'] for r in audit}
    # Clinical validity, within-ancestry heterozygosity and exact-copy decisions
    # are explicitly inherited from the existing R sample audit. PLINK tests missingness.
    remove = [r['id'] for r in audit if r['heterozygosity_outlier']=='TRUE' or r['exact_duplicate']=='TRUE' or r['valid_phenotype']!='TRUE']
    (out/'remove.tsv').write_text(''.join(f'{s}\t{s}\n' for s in remove))
    version = subprocess.check_output([a.plink, '--version'], text=True).strip()
    (out/'version.txt').write_text(version+'\n')
    with open(out/'commands.sh','w') as commands:
        commands.write('#!/usr/bin/env bash\nset -euo pipefail\n')
        def run(name, *args):
            cmd=[a.plink, '--threads',str(a.threads),'--memory','2048',*map(str,args),'--out',str(out/name)]
            commands.write(shlex.join(cmd)+'\n'); commands.flush()
            with open(out/(name+'.stdout.txt'),'w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        run('raw','--vcf',f'data/raw/{c}.vcf.gz','--double-id','--make-pgen')
        run('raw_missing','--pfile',out/'raw','--missing')
        args=['--pfile',out/'raw']
        if remove: args += ['--remove',out/'remove.tsv']
        run('samples',*args,'--mind','0.03','--make-pgen')
        samples=read(out/'samples.psam'); ids=[r['IID'] for r in samples]
        assert set(ids)=={r['id'] for r in audit if r['keep']=='TRUE'}, 'Sample QC differs from R audit'
        phen={r['id']:r for r in ph}
        run('base','--pfile',out/'samples','--geno','0.02','--maf','0.05','--mac','100','--make-pgen')
        failed=set()
        for i, group in enumerate(sorted({phen[s]['ancestry_group'] for s in ids})):
            keep=out/f'ancestry_{i}.keep'
            keep.write_text(''.join(f'{s}\t{s}\n' for s in ids if phen[s]['ancestry_group']==group))
            run(f'hwe_{i}','--pfile',out/'base','--keep',keep,'--hardy')
            for r in read(out/f'hwe_{i}.hardy'):
                p=float(r['P']); assert math.isfinite(p)
                if p < 1e-6: failed.add(r['ID'])
        (out/'hwe_exclude.txt').write_text(''.join(s+'\n' for s in sorted(failed)))
        args=['--pfile',out/'base']
        if failed: args+=['--exclude',out/'hwe_exclude.txt']
        run('qc',*args,'--make-pgen')
        run('qc_stats','--pfile',out/'qc','--freq','--missing')
        run('ld','--pfile',out/'qc','--indep-pairwise','250kb','1','0.2')
        run('pcs','--pfile',out/'qc','--extract',out/'ld.prune.in','--pca','3')
        run('gate','--pfile',out/'qc','--make-bed')
        run('dosage','--pfile',out/'qc','--export','A')
    pvar=read(out/'qc.pvar'); vids=[r['ID'] for r in pvar]
    assert set(vids)=={r['variant'] for r in variants if r['keep']=='TRUE'}, 'Variant QC differs from R audit'
    af={r['ID']:float(r['ALT_FREQS']) for r in read(out/'qc_stats.afreq')}
    expected={r['variant']:r for r in variants}
    maxerr=max(abs(af[v]-float(expected[v]['effect_allele_frequency'])) for v in vids)
    assert maxerr < 1e-5
    # Validate every exported dosage and its counted allele against the original VCF.
    raw=read(out/'dosage.raw'); assert [r['IID'] for r in raw]==ids
    checks=0; vset=set(vids)
    with gzip.open(f'data/raw/{c}.vcf.gz','rt') as f:
        for line in f:
            if line.startswith('##'): continue
            row=line.rstrip().split('\t')
            if row[0]=='#CHROM': indexes={s:i+9 for i,s in enumerate(row[9:])}; continue
            if row[2] not in vset: continue
            v=row[2]; ref,alt=row[3:5]
            assert ref==expected[v]['ref'] and alt==expected[v]['alt']
            col=next(k for k in raw[0] if k.startswith(v+'_'))
            allele=col[len(v)+1:]; assert allele in [ref,alt]
            for j,s in enumerate(ids):
                gt=row[indexes[s]].split(':')[0].replace('|','/')
                val=raw[j][col]
                if '.' in gt: assert val=='NA'
                else:
                    dosage=sum(int(g) for g in gt.split('/'))
                    assert float(val)==(dosage if allele==alt else 2-dosage)
                checks+=1
    pcs=read(out/'pcs.eigenvec'); assert [r['IID'] for r in pcs]==ids
    for r in pcs: assert all(math.isfinite(float(r[k])) for k in ['PC1','PC2','PC3'])
    joined=[]
    for r in pcs:
        d=dict(phen[r['IID']]); d.update({k:r[k] for k in ['PC1','PC2','PC3']})
        joined.append(d)
    write(out/'survival_covariates.tsv',joined,list(ph[0])+['PC1','PC2','PC3'])
    # Decode the PLINK-produced BED as an independent check on the GATE handoff.
    bim=[line.split() for line in (out/'gate.bim').read_text().splitlines()]
    fam=[line.split()[1] for line in (out/'gate.fam').read_text().splitlines()]
    assert fam==ids and [r[1] for r in bim]==vids
    bed=(out/'gate.bed').read_bytes(); assert bed[:3]==bytes([108,27,1])
    stride=(len(ids)+3)//4; assert len(bed)==3+stride*len(vids)
    for j,b in enumerate(bim):
        v=b[1]; col=next(k for k in raw[0] if k.startswith(v+'_'))
        counted=col[len(v)+1:]
        for i,r in enumerate(raw):
            code=(bed[3+j*stride+i//4] >> (2*(i%4))) & 3
            val=[2,None,1,0][code] # BED counts BIM A1
            if r[col]=='NA': assert val is None
            else: assert val==(float(r[col]) if counted==b[4] else 2-float(r[col]))
    summary=dict(cohort=c,plink_version=version,raw_samples=len(ph),raw_variants=len(variants),retained_samples=len(ids),retained_variants=len(vids),pruned_variants=len((out/'ld.prune.in').read_text().splitlines()),hwe_excluded=len(failed),genotype_cells_checked=checks,max_alt_frequency_error=maxerr,sample_ids_match_R=True,variant_ids_match_R=True,vcf_dosages_match=True,bed_roundtrip_matches=True,finite_pcs=True,survival_join_matches=True)
    (out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
