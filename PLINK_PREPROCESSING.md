# PLINK preprocessing for survival GWAS

This branch prepares the synthetic discovery and replication cohorts with PLINK 2, then validates the resulting files against the existing R QC analysis. It is an explicit Snakemake target, `plink_preprocessing`. The previously reported Cox and GATE results continue to use their original R-derived PCs; the new branch supplies inputs for a separately specified survival analysis.

## Inputs and analysis boundary

Run from the repository root. The branch requires `data/raw/{cohort}.vcf.gz`, the corresponding phenotype TSV, and `results/qc/{cohort}_{samples,variants}.tsv`. Snakemake builds missing inputs using the existing simulation and R QC rules. The small study has 2,200 autosomal variants and 3,604 discovery / 2,404 replication records after fault injection. Coordinates are synthetic, not a human genome reference.

The existing sample audit supplies exact-duplicate, within-ancestry heterozygosity-outlier and invalid-phenotype exclusions. PLINK independently applies the sample-missingness filter. Thus this is an R-assisted PLINK branch, not a replacement for those audit decisions. Exact copied-record detection is not a general relatedness screen. Real-cohort applications also require assay-specific sex, ancestry, relatedness, batch, imputation and reference-build checks.

## Execution

Requirements: Python 3.8 or newer (standard library only), PLINK 2, and the existing Snakemake/R environment when upstream outputs need generation. `PLINK2` may specify an executable outside PATH. No installation or download is performed by this branch.

```bash
# WSL/Linux, from the repository root
export PLINK2=/home/vusi/miniforge3/bin/plink2  # recorded workstation; adapt locally
bash run_snakemake.sh plink_preprocessing --dry-run
bash run_snakemake.sh plink_preprocessing

# Direct execution when R QC and raw inputs already exist
python3 scripts/plink_preprocess.py discovery --plink "$PLINK2" --threads 2
python3 scripts/plink_preprocess.py replication --plink "$PLINK2" --threads 2
```

On Windows, enter WSL first, then `cd /mnt/d/survival_gwas_simulation`. On Linux/HPC, set `R_SCRIPT=Rscript` if upstream rules run. The launcher supports additional Snakemake arguments. Outputs are isolated in `results/plink/discovery/` and `results/plink/replication/`; a rerun replaces that branch's outputs. A successful run ends with `validation.json`. The previous success marker is removed before every direct rerun.

## Processing stages and exact commands

The examples below use `c=discovery` and `o=results/plink/$c`. The Python driver creates exclusion lists and calls these stages in order. Each actual command, including executable and resource flags, is saved to `$o/commands.sh`; stdout and PLINK logs are retained per stage. Thresholds are explicitly fixed to the current study design: missingness 0.03/0.02, MAF 0.05, MAC 100, HWE 1e-6 and three PCs. Changes to `config.R` require reviewing this branch's thresholds as well; disagreement with the R QC reference stops validation.

### 1. VCF import and missingness audit

```bash
plink2 --vcf data/raw/$c.vcf.gz --double-id --make-pgen --out $o/raw
plink2 --pfile $o/raw --missing --out $o/raw_missing
```

`--double-id` preserves the full VCF sample name as both FID and IID, including underscores. PGEN stores genotypes; PVAR records variant coordinates and alleles; PSAM identifies samples. `.smiss` and `.vmiss` describe sample and variant missingness. The synthetic VCF contains autosomal hard calls, so there is no dosage-quality or sex-chromosome filter in this run.

### 2. Sample exclusions before variant filtering

```bash
plink2 --pfile $o/raw --remove $o/remove.tsv --mind 0.03 \
  --make-pgen --out $o/samples
```

The driver writes a two-column FID/IID removal list from the existing audit and omits `--remove` if the list is empty. PLINK removes samples exceeding 3% missing genotypes. The driver checks that the remaining sample IDs exactly match the R QC keep set before proceeding. This ordering ensures variant statistics use retained participants.

### 3. Missingness, minor allele frequency and count

```bash
plink2 --pfile $o/samples --geno 0.02 --maf 0.05 --mac 100 \
  --make-pgen --out $o/base
```

Variants with more than 2% missing genotypes, MAF below 5%, or MAC below 100 are removed. These thresholds reproduce the current common-variant simulation design; they are not universal thresholds for a real study.

### 4. Hardy–Weinberg equilibrium within ancestry

```bash
plink2 --pfile $o/base --keep $o/ancestry_0.keep --hardy --out $o/hwe_0
plink2 --pfile $o/base --keep $o/ancestry_1.keep --hardy --out $o/hwe_1
plink2 --pfile $o/base --exclude $o/hwe_exclude.txt --make-pgen --out $o/qc
plink2 --pfile $o/qc --freq --missing --out $o/qc_stats
```

Ancestry keep lists are derived from the retained participants' metadata, sorted by ancestry label. The driver excludes the union of variants with exact HWE p < 1e-6 in either group; it omits `--exclude` when there are no failures. All synthetic participants are founders. HWE is not tested on the pooled ancestry sample for filtering, because population mixture can create departures unrelated to assay error. The resulting variant-ID set and ALT frequencies are compared against the R QC reference.

### 5. LD pruning and principal components

```bash
plink2 --pfile $o/qc --indep-pairwise 250kb 1 0.2 --out $o/ld
plink2 --pfile $o/qc --extract $o/ld.prune.in --pca 3 --out $o/pcs
```

PLINK prunes using a 250 kb sliding window, one-variant step and r-squared threshold 0.2. `.prune.in` identifies variants used for PCA; the full QC marker set is retained for association testing. The exact PCA computation produces `.eigenvec` and `.eigenval`; the driver requires three finite PCs and complete sample-ID alignment.

PLINK's pooled-genotype LD pruning differs from the original R ancestry-centred pruning. The retained pruning sets and PCs need not agree. PC signs are arbitrary. Using these PCs changes the model specification and requires rerunning and evaluating Cox/GATE models; the published association statistics must not be relabelled as PLINK-PC results.

### 6. Survival covariates and GATE genotype export

```bash
plink2 --pfile $o/qc --make-bed --out $o/gate
plink2 --pfile $o/qc --export A --out $o/dosage
```

`survival_covariates.tsv` joins PCs to phenotype records by IID and preserves site, ancestry group, age, centred age, sex indicator, signal/null follow-up time and event status. Time is in years; event coding is 1 for an event and 0 for censoring. Survival outcomes are not taken from the BED companion FAM phenotype field.

The PLINK-produced BED/BIM/FAM files provide the genotype input format for GATE. BIM allele ordering must be read explicitly: the driver decodes every BED call and compares it with the allele-labelled dosage export. It also checks every exported dosage and missing call against the original VCF and checks ALT-frequency agreement within 1e-5. Do not assume BIM A1 or the `.raw` counted allele is always ALT.

For Cox analysis, join genotype dosages, survival outcomes and PCs by sample ID and preserve the desired baseline-hazard strata. For GATE, regenerate its phenotype/covariate design and indexed VCF from a consistently aligned branch, fit the survival null model, and run the association stage. The existing GATE scripts use `data/gate/` and do not automatically consume `results/plink/`.

`plink2 --glm ... firth-fallback` is a case–control association command in the supplied practice transcript. It is not the survival association stage used here: the survival model must account for follow-up time and censoring.

## Validation outputs

`validation.json` is written only after sample and variant keep-set concordance, ALT-frequency concordance, exhaustive VCF dosage comparison, BED decoding, finite PCs and the ordered survival-covariate join all pass. This certifies file integrity and the tested preprocessing specification, not clinical validity or calibrated survival inference. Compact validation summaries are versioned; sample-level files, genotype exports and verbose logs remain local.

### Recorded small-dataset test

Both cohorts completed through Snakemake on 22 October 2026 using PLINK v2.0.0-a.7.8 (19 September 2026), two threads and 2,048 MB requested PLINK memory. Sample and variant sets matched the existing R QC outputs, and all dosage, BED, PCA and covariate-join checks passed.

| Cohort | Input samples | QC samples | QC variants | PCA markers | Genotype cells verified |
|---|---:|---:|---:|---:|---:|
| Discovery | 3,604 | 3,588 | 2,160 | 216 | 7,750,080 |
| Replication | 2,404 | 2,388 | 2,160 | 216 | 5,158,080 |

Each input contained 2,200 markers; 10 variants per cohort failed the ancestry-specific HWE stage after the preceding variant filters. The maximum ALT-frequency discrepancy was below 0.000001, consistent with printed-frequency rounding. Validation summaries: [discovery](results/plink/discovery/validation.json) and [replication](results/plink/replication/validation.json).

## References

- [PLINK 2 sample and variant filtering](https://www.cog-genomics.org/plink/2.0/filter)
- [PLINK 2 Hardy–Weinberg reporting](https://www.cog-genomics.org/plink/2.0/basic_stats#hardy)
- [PLINK 2 LD pruning](https://www.cog-genomics.org/plink/2.0/ld)
- [PLINK 2 principal components](https://www.cog-genomics.org/plink/2.0/strat)
- [PLINK 2 file formats](https://www.cog-genomics.org/plink/2.0/formats)
