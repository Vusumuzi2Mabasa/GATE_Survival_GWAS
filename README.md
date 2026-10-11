# Simulated survival GWAS in R

A complete **training workflow** from synthetic genotypes and survival records
through QC, association discovery, independent replication and an R-generated
report. This is not AWI-Gen data, a dense human GWAS, or a production pipeline for
restricted participant data. The earlier one-variant example remains separate.

## Open the result

Open `report/index.html` in a browser after the pipeline finishes. The report uses
**gtsummary**, **flextable** and **ggplot2**. R table objects are saved in
`report/tables/` so you can restyle them without refitting the models.

## Run on this computer

All project files are on **D:\survival_gwas_simulation**, the Toshiba external
drive. The drive must remain connected throughout the run.

```powershell
Set-Location D:\survival_gwas_simulation
.\run_pipeline.ps1 -DryRun
.\run_pipeline.ps1
```

The launcher uses the existing Ubuntu-24.04 WSL installation to run Snakemake and
the existing Windows R 4.5.1 installation to execute every statistical and report
script. It uses **two concurrent jobs** and a **4 GB scheduling budget**, not an
enforced OS memory limit. Snakemake's isolated Python environment is under
`.venv/` on D:. It contains orchestration dependencies, not statistical analyses.
The report package versions are recorded in `logs/environment.tsv` and
`report/tables/package_versions.tsv`.

Alternatively, from WSL:

```bash
cd /mnt/d/survival_gwas_simulation
bash run_snakemake.sh --dry-run
bash run_snakemake.sh
```

To run the identical stages sequentially without Snakemake, from the project root:

```powershell
& 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe' --vanilla scripts/run_all.R
```

That alternative reruns every stage. The Snakemake launcher runs only stages
requiring updates. A changed config causes downstream regeneration. Copy the
project or archive outputs before trying a new configuration: result filenames
are stable and rerunning replaces them. All R entrypoints use `--vanilla` so
personal startup files do not affect execution.

## Workflow

```text
Environment check -> simulate independent discovery and replication cohorts
                      |                       |
                 sample/variant QC       sample/variant QC
                      |                       |
                 LD pruning + PCA        LD pruning + PCA
                      |
              22-chromosome Cox scans
              (signal and null outcomes)
                      |
              multiple-testing summaries
                      |
              LD clumping -> test discovery leads in replication
                      |
              PH diagnostics + site sensitivity
                      |
              computational checks + reveal simulation truth
              + 20 repeated null-outcome calibration runs
                      |
              gtsummary / flextable / ggplot2 HTML report
```

## Prespecified simulation

The settings are in `config.R`.

* 3,600 discovery and 2,400 replication participants before fault injection.
* Six fictional sites, two fictional ancestry groups within sites; all people
  are generated independently. These do not represent named African populations.
* 2,200 autosomal markers across 22 chromosomes; 10 LD blocks per chromosome,
  10 markers per block. Both haplotypes are generated independently, with markers
  copying a block anchor with probability 0.93 or redrawing from its population
  frequency. Chromosomes and positions are illustrative, not a real assembly.
* Three direct-effect variants on chromosomes 3, 10 and 18, with conditional HRs
  1.50, 0.67 and 1.40. These are deliberately large teaching effects. The other
  variants in their blocks may associate through LD without a direct effect.
* Weibull proportional hazards with site, ancestry, age and sex effects. Exact
  event times and independent dropout, with administrative censoring at 5 years.
* A separate null outcome has the same covariate structure but no genetic effect.
* Technical faults: high missingness, excessive heterozygosity, exact duplicate
  records, monomorphic/rare markers, and HWE-deviating markers.
* The two cohorts share the generating parameters but have independent samples,
  genotypes, outcomes and random seeds. No tuning uses replication outcomes.

The simulator exports RDS matrices, phenotype TSVs and ordinary gzip-compressed
VCFs (not bgzip/tabix indexed). The analysis uses the matching RDS matrices for
speed and to avoid an extra VCF parser dependency. This starts from simulated
genotypes, not raw sequencing reads. No imputation or real reference panel is
claimed.

## Statistical specification

For marker j, the scan fits:

```r
coxph(Surv(time, event) ~ dosage_j + age10 + male +
        PC1 + PC2 + PC3 + strata(site, ancestry_group),
      ties = "efron")
```

The effect allele is always ALT (`G`) and dosage counts 0, 1 or 2 copies.
`age10` is centred age per ten years. Population group and site enter as separate
baseline hazard strata, not as a causal biological interpretation of ancestry.
Genotype PCs are obtained without outcome information. Mean-filled genotypes are
used only for PCA/LD calculations. Association models use per-variant complete
cases and report their individual sample and event counts.

The ordinary Cox **two-sided Wald test** is used for common variants (MAF >= 5%,
MAC >= 100) with abundant events in unrelated simulated participants. Warnings,
singular fits and nonfinite estimates are reported rather than silently accepted.
The scanner is not designed for rare variants, severe sparsity, general familial
relatedness or millions of dense markers. At that scale, select and validate an
appropriate tool such as SPACox or GATE, rather than assuming ordinary full Cox
fits remain appropriate.

### QC and population structure

1. Validate sample/variant identifiers and genotype coding.
2. Remove samples with >3% missing genotypes, robust within-population
   heterozygosity outliers, exact copied records or invalid phenotype records.
3. Remove markers with >2% missingness, MAF <5%, MAC <100, or exact HWE p <1e-6
   in either homogeneous simulated population.
4. Prune local LD at r-squared 0.2 within 250 kb. Genotypes are centred within
   population for LD calculations, so between-population frequency differences
   are not mistaken for local linkage.
5. Compute three PCs from globally variance-standardised pruned genotypes in
   each cohort separately. Inspect the PCA plots.

These are fixed teaching thresholds, not universal prescriptions. Sex checks,
array batch effects, reference-allele alignment, imputation quality and general
kinship are not simulated. The duplicate detector finds exact copies only;
production data need a validated relatedness method. Broad group strata do not
replace modelling fine-scale structure in real African cohorts.

### Multiple testing and replication

* The primary threshold is **0.05 / all QC-passing variants**, including any
  variants whose model subsequently fails. It controls the testing burden of
  this simulated marker panel if its p-values are calibrated.
* The conventional **5e-8** threshold is separately labelled. Do not equate a
  2,200-marker training scan with a densely imputed human GWAS.
* BH q-values are exploratory supplementary results.
* Select leads by discovery p-value, then clump within 250 kb at ancestry-centred
  r-squared >=0.2. Clumping is not fine-mapping or proof of independence.
* Test only these leads in independently generated replication data. Require the
  same effect direction and two-sided p <0.05 / number of discovery leads.
  Leads failing replication QC still count toward that denominator.
* Record global and term-specific Schoenfeld-residual diagnostics for lead fits,
  and discovery site-specific sensitivity estimates. Their nominal p-values are
  exploratory; selected-model diagnostics are not confirmatory tests.

### Validation and interpretation

`scripts/validate.R` checks known HWE examples, allele-flip invariance, direct
Cox-fit agreement, sample/variant alignment, QC exclusions, full scan coverage,
multiple-testing arithmetic and the absence of ground-truth reads in discovery
scripts. The simulator's causal labels are only used afterward in evaluation.

The signal scan contains omitted causal markers in each one-variant model;
marginal HR estimates can differ from conditional generating coefficients and
may show nonproportionality. A significant tag SNP is not proof of a causal SNP.
The null scan is the cleaner calibration check. Lambda GC and nominal rejection
rates are descriptive: LD makes tests dependent, and one scan cannot establish
repeated-sampling type-I error or power. Outcomes are never permuted across
incompatible strata to manufacture a null.

A supplementary calibration experiment regenerates 20 independent null outcomes
on the fixed discovery genotype sample and scans the LD-pruned marker panel.
It reports each repetition's nominal rejection fraction, the mean and its Monte
Carlo standard error, and a t-based interval across repetitions. It also reports
family-wise error and an exact binomial interval for that smaller panel, using
its own Bonferroni threshold. This is a coarse, conditional calibration check;
it does not estimate full-GWAS power or precisely establish full-GWAS family-wise
error. The number of repetitions is set in the Snakefile null_calibration rule.

Independent replication here means an independent sample from the same
simulation model. It is not validation across new real-world populations.

## Main outputs

| Location | Contents |
|---|---|
| `report/index.html` | Main R-generated report |
| `report/tables/` | HTML tables, gtsummary/flextable R objects, package versions |
| `report/plots/` | ggplot2 PNG figures |
| `results/qc/` | Every sample and variant QC decision |
| `results/pca/` | PCs, pruned markers, PCA objects |
| `results/summary/` | Full association results and calibration summaries |
| `results/replication/` | Leads, LD clumps, replication, PH and site diagnostics |
| `results/validation/` | Checks, known-truth recovery and source checksums |
| `logs/` | Per-stage logs and R session information |
| `data/truth/` | Simulator truth, excluded from discovery analysis |

## Resources and scaling

This pilot runs on the laptop. A real whole-genome scan or a large Monte Carlo
study may need a cluster. The rules already split association scans by chromosome;
adapting a cluster executor requires checking the institute's scheduler, software
environment, storage policy and approved analysis method. No cluster access or
institutional data is used here.

On a different Linux machine, install the R packages listed above, create an
isolated Snakemake environment, and set `R_SCRIPT=Rscript`. The local wrapper's
default Windows R path is specific to this computer. Installed Python package
versions are recorded in `environment-lock.txt`; R versions in the report.

## Word results report

The pipeline also creates `report/Survival_GWAS_Results.docx` in R using
`officer`, `gtsummary`, and `flextable`, with figures generated by `ggplot2`.
It contains editable text and tables; embedded plots can be regenerated from R.
No Pandoc is needed. Install `officer` alongside the other reporting packages.

To rebuild only the Word document after the reporting outputs exist, run from
the project root:

```powershell
& 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe' --vanilla scripts/report_word.R
```

The Snakemake workflow and `scripts/run_all.R` both include the Word stage.

## GATE extension

GATE 0.45 is installed in an isolated Linux runtime image on this Toshiba drive.
Both 2160-variant survival scans completed and passed output validation. See
[GATE_PILOT.md](GATE_PILOT.md) for execution commands, build changes, and limits.
The R-authored supplement is `report/GATE_Survival_GWAS_Results.docx`.
Signal: 21 panel Bonferroni hits. Null: zero hits. Both scans show deflation;
the non-LOCO run is an implementation trial, not a production calibration claim.

## Method references

* [survival::coxph](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html)
* [gwasurvivr](https://bioconductor.org/packages/release/bioc/html/gwasurvivr.html)
* [GATE survival GWAS method](https://www.nature.com/articles/s41467-022-32885-x)
* [PLINK HWE and population structure](https://www.cog-genomics.org/plink/2.0/tutorials/qc1)
* [PLINK PCA](https://www.cog-genomics.org/plink/2.0/strat)
* [Snakemake installation and Windows/WSL support](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* [gtsummary conversion to flextable](https://www.danieldsjoberg.com/gtsummary/reference/as_flex_table.html)
