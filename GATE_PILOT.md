# GATE pilot setup

Status: GATE 0.45 installed and passed its package-load check. The chromosome 3
pilot and both full 2160-variant scans completed on synthetic data. Both null
models converged and each scan passed 13 output checks. Signal: 21 panel
Bonferroni hits, 13 at 5e-8, lambda GC 0.713. Null: zero hits, lambda GC 0.874.
Required SPA tests converged (56 signal, 61 null). Deflation is present; this is
an implementation trial, not proof of calibration. The Word supplement is
`report/GATE_Survival_GWAS_Results.docx`; R source is `scripts/report_gate.R`.

## Prepared inputs

Run `Rscript --vanilla scripts/prepare_gate.R` from this project directory.
The script exports discovery QC genotypes to `data/gate/discovery.bed`,
`.bim`, `.fam`, and `.vcf`, plus `discovery_phenotypes.tsv`.
The BED export is decoded and compared with every original genotype,
including missingness. See `data/gate/validation.txt` for validation status.
These are synthetic data, not AWI Gen data.

## Installation plan

Use Ubuntu WSL and a dedicated library/environment under this Toshiba project.
Do not replace the system R library. The preflight found R and a C++ compiler,
but no GATE, Docker, Conda or PLINK on PATH. Missing R dependencies include
RcppParallel, SPAtest, RcppEigen, SKAT and MetaSKAT. Also check for bgzip/tabix
and required C++ libraries before compilation. Pin the exact GATE source commit
and record all package versions. Compilation and load checks succeeded. The GitHub source archive is tools/GATE-main.zip,
pinned to commit 962f1028120c73e24ec34f29e02eeae3bdc3bfb3.

Direct installation through WSL onto NTFS was prohibitively slow. The runtime
now lives inside `tools/gate-runtime.ext4`, an 8 GB Linux filesystem image on
the Toshiba drive, mounted in WSL at `/mnt/gate-runtime`. This does not move the
project to the internal drive. After restarting WSL, mount it before running:

```powershell
wsl -d Ubuntu-24.04 -u root --exec bash /mnt/d/survival_gwas_simulation/scripts/resume_gate_build.sh
```

The wrapper mounts the image and keeps WSL running during compilation. For the
pilot after a successful build, use the same command with `resume_gate_pilot.sh`.
An isolated mount command can be lost when WSL shuts down between commands.

The environment uses R 4.2 and Python 3.10. The build-only compatibility patch
aligns GATE's DLL lookup with its package name and selects GNU C11 for legacy C;
completed native dependencies are retained. BGEN build scripts are converted
to Python 3 syntax, its Waf generator terminates using return, and the obsolete
stream position typedef uses std::streampos. Legacy implicit C declarations
are warnings. Backups and change records are in `logs/gate`. Statistical
routines are not intentionally altered. Installation is only successful after
the package-load check completes, not merely because the source is present.

The first trial is configured in scripts/run_gate_pilot.sh for chromosome 3.
The matched fixed-effect Cox comparator is scripts/gate_cox_comparator.R;
its 100-variant output is results/gate/chr3_matched_cox.tsv. It uses the same
site and ancestry indicators and mean-imputed missing genotypes as the planned
GATE trial. It is not itself a GATE result.

Full execution: `scripts/resume_gate_full.sh` mounts the runtime and executes
`scripts/run_gate_full.sh`. Summarise and validate in R with
`Rscript --vanilla scripts/summarise_gate.R`, then render the Word supplement
with `Rscript --vanilla scripts/report_gate.R`. Scan files use whitespace;
validated combined files and summary tables in results/gate use tabs.
The VCF regional reader requires a CSI index; TBI indexing failed to open the
region. Sample IDs, positions, alleles and frequencies passed full output checks.

## Analysis plan after installation

1. Verify GATE's packaged example and installed function arguments.
2. Compress and index the exported VCF with bgzip and tabix.
3. Fit the survival null model using event=1, censored=0, and observed time
   in years. Preserve raw event times rather than binning at one year.
4. Adjust for age, sex, three PCs and explicit site/ancestry indicators.
   Construct the genetic relationship matrix from suitable QC markers;
   inspect variance-component estimation and convergence before testing.
5. Test a small chromosome subset, then all eligible variants if successful.
   Verify effect-allele orientation and record SPA and convergence diagnostics.
6. Compare with Cox using matched covariates, alongside the original stratified
   Cox results. Repeat for the null outcome before assessing calibration.
7. Add verified results and limitations to R tables and the Word report.

The current simulation contains unrelated individuals and a small marker panel.
It can test compatibility and implementation, but cannot demonstrate GATE's
biobank-scale performance or its advantages with related individuals. GATE
is a survival mixed-model association method, not a Bayesian posterior sampler.

## Primary sources

- https://github.com/weizhou0/GATE
- https://github.com/weizhou0/GATE/blob/main/DESCRIPTION
- https://www.nature.com/articles/s41467-022-32885-x
