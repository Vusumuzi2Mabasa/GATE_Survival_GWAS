import os
import glob

# Native Windows R is used from WSL here so existing R report packages are reused.
# Override R_SCRIPT on Linux/HPC, where the usual value is Rscript.
R = os.environ.get("R_SCRIPT", "Rscript")
shell.executable("/bin/bash")
shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; ")
COHORTS = ["discovery", "replication"]
SCENARIOS = ["signal", "null"]
CHRS = [f"{i:02d}" for i in range(1, 23)]
BASE = ["config.R", "scripts/common.R"]

rule all:
    input:
        "report/Survival_GWAS_Results.docx",
        "report/index.html",
        "report/tables/tables.html",
        "results/validation/checks.tsv",
        "results/validation/recovery.tsv",
        "results/validation/null_calibration_summary.tsv"

rule preflight:
    input: "scripts/preflight.R"
    output: "logs/environment.tsv"
    log: "logs/preflight.log"
    shell: "{R:q} --vanilla scripts/preflight.R > {log:q} 2>&1"

rule simulate:
    input: BASE + ["scripts/simulate.R", "logs/environment.tsv"]
    output:
        expand("data/raw/{cohort}.rds", cohort=COHORTS),
        expand("data/raw/{cohort}.vcf.gz", cohort=COHORTS),
        expand("data/raw/{cohort}_phenotypes.tsv", cohort=COHORTS),
        "data/raw/variants.tsv", "data/truth/variants.tsv",
        expand("data/truth/{cohort}_samples.tsv", cohort=COHORTS)
    log: "logs/simulate.log"
    resources: mem_mb=2200
    shell: "{R:q} --vanilla scripts/simulate.R > {log:q} 2>&1"

rule qc:
    input: raw="data/raw/{cohort}.rds", code=BASE + ["scripts/qc.R"]
    output:
        "results/qc/{cohort}.rds", "results/qc/{cohort}_samples.tsv",
        "results/qc/{cohort}_variants.tsv", "results/qc/{cohort}_summary.tsv"
    log: "logs/qc_{cohort}.log"
    resources: mem_mb=2200
    shell: "{R:q} --vanilla scripts/qc.R {wildcards.cohort:q} > {log:q} 2>&1"

rule pca:
    input: qc="results/qc/{cohort}.rds", code=BASE + ["scripts/pca.R"]
    output:
        "results/pca/{cohort}.tsv", "results/pca/{cohort}_pruned_variants.tsv",
        "results/pca/{cohort}_model.rds", "results/pca/{cohort}.png"
    log: "logs/pca_{cohort}.log"
    resources: mem_mb=2200
    shell: "{R:q} --vanilla scripts/pca.R {wildcards.cohort:q} > {log:q} 2>&1"

rule scan:
    input:
        qc="results/qc/discovery.rds", pca="results/pca/discovery.tsv",
        code=BASE + ["scripts/scan.R"]
    output: "results/scan/{scenario}/chr{chrom}.tsv"
    wildcard_constraints: scenario="signal|null", chrom="[0-9]{2}"
    log: "logs/scan_{scenario}_{chrom}.log"
    resources: mem_mb=800
    shell: "{R:q} --vanilla scripts/scan.R {wildcards.scenario:q} {wildcards.chrom:q} > {log:q} 2>&1"

rule summarise:
    input:
        scans=lambda w: expand("results/scan/{scenario}/chr{chrom}.tsv", scenario=w.scenario, chrom=CHRS),
        code=BASE + ["scripts/summarise.R"]
    output:
        "results/summary/{scenario}_associations.tsv", "results/summary/{scenario}_metrics.tsv",
        "results/summary/{scenario}_manhattan.png", "results/summary/{scenario}_qq.png"
    log: "logs/summary_{scenario}.log"
    resources: mem_mb=800
    shell: "{R:q} --vanilla scripts/summarise.R {wildcards.scenario:q} > {log:q} 2>&1"

rule replicate:
    input:
        "results/summary/signal_associations.tsv",
        expand("results/qc/{cohort}.rds", cohort=COHORTS),
        expand("results/pca/{cohort}.tsv", cohort=COHORTS),
        BASE + ["scripts/replicate.R"]
    output:
        "results/replication/lead_results.tsv", "results/replication/clumps.tsv",
        "results/replication/ph_diagnostics.tsv", "results/replication/site_sensitivity.tsv"
    log: "logs/replication.log"
    resources: mem_mb=1600
    shell: "{R:q} --vanilla scripts/replicate.R > {log:q} 2>&1"

rule null_calibration:
    input:
        "results/qc/discovery.rds", "results/pca/discovery.tsv",
        "results/pca/discovery_pruned_variants.tsv", BASE + ["scripts/null_calibration.R"]
    output:
        "results/validation/null_repetitions.tsv", "results/validation/null_calibration_summary.tsv"
    log: "logs/null_calibration.log"
    params: repetitions=20
    resources: mem_mb=1000
    shell: "{R:q} --vanilla scripts/null_calibration.R {params.repetitions} > {log:q} 2>&1"

rule validate:
    input:
        expand("results/summary/{scenario}_associations.tsv", scenario=SCENARIOS),
        "results/replication/lead_results.tsv",
        expand("results/qc/{cohort}_samples.tsv", cohort=COHORTS),
        expand("results/qc/{cohort}_variants.tsv", cohort=COHORTS),
        expand("data/truth/{cohort}_samples.tsv", cohort=COHORTS),
        "data/truth/variants.tsv", "Snakefile", "config.R", sorted(glob.glob("scripts/*.R"))
    output:
        "results/validation/checks.tsv", "results/validation/truth_comparison.tsv",
        "results/validation/recovery.tsv", "results/validation/source_checksums.tsv"
    log: "logs/validate.log"
    resources: mem_mb=1600
    shell: "{R:q} --vanilla scripts/validate.R > {log:q} 2>&1"

rule report:
    input:
        "results/validation/checks.tsv", "results/validation/recovery.tsv",
        "results/validation/null_calibration_summary.tsv", "results/validation/null_repetitions.tsv",
        expand("results/summary/{scenario}_metrics.tsv", scenario=SCENARIOS),
        expand("results/summary/{scenario}_associations.tsv", scenario=SCENARIOS),
        expand("results/qc/{cohort}_summary.tsv", cohort=COHORTS),
        expand("results/pca/{cohort}.tsv", cohort=COHORTS),
        "results/replication/lead_results.tsv", "results/replication/ph_diagnostics.tsv",
        BASE + ["scripts/report.R"]
    output:
        "report/index.html", "report/tables/tables.html",
        "report/tables/cohort_gtsummary.rds", "report/tables/flextables.rds",
        "report/tables/package_versions.tsv",
        "report/plots/signal_manhattan.png", "report/plots/null_manhattan.png",
        "report/plots/qq.png", "report/plots/pca.png", "report/plots/survival.png",
        "report/plots/replication.png",
        "report/plots/null_calibration.png", "report/tables/numbers_at_risk.tsv"
    log: "logs/report.log"
    resources: mem_mb=1600
    shell: "{R:q} --vanilla scripts/report.R > {log:q} 2>&1"

rule word_report:
    input:
        "report/index.html", "report/tables/flextables.rds",
        "report/tables/cohort_gtsummary.rds", "scripts/report_word.R",
        expand("report/plots/{plot}.png", plot=["pca", "survival", "signal_manhattan", "replication", "qq", "null_calibration"]),
        expand("results/summary/{scenario}_metrics.tsv", scenario=SCENARIOS),
        "results/replication/lead_results.tsv", "results/replication/ph_diagnostics.tsv",
        "results/validation/checks.tsv", "results/validation/null_calibration_summary.tsv"
    output: "report/Survival_GWAS_Results.docx"
    log: "logs/report_word.log"
    resources: mem_mb=1600
    shell: "{R:q} --vanilla scripts/report_word.R > {log:q} 2>&1"

# Optional PLINK preprocessing branch. Existing Cox/GATE analyses retain their R PCs.
rule plink_preprocessing:
    input:
        expand("results/plink/{cohort}/validation.json", cohort=COHORTS)

rule plink_cohort:
    input:
        vcf="data/raw/{cohort}.vcf.gz",
        pheno="data/raw/{cohort}_phenotypes.tsv",
        samples="results/qc/{cohort}_samples.tsv",
        variants="results/qc/{cohort}_variants.tsv",
        code="scripts/plink_preprocess.py"
    output:
        "results/plink/{cohort}/validation.json",
        "results/plink/{cohort}/survival_covariates.tsv",
        "results/plink/{cohort}/qc.pgen",
        "results/plink/{cohort}/qc.pvar",
        "results/plink/{cohort}/qc.psam",
        "results/plink/{cohort}/gate.bed",
        "results/plink/{cohort}/gate.bim",
        "results/plink/{cohort}/gate.fam",
        "results/plink/{cohort}/pcs.eigenvec",
        "results/plink/{cohort}/pcs.eigenval",
        "results/plink/{cohort}/ld.prune.in",
        "results/plink/{cohort}/commands.sh",
        "results/plink/{cohort}/version.txt"
    wildcard_constraints: cohort="discovery|replication"
    log: "logs/plink_{cohort}.log"
    threads: 2
    resources: mem_mb=2400
    shell: "python3 scripts/plink_preprocess.py {wildcards.cohort:q} --threads {threads} > {log:q} 2>&1"
