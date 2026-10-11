#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mountpoint -q /mnt/gate-runtime
export MAMBA_ROOT_PREFIX=/mnt/gate-runtime/mamba
eval "$(tools/micromamba/bin/micromamba shell hook --shell bash)"
micromamba activate /mnt/gate-runtime/env
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
Rscript --vanilla scripts/check_gate_null.R signal
Rscript --vanilla tools/GATE-main/extdata/step1_fitNULLGLMM.R \
  --plinkFile=data/gate/discovery --phenoFile=data/gate/discovery_phenotypes.tsv \
  --phenoCol=null_event --eventTimeCol=null_time --traitType=survival \
  --sampleIDColinphenoFile=id \
  --covarColList=age10,male,PC1,PC2,PC3,siteSite_2,siteSite_3,siteSite_4,siteSite_5,siteSite_6,ancestry_groupPop_B \
  --nThreads=2 --LOCO=FALSE --maxiter=100 \
  --outputPrefix=results/gate/discovery_null > logs/gate/step1_null.log 2>&1
Rscript --vanilla scripts/check_gate_null.R null
for outcome in signal null; do
  for chrom in {1..22}; do
    Rscript --vanilla tools/GATE-main/extdata/step2_SPAtests.R \
      --vcfFile=data/gate/discovery.vcf.gz --vcfFileIndex=data/gate/discovery.vcf.gz.csi \
      --vcfField=GT --sampleFile=data/gate/discovery.samples --chrom="$chrom" \
      --GMMATmodelFile="results/gate/discovery_${outcome}.rda" \
      --varianceRatioFile="results/gate/discovery_${outcome}.varianceRatio.txt" \
      --SAIGEOutputFile="results/gate/chr${chrom}_${outcome}_full.tsv" \
      --minMAF=0.05 --minMAC=100 --LOCO=FALSE \
      > "logs/gate/step2_chr${chrom}_${outcome}_full.log" 2>&1
  done
done
