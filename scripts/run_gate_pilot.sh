#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
mountpoint -q /mnt/gate-runtime
ENV=/mnt/gate-runtime/env
export PATH="$ENV/bin:$PATH"
export MAMBA_ROOT_PREFIX=/mnt/gate-runtime/mamba
eval "$(tools/micromamba/bin/micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
mkdir -p results/gate logs/gate
Rscript --vanilla -e 'library(GATE); print(packageVersion("GATE")); print(sessionInfo())' > logs/gate/session.txt
bgzip -c data/gate/discovery.vcf > data/gate/discovery.vcf.gz
tabix -f -C -p vcf data/gate/discovery.vcf.gz
cut -f1 data/gate/discovery_phenotypes.tsv | tail -n +2 > data/gate/discovery.samples
# Non-LOCO pilot: inspect variance component and calibration before expansion.
Rscript --vanilla tools/GATE-main/extdata/step1_fitNULLGLMM.R \
  --plinkFile=data/gate/discovery \
  --phenoFile=data/gate/discovery_phenotypes.tsv \
  --phenoCol=event --eventTimeCol=time --traitType=survival \
  --sampleIDColinphenoFile=id \
  --covarColList=age10,male,PC1,PC2,PC3,siteSite_2,siteSite_3,siteSite_4,siteSite_5,siteSite_6,ancestry_groupPop_B \
  --nThreads=2 --LOCO=FALSE --maxiter=100 \
  --outputPrefix=results/gate/discovery_signal > logs/gate/step1.log 2>&1
Rscript --vanilla scripts/check_gate_null.R > logs/gate/null_checks.log 2>&1
Rscript --vanilla tools/GATE-main/extdata/step2_SPAtests.R \
  --vcfFile=data/gate/discovery.vcf.gz --vcfFileIndex=data/gate/discovery.vcf.gz.csi \
  --vcfField=GT --sampleFile=data/gate/discovery.samples --chrom=3 \
  --GMMATmodelFile=results/gate/discovery_signal.rda \
  --varianceRatioFile=results/gate/discovery_signal.varianceRatio.txt \
  --SAIGEOutputFile=results/gate/chr3_signal.tsv \
  --minMAF=0.05 --minMAC=100 --LOCO=FALSE \
  > logs/gate/step2_chr3.log 2>&1
