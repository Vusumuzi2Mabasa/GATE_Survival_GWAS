#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export R_SCRIPT="${R_SCRIPT:-/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec .venv/bin/snakemake --cores 2 --resources mem_mb=4000 --rerun-incomplete "$@"
