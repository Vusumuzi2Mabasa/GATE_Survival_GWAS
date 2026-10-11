#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
mountpoint -q /mnt/gate-runtime
export MAMBA_ROOT_PREFIX=/mnt/gate-runtime/mamba
ENV=/mnt/gate-runtime/env
SRC=/mnt/gate-runtime/source/GATE-main
export TMPDIR=/mnt/gate-runtime/tmp PIP_CACHE_DIR=/mnt/gate-runtime/pip-cache
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"
eval "$(tools/micromamba/bin/micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null MAKEFLAGS=-j1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p logs/gate config
printf 'CXXFLAGS=-O1 -g0\nCXX11FLAGS=-O1 -g0\nCXX14FLAGS=-O1 -g0\n' > config/gate.Makevars
export R_MAKEVARS_USER="$ROOT/config/gate.Makevars"
if ! python -c 'import cget' >/dev/null 2>&1; then
  python -m pip install cget > logs/gate/cget-install.log 2>&1
fi
Rscript --vanilla -e 'if (!requireNamespace("lintools", quietly=TRUE)) {options(timeout=600); install.packages("lintools", repos="https://cloud.r-project.org")}; stopifnot(requireNamespace("lintools", quietly=TRUE))' > logs/gate/lintools-install.log 2>&1
python scripts/gate_build_compat.py "$SRC"
find "$SRC/thirdParty/bgen" -name wscript -print0 | xargs -0 python -m lib2to3 -w -n > logs/gate/waf-python3.log 2>&1
chmod +x "$SRC/configure" "$SRC/thirdParty/bgen/waf"
R CMD INSTALL "$SRC" > logs/gate/package-build.log 2>&1
Rscript --vanilla -e 'library(GATE); print(packageVersion("GATE")); print(sessionInfo())' > logs/gate/package-load.txt 2>&1
micromamba list -p "$ENV" --explicit > logs/gate/environment-explicit.txt
printf 'GATE build and package-load check succeeded.\n'
