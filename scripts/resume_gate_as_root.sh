#!/usr/bin/env bash
set -euo pipefail
cd /mnt/d/survival_gwas_simulation
bash scripts/mount_gate_runtime.sh
runuser -u vusi -- bash scripts/install_gate.sh
