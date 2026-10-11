#!/usr/bin/env bash
# Resume installation from the downloaded cache after mounting the Toshiba image.
set -euo pipefail
cd "$(dirname "$0")/.."
mountpoint -q /mnt/gate-runtime
python3 -u scripts/install_gate_image.py
bash scripts/build_gate.sh
