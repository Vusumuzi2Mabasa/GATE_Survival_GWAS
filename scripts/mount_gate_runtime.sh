#!/usr/bin/env bash
# Run as root under WSL; all runtime contents remain in an image on Toshiba.
set -euo pipefail
image=/mnt/d/survival_gwas_simulation/tools/gate-runtime.ext4
target=/mnt/gate-runtime
if ! test -e "$image"; then
  truncate -s 8G "$image"
  mkfs.ext4 -q -F -L GATE_RUNTIME "$image"
fi
mkdir -p "$target"
if ! mountpoint -q "$target"; then
  mount -o loop,noatime "$image" "$target"
fi
test "$(findmnt -n -o LABEL "$target")" = GATE_RUNTIME
chown 1000:1000 "$target"
df -h "$target"
