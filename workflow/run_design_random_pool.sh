#!/usr/bin/env bash
# Build one reusable barcode-to-insertion pool from the standard workspace references.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 WORKSPACE_DIRECTORY" >&2; exit 2
fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_dir="$1"
map_tab="$workspace_dir/reference/mapping/Xylella_ML3.tab"
genes="$workspace_dir/reference/genes/genes"
pool="$workspace_dir/reference/pool/ML3.pool"
for file in "$map_tab" "$genes"; do [[ -f "$file" ]] || { echo "Missing: $file" >&2; exit 1; }; done
mkdir -p "$(dirname "$pool")" "$workspace_dir/logs"
log="$workspace_dir/logs/DesignRandomPool.log"
/usr/bin/perl "$repo_dir/feba/bin/DesignRandomPool.pl" -pool "$pool" -genes "$genes" "$map_tab" \
  > "$log" 2>&1
printf 'DesignRandomPool complete:\n  Results: %s\n  Log: %s\n' "$pool" "$log"
