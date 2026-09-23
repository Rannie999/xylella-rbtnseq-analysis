#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: $0 RUN_DIRECTORY" >&2; exit 2; fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; run_dir="$1"
out_dir="$run_dir/analyses/barcode_qc"; mkdir -p "$out_dir"; log="$out_dir/run_barcode_qc.log"
python "$repo_dir/scripts/barcode_qc.py" --input "$run_dir/feba/all.poolcount" \
  --name-map "$run_dir/input/sample_name_map.csv" --output "$out_dir" > "$log" 2>&1
printf 'Barcode QC complete:\n  Results: %s\n  Log: %s\n' "$out_dir" "$log"
