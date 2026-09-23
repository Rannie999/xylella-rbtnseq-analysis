#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: $0 RUN_DIRECTORY" >&2; exit 2; fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; run_dir="$1"
out_dir="$run_dir/analyses/pairwise_ttest"; mkdir -p "$out_dir"; log="$out_dir/run_pairwise_ttest.log"
python "$repo_dir/scripts/pairwise_ttest.py" --input "$run_dir/feba" \
  --name-map "$run_dir/input/sample_name_map.csv" --groups "$run_dir/input/groups.json" \
  --filtered "$run_dir/analyses/fitness/csv/fit_filtered.csv" \
  --output "$out_dir" > "$log" 2>&1
printf 'Pairwise condition tests complete:\n  Results: %s\n  Log: %s\n' "$out_dir" "$log"
