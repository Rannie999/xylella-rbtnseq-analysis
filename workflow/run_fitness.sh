#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then echo "Usage: $0 RUN_DIRECTORY [T_THRESHOLD]" >&2; exit 2; fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; run_dir="$1"; threshold="${2:-2.5}"
out_dir="$run_dir/analyses/fitness"; mkdir -p "$out_dir"; log="$out_dir/run_fitness.log"
python "$repo_dir/scripts/fit_analysis.py" --input "$run_dir/feba" --output "$out_dir" \
  --name-map "$run_dir/input/sample_name_map.csv" --groups "$run_dir/input/groups.json" --threshold "$threshold" > "$log" 2>&1
printf 'Fitness analysis complete:\n  Results: %s\n  Log: %s\n' "$out_dir" "$log"
