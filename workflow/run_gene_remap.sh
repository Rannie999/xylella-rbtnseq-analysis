#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 RUN_DIRECTORY" >&2; exit 2
fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; run_dir="$1"
workspace_dir="$(cd "$run_dir/../.." && pwd)"
input_dir="$run_dir/analyses/fitness/csv"
[[ -d "$input_dir" ]] || input_dir="$run_dir/fit_analysis/csv"
[[ -d "$input_dir" ]] || { echo "Fitness CSV directory not found under: $run_dir" >&2; exit 1; }
out_dir="$run_dir/analyses/gene_remap"; mkdir -p "$out_dir"; log="$out_dir/run_gene_remap.log"
python "$repo_dir/scripts/gene_remap.py" \
  --input-dir "$input_dir" \
  --genome-dir "$workspace_dir/reference/genomes" \
  --output-dir "$out_dir" > "$log" 2>&1
printf 'Gene remapping complete:\n  Results: %s\n  Log: %s\n' "$out_dir/new_locus_info.csv" "$log"
