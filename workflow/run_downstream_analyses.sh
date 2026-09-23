#!/usr/bin/env bash
# Run the complete custom analysis suite after FEBA has finished.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 RUN_DIRECTORY [T_THRESHOLD]" >&2
  exit 2
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
run_dir="$1"
threshold="${2:-2.5}"
analyses_dir="$run_dir/analyses"
log="$analyses_dir/run_downstream_analyses.log"

[[ -f "$run_dir/feba/fit_logratios.tab" ]] || {
  echo "Missing FEBA output: $run_dir/feba/fit_logratios.tab" >&2
  exit 1
}
mkdir -p "$analyses_dir"

{
  printf 'Running Barcode QC...\n'
  "$repo_dir/workflow/run_barcode_qc.sh" "$run_dir"
  printf 'Running fitness analysis with t-like threshold %s...\n' "$threshold"
  "$repo_dir/workflow/run_fitness.sh" "$run_dir" "$threshold"
  printf 'Running pairwise condition tests...\n'
  "$repo_dir/workflow/run_pairwise_ttest.sh" "$run_dir"
  printf 'Running gene remapping...\n'
  "$repo_dir/workflow/run_gene_remap.sh" "$run_dir"
} > "$log" 2>&1

printf 'Downstream analyses complete:\n  Results: %s\n  Log: %s\n' "$analyses_dir" "$log"
