#!/usr/bin/env bash
# Estimate strain and gene fitness from a prepared FEBA directory.
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: $0 RUN_DIRECTORY" >&2; exit 2; fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; run_dir="$1"; feba_dir="$run_dir/feba"
Rscript_bin="${RSCRIPT:-Rscript}"
if ! command -v "$Rscript_bin" >/dev/null 2>&1; then
  echo "Rscript not found: $Rscript_bin. Activate the R environment or set RSCRIPT to its full path." >&2
  exit 1
fi
[[ -f "$feba_dir/all.poolcount" ]] || { echo "Missing: $feba_dir/all.poolcount" >&2; exit 1; }
log="$feba_dir/runFEBA.log"
"$Rscript_bin" "$repo_dir/feba/bin/RunFEBA.R" Xylella "$feba_dir" "$repo_dir/feba" > "$log" 2>&1
python "$repo_dir/scripts/summarize_feba_run.py" --run_dir "$run_dir" >> "$log" 2>&1
printf 'RunFEBA complete:\n  Results: %s\n  Log: %s\n' "$feba_dir" "$log"
