#!/usr/bin/env bash
# Build all.poolcount and other FEBA inputs for one biological run.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 POOLCOUNT_DIRECTORY INPUT_DIRECTORY" >&2; exit 2
fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
poolcount_dir="$(cd "$1" && pwd)"
workspace_dir="$(cd "$poolcount_dir/../.." && pwd)"
experiment="$(basename "$poolcount_dir")"
input_dir="$(cd "$2" && pwd)"
run_dir="$(cd "$input_dir/.." && pwd)"
exps="$input_dir/exps.tsv"; out_dir="$run_dir/feba"
set_name="$(awk -F '\t' '
  NR > 1 && $1 != "" { names[$1] = 1 }
  END {
    for (name in names) { value = name; count++ }
    if (count == 1) print value; else exit 1
  }
' "$exps")" || { echo "exps.tsv must contain exactly one non-empty SetName" >&2; exit 1; }
[[ "$set_name" == "$experiment" ]] || {
  echo "SetName '$set_name' in $exps must match poolcount prefix '$experiment'" >&2
  exit 1
}
genes="$workspace_dir/reference/genes/genes.GC"
pool="$workspace_dir/reference/pool/ML3.pool"
for file in "$poolcount_dir/$experiment.poolcount" "$genes" "$pool"; do [[ -f "$file" ]] || { echo "Missing: $file" >&2; exit 1; }; done
mkdir -p "$out_dir"
log="$out_dir/BarSeqR.log"
/usr/bin/perl "$repo_dir/feba/bin/BarSeqR.pl" -org Xylella -indir "$poolcount_dir" \
  -exps "$exps" -genesfile "$genes" -poolfile "$pool" -outdir "$out_dir" -noR \
  > "$log" 2>&1
printf 'BarSeqR complete:\n  Results: %s\n  Log: %s\n' "$out_dir/all.poolcount" "$log"
