#!/usr/bin/env bash
# Combine all .codes files for one sequencing experiment into one poolcount matrix.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 CODES_DIRECTORY" >&2; exit 2
fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
codes_dir="$1"
workspace_dir="$(cd "$codes_dir/../.." && pwd)"
experiment="$(basename "$codes_dir")"
pool="$workspace_dir/reference/pool/ML3.pool"
prefix="$workspace_dir/poolcounts/$experiment/$experiment"
[[ -f "$pool" ]] || { echo "Missing: $pool" >&2; exit 1; }
shopt -s nullglob
codes=("$codes_dir"/*.codes)
(( ${#codes[@]} )) || { echo "No .codes files in $codes_dir" >&2; exit 1; }
mkdir -p "$(dirname "$prefix")"
log="$(dirname "$prefix")/combineBarSeq.log"
/usr/bin/perl "$repo_dir/feba/bin/combineBarSeq.pl" -all "$prefix" "$pool" "${codes[@]}" \
  > "$log" 2>&1
printf 'combineBarSeq complete:\n  Results: %s\n  Log: %s\n' "$prefix.poolcount" "$log"
