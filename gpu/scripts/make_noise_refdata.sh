#!/usr/bin/env bash
# GS2's own round-off noise for every case: rerun the reference bin/gs2 with a
# different data layout ('yxels', or 'lxyes' for a case that already uses
# yxels). The layout changes only the order of GS2's velocity-space sums, so
# the difference between this run and gpu/refdata/<case>/gs2/<case>.out.nc is
# GS2's round-off spread. tests/test_end_to_end.py uses it as the tolerance for
# variables that are zero up to round-off (e.g. the particle flux with
# adiabatic electrons).
#
# Usage: gpu/scripts/make_noise_refdata.sh [case ...]   (default: all gpu/cases)
# Output: gpu/refdata/<case>/gs2_layout/<case>.out.nc
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GS2="$ROOT/bin/gs2"
CASES_DIR="$ROOT/gpu/cases"
REF="$ROOT/gpu/refdata"

if [ $# -eq 0 ]; then
  set -- $(cd "$CASES_DIR" && ls *.in | sed 's/\.in$//')
fi

for case in "$@"; do
  dir="$REF/$case/gs2_layout"
  rm -rf "$dir"
  mkdir -p "$dir"
  layout=yxels
  if grep -q 'layout = "yxels"' "$CASES_DIR/$case.in"; then
    layout=lxyes
  fi
  sed -E "s/^ layout = \"[a-z]+\"/ layout = \"$layout\"/" "$CASES_DIR/$case.in" > "$dir/$case.in"
  grep -q "layout = \"$layout\"" "$dir/$case.in" || { echo "$case: could not set layout"; exit 1; }
  ( cd "$dir" && mpirun -np 1 "$GS2" "$case.in" > gs2.log 2>&1 ) \
    && echo "$case: ok (layout $layout)" || echo "$case: FAILED (see $dir/gs2.log)"
done
