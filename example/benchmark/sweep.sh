#!/usr/bin/env bash
# =============================================================================
# sweep.sh — compare.sh across a few shapes, so one result is not mistaken for
#            a general one
# =============================================================================
#
# Two axes decide which regime a run is in, and they matter more than size:
#
#   orientation  bins > samples takes one code path, samples > bins the other
#                (the input is transposed, and the products go the other way)
#   width        numPC + oversample. The QR costs rows x width squared and the
#                products cost bins x samples x width, so widening the
#                decomposition shifts which one dominates
#
# The shapes below cross those two. They are small enough to run on a laptop in
# a few minutes, which also means none of them is a production cohort - see the
# README for what that does and does not tell you.
#
# Usage: sweep.sh <baselineJar> <currentJar>
# =============================================================================

set -euo pipefail

BASELINE_JAR="${1:?usage: sweep.sh <baselineJar> <currentJar>}"
CURRENT_JAR="${2:?usage: sweep.sh <baselineJar> <currentJar>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${THREADS:-8}"

# label is one field - read splits on whitespace
# label              samples  bins  numPC  oversample
SHAPES="
bins>samples          1200    3000    100    400
samples>bins          3000    1200    100    400
near-square           2000    2000    100    400
samples>bins/narrow   3000    1200     20     80
bins>samples/narrow   1200    3000     20     80
"

printf '%-20s %12s %7s %9s %9s %9s %9s %8s\n' \
  shape "samples x bins" width "base s" "now s" "base MB" "now MB" speedup
printf '%s\n' "--------------------------------------------------------------------------------------------"

echo "${SHAPES}" | while read -r label samples bins numpc oversample; do
  [ -z "${label:-}" ] && continue
  out="$(NUMPC="${numpc}" OVERSAMPLE="${oversample}" THREADS="${THREADS}" \
         "${HERE}/compare.sh" "${BASELINE_JAR}" "${CURRENT_JAR}" "${samples}" "${bins}" 2>&1)" || {
    printf '%-20s %12s %7s   FAILED\n' "${label}" "${samples}x${bins}" "$((numpc + oversample))"
    echo "${out}" | tail -3
    continue
  }
  bs=$(echo "${out}" | awk '$1 == "baseline" { print $2 }')
  br=$(echo "${out}" | awk '$1 == "baseline" { print $3 }')
  cs=$(echo "${out}" | awk '$1 == "current"  { print $2 }')
  cr=$(echo "${out}" | awk '$1 == "current"  { print $3 }')
  printf '%-20s %12s %7s %9s %9s %9s %9s %7sx\n' \
    "${label}" "${samples}x${bins}" "$((numpc + oversample))" "${bs}" "${cs}" "${br}" "${cr}" \
    "$(awk -v a="${bs}" -v b="${cs}" 'BEGIN { if (b > 0) printf "%.1f", a / b; else print "-" }')"
done
