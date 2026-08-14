#!/usr/bin/env bash
# =============================================================================
# generate.sh — synthesize a cohort with more samples than bins
# =============================================================================
#
# RandomizedSVD transposes its input when there are fewer rows than columns,
# and then takes A_t from the matrix as given rather than copying it. Neither
# happens for the 1000G example, which has 15,593 bins against 18 samples, so
# without this the orientation that a selected-bin analysis of a large cohort
# actually runs in has no coverage at all.
#
# The cohort is generated rather than committed: 2,000 files are a lot to carry
# in a repository and nothing about them needs to be inspected by hand. What is
# committed is this script and the checksums of what NGS-PCA makes of it, so a
# change in the numbers fails rather than passing quietly.
#
# Depths come from a MINSTD generator and are formatted by integer arithmetic,
# so every awk produces the same bytes: the products stay under 2^53, where a
# double is exact, and no float is ever rounded for printing.
#
# Usage: generate.sh <outputDir> [samples] [bins]
# =============================================================================

set -euo pipefail

OUT_DIR="${1:?usage: generate.sh <outputDir> [samples] [bins]}"
SAMPLES="${2:-2000}"
BINS="${3:-1500}"

if [ "$BINS" -ge "$SAMPLES" ]; then
  echo "ERROR: this fixture exists to test bins < samples; got $BINS bins, $SAMPLES samples" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

for i in $(seq 1 "$SAMPLES"); do
  awk -v seed="$i" -v bins="$BINS" 'BEGIN {
    # MINSTD: every product stays below 2^53, so double arithmetic is exact here
    s = (seed * 16807) % 2147483647
    for (b = 0; b < bins; b++) {
      s = (s * 16807) % 2147483647
      # 2.00 to 41.99, printed as two integers so no float is rounded
      v = s % 4000 + 200
      printf "chr%d\t%d\t%d\t%d.%02d\n", 1 + int(b * 22 / bins), b * 1000, b * 1000 + 1000,
             int(v / 100), v % 100
    }
  }' | gzip -c > "$OUT_DIR/$(printf 'sample%05d' "$i").regions.bed.gz"
done

echo "wrote $SAMPLES samples x $BINS bins to $OUT_DIR"
