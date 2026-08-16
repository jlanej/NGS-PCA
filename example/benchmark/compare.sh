#!/usr/bin/env bash
# =============================================================================
# compare.sh — wall time and peak RSS, this build against another
# =============================================================================
#
# Indicative, not authoritative. Both numbers move with the machine, the JVM,
# the heap you allow, and above all the shape of the cohort: the decomposition
# is dominated by products whose cost is bins x samples x (numPC + oversample),
# so a result at one shape says little about another. Run it on the hardware
# you care about, at a shape near the one you run.
#
# It checks that both builds produced identical output before reporting any
# timing. A speed comparison between two different answers is not a comparison.
#
# Usage:
#   compare.sh <baselineJar> <currentJar> [samples] [bins]
#
# THREADS, HEAP, NUMPC and OVERSAMPLE come from the environment. NUMPC plus
# OVERSAMPLE is the width of the decomposition, and it is what the QR cost
# scales with ((rows + columns) x width squared) - so a run left at the small
# default says little about a production one at -numPC 500.
#
# To get a baseline jar for upstream:
#   git clone https://github.com/PankratzLab/NGS-PCA.git /tmp/upstream
#   git -C /tmp/upstream checkout 2ffbcfc
#   (cd /tmp/upstream/ngspca && mvn -B package -DskipTests)
#
# Upstream is fed a manifest rather than a directory, because it does not sort
# its directory listing and a different sample order is a different answer.
# Its -outputDir needs the trailing slash this script passes: it builds output
# paths by concatenation.
# =============================================================================

set -euo pipefail

BASELINE_JAR="${1:?usage: compare.sh <baselineJar> <currentJar> [samples] [bins]}"
CURRENT_JAR="${2:?usage: compare.sh <baselineJar> <currentJar> [samples] [bins]}"
SAMPLES="${3:-2000}"
BINS="${4:-1500}"
THREADS="${THREADS:-4}"
HEAP="${HEAP:-4g}"
NUMPC="${NUMPC:-20}"
OVERSAMPLE="${OVERSAMPLE:-60}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

# peak RSS: GNU time reports kilobytes, BSD/macOS time reports bytes
case "$(uname -s)" in
  Darwin) TIME_FLAG="-l" ;;
  *)      TIME_FLAG="-v" ;;
esac

peak_rss_mb() {  # reads a /usr/bin/time log on stdin
  awk '/[Mm]aximum resident set size/ {
         for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+$/) { v = $i; break }
       }
       END { if (v == "") print "?"; else if (v > 100000000) printf "%.0f", v / 1048576;
             else printf "%.0f", v / 1024 }'
}

echo "cohort: ${SAMPLES} samples x ${BINS} bins, width ${NUMPC}+${OVERSAMPLE}, -threads ${THREADS}, -Xmx${HEAP}"
"${HERE}/../wideOrientation/generate.sh" "${WORK}/mosdepth" "${SAMPLES}" "${BINS}" > /dev/null
ls -1 "${WORK}/mosdepth"/*.regions.bed.gz | sort > "${WORK}/manifest.txt"

for side in baseline current; do
  jar="${BASELINE_JAR}"
  [ "${side}" = current ] && jar="${CURRENT_JAR}"
  mkdir -p "${WORK}/${side}"
  start=$(date +%s)
  /usr/bin/time ${TIME_FLAG} java "-Xmx${HEAP}" -jar "${jar}" \
    -input "${WORK}/manifest.txt" -outputDir "${WORK}/${side}/" \
    -numPC "${NUMPC}" -iters 5 -oversample "${OVERSAMPLE}" -randomSeed 42 -threads "${THREADS}" \
    > "${WORK}/${side}.log" 2> "${WORK}/${side}.time"
  printf '%s\t%s\t%s\n' "${side}" "$(( $(date +%s) - start ))" \
    "$(peak_rss_mb < "${WORK}/${side}.time")" >> "${WORK}/results.tsv"
done

echo
for f in svd.pcs.txt svd.loadings.txt svd.singularvalues.txt svd.bins.txt svd.samples.txt; do
  if ! cmp -s "${WORK}/baseline/${f}" "${WORK}/current/${f}"; then
    echo "STOP: ${f} differs between the two builds - the timings below would compare two answers"
    exit 1
  fi
done
echo "output identical across all five files"
echo
printf '%-10s %10s %12s\n' build seconds "peak RSS MB"
awk -F'\t' '{ printf "%-10s %10s %12s\n", $1, $2, $3 }' "${WORK}/results.tsv"
