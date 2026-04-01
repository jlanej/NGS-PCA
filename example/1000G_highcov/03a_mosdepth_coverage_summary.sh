#!/usr/bin/env bash
# =============================================================================
# 03a_mosdepth_coverage_summary.sh — Compute per-sample autosomal coverage
#                                      statistics from mosdepth output
# =============================================================================
#
# Reads mosdepth *.regions.bed.gz files and computes summary statistics
# (mean, median, SD, MAD, IQR) of per-bin coverage across autosomal
# chromosomes (chr1–chr22) for each sample.
#
# If bedtools is available and BED_EXCLUDE points to a valid file, the script
# also computes "HQ" (high-quality) variants of each statistic using only
# autosomal bins that do NOT overlap the exclusion regions.  The HQ median
# can be used downstream to estimate mitochondrial DNA copy number (mtDNA CN).
#
# Output:
#   $QC_OUTPUT/mosdepth_coverage_summary.tsv
#   Columns: SAMPLE_ID  MEAN_COV  MEDIAN_COV  SD_COV  MAD_COV  IQR_COV
#            HQ_MEAN_COV  HQ_MEDIAN_COV  HQ_SD_COV  HQ_MAD_COV  HQ_IQR_COV
#
# The script auto-detects MOSDEPTH_DIR from config.sh and parallelizes across
# all available CPU cores using xargs.  Each sample is independent, so this
# scales linearly with the number of cores.
#
# Usage:
#   bash 03a_mosdepth_coverage_summary.sh
#
#   # Or submit as a SLURM job (benefits from many cores):
#   sbatch --cpus-per-task=$(nproc) --mem=8G --time=00:30:00 \
#     --wrap="bash 03a_mosdepth_coverage_summary.sh"
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SLURM_SUBMIT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SCRIPT_DIR}/config.sh"
elif [[ -f "$(pwd)/config.sh" ]]; then
  CONFIG_FILE="$(pwd)/config.sh"
else
  echo "ERROR: Could not find config.sh (checked \$SLURM_SUBMIT_DIR, script directory, and current directory)."
  exit 1
fi
source "${CONFIG_FILE}"

mkdir -p "${QC_OUTPUT}"

# ── Auto-detect mosdepth directory ──────────────────────────────────────────
if [[ ! -d "${MOSDEPTH_DIR}" ]]; then
  echo "ERROR: Mosdepth output directory not found: ${MOSDEPTH_DIR}"
  echo "  Set MOSDEPTH_DIR or ensure mosdepth was run in step 01."
  exit 1
fi

REGIONS_FILES=( "${MOSDEPTH_DIR}"/*.regions.bed.gz )
if [[ ! -f "${REGIONS_FILES[0]}" ]]; then
  echo "ERROR: No *.regions.bed.gz files found in ${MOSDEPTH_DIR}"
  exit 1
fi

NSAMPLES=${#REGIONS_FILES[@]}
NCPUS="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 4)}"
OUT_TSV="${QC_OUTPUT}/mosdepth_coverage_summary.tsv"

# ── Check for HQ (high-quality) stats capability ───────────────────────────
# HQ stats use only autosomal bins that do NOT overlap the exclusion BED,
# providing a cleaner baseline for downstream metrics like mtDNA CN.
HQ_ENABLED=0
if command -v bedtools &>/dev/null; then
  if [[ -n "${BED_EXCLUDE:-}" && -f "${BED_EXCLUDE}" ]]; then
    HQ_ENABLED=1
  fi
fi
export HQ_ENABLED BED_EXCLUDE

# ── Progress tracking ───────────────────────────────────────────────────────
PROGRESS_DIR=$(mktemp -d)
export PROGRESS_DIR NSAMPLES
trap 'rm -rf "${PROGRESS_DIR}"' EXIT

echo "============================================================"
echo " Computing autosomal coverage summary statistics"
echo "============================================================"
echo ""
echo " Mosdepth dir:  ${MOSDEPTH_DIR}"
echo " Samples:       ${NSAMPLES}"
echo " Parallel jobs: ${NCPUS}"
echo " Output:        ${OUT_TSV}"
if [[ "${HQ_ENABLED}" -eq 1 ]]; then
  echo " Exclude BED:   ${BED_EXCLUDE}"
  echo " HQ stats:      enabled (bins not overlapping exclude BED)"
else
  echo " HQ stats:      disabled (requires bedtools + BED_EXCLUDE)"
fi
echo ""

# ── Helper: compute summary stats from a sorted coverage file ───────────────
# Input:  path to a file with one coverage value per line, numerically sorted.
# Output: space-separated: mean median sd mad iqr
#         Returns "NA NA NA NA NA" when the file is empty.
compute_stats_for_sorted() {
  local sorted_file="$1"
  local n
  n=$(wc -l < "${sorted_file}")

  if [[ "${n}" -eq 0 ]]; then
    echo "NA NA NA NA NA"
    return
  fi

  # Compute quantile positions (1-based indexing on sorted data)
  local mp1 mp2 qp1 qp3
  if (( n % 2 == 1 )); then
    mp1=$(( n / 2 + 1 )); mp2=${mp1}
  else
    mp1=$(( n / 2 )); mp2=$(( n / 2 + 1 ))
  fi
  qp1=$(( n / 4 + 1 ))
  qp3=$(( 3 * n / 4 + 1 ))

  # Single pass over sorted data → mean, SD, median, Q1, Q3
  local stats
  stats=$(awk -v n="${n}" -v mp1="${mp1}" -v mp2="${mp2}" \
              -v qp1="${qp1}" -v qp3="${qp3}" '
    { sum += $1; sumsq += $1 * $1 }
    NR == mp1 { m1 = $1 }
    NR == mp2 { m2 = $1 }
    NR == qp1 { q1v = $1 }
    NR == qp3 { q3v = $1 }
    END {
      mean = sum / n
      var = (n > 1) ? (sumsq - sum*sum/n) / (n - 1) : 0
      sd = sqrt(var > 0 ? var : 0)
      median = (mp1 == mp2) ? m1 : (m1 + m2) / 2
      printf "%.4f %.4f %.4f %.4f %.4f\n", mean, median, sd, q1v, q3v
    }
  ' "${sorted_file}")

  local mean median sd q1 q3
  read -r mean median sd q1 q3 <<< "${stats}"
  local iqr
  iqr=$(awk "BEGIN { printf \"%.4f\", ${q3} - ${q1} }")

  # MAD = median of |xi - median|
  local mad
  mad=$(awk -v med="${median}" '{ d = $1 - med; print (d < 0) ? -d : d }' "${sorted_file}" \
    | sort -g \
    | awk -v n="${n}" -v mp1="${mp1}" -v mp2="${mp2}" '
        NR == mp1 { m1 = $1 }
        NR == mp2 { m2 = $1; exit }
        END { printf "%.4f\n", (mp1 == mp2) ? m1 : (m1 + m2) / 2 }')

  echo "${mean} ${median} ${sd} ${mad} ${iqr}"
}
export -f compute_stats_for_sorted

# ── Per-sample stats function ───────────────────────────────────────────────
# For each sample:
#   1. Extract autosomal (chr1–22) coverage values from regions.bed.gz
#   2. Sort them and compute mean, median, SD, MAD, IQR
#   3. (If HQ enabled) filter out bins overlapping the exclude BED, recompute
#
# This approach avoids storing millions of values in awk arrays.  The external
# sort handles memory management efficiently even for >3M bins.
compute_sample_stats() {
  local bed_gz="$1"
  local fname
  fname="$(basename "${bed_gz}")"
  local sample_id="${fname%%.*}"

  # Extract autosomal BED lines (chr1–22)
  local tmp_auto_bed
  tmp_auto_bed="$(mktemp)"
  zcat "${bed_gz}" \
    | awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/' \
    > "${tmp_auto_bed}"

  # ALL autosomal stats — sort coverage values and compute
  local tmp_sorted
  tmp_sorted="$(mktemp)"
  awk -F'\t' '{ print $4 }' "${tmp_auto_bed}" | sort -g > "${tmp_sorted}"

  local all_stats
  all_stats=$(compute_stats_for_sorted "${tmp_sorted}")
  rm -f "${tmp_sorted}"

  # HQ autosomal stats — exclude bins overlapping the exclusion BED
  local hq_stats="NA NA NA NA NA"
  if [[ "${HQ_ENABLED}" -eq 1 ]]; then
    local tmp_hq_sorted
    tmp_hq_sorted="$(mktemp)"
    bedtools intersect -v -a "${tmp_auto_bed}" -b "${BED_EXCLUDE}" \
      | awk -F'\t' '{ print $4 }' \
      | sort -g > "${tmp_hq_sorted}"
    hq_stats=$(compute_stats_for_sorted "${tmp_hq_sorted}")
    rm -f "${tmp_hq_sorted}"
  fi

  rm -f "${tmp_auto_bed}"

  # Format output: tab-separated fields
  local mean median sd mad iqr
  read -r mean median sd mad iqr <<< "${all_stats}"
  local hq_mean hq_median hq_sd hq_mad hq_iqr
  read -r hq_mean hq_median hq_sd hq_mad hq_iqr <<< "${hq_stats}"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample_id}" "${mean}" "${median}" "${sd}" "${mad}" "${iqr}" \
    "${hq_mean}" "${hq_median}" "${hq_sd}" "${hq_mad}" "${hq_iqr}"

  # Progress reporting
  touch "${PROGRESS_DIR}/${sample_id}"
  local ndone
  ndone=$(find "${PROGRESS_DIR}" -maxdepth 1 -type f | wc -l)
  printf "\r  [%d/%d] samples processed" "${ndone}" "${NSAMPLES}" >&2
}
export -f compute_sample_stats

# ── Write header ────────────────────────────────────────────────────────────
echo -e "SAMPLE_ID\tMEAN_COV\tMEDIAN_COV\tSD_COV\tMAD_COV\tIQR_COV\tHQ_MEAN_COV\tHQ_MEDIAN_COV\tHQ_SD_COV\tHQ_MAD_COV\tHQ_IQR_COV" > "${OUT_TSV}"

# ── Parallel execution across all samples ───────────────────────────────────
# xargs -P parallelizes across NCPUS cores.  Each sample is fully independent.
printf '%s\n' "${REGIONS_FILES[@]}" \
  | xargs -P "${NCPUS}" -I{} bash -c 'compute_sample_stats "$@"' _ {} \
  | sort -t$'\t' -k1,1 \
  >> "${OUT_TSV}"

echo "" >&2

NROWS=$(tail -n +2 "${OUT_TSV}" | wc -l)

echo ""
echo "============================================================"
echo " Coverage summary complete: ${OUT_TSV}"
echo " Samples processed: ${NROWS}"
echo ""
echo " Preview (first 5 rows):"
head -6 "${OUT_TSV}" | column -t -s$'\t'
echo ""
if [[ "${HQ_ENABLED}" -eq 1 ]]; then
  echo " Columns: SAMPLE_ID  MEAN_COV  MEDIAN_COV  SD_COV  MAD_COV  IQR_COV"
  echo "          HQ_MEAN_COV  HQ_MEDIAN_COV  HQ_SD_COV  HQ_MAD_COV  HQ_IQR_COV"
else
  echo " Columns: SAMPLE_ID  MEAN_COV  MEDIAN_COV  SD_COV  MAD_COV  IQR_COV"
  echo "          HQ_* columns are NA (bedtools or BED_EXCLUDE not available)"
fi
echo "============================================================"
