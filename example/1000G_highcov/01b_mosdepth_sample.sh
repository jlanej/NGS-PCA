#!/usr/bin/env bash
# =============================================================================
# 01b_mosdepth_sample.sh — mosdepth for one sample whose CRAM is already local
# =============================================================================
#
# Spawned by 01_download_and_mosdepth.sh the moment a sample's download passes
# its MD5 check - one job per sample, so mosdepth concurrency is whatever the
# scheduler grants and never waits on a download slot. Runs mosdepth (twice,
# with timing, when COMPARE_FAST_MODE=1), verifies the output, and removes the
# CRAM/CRAI.
#
# Usage (normally via the manager, but safe by hand for one sample):
#   sbatch 01b_mosdepth_sample.sh SAMPLE_ID MANIFEST_LINE_NUM
# =============================================================================

#SBATCH --job-name=1kG_mosdepth
#SBATCH --output=logs/mosdepth_%j.out
#SBATCH --error=logs/mosdepth_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00

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

SAMPLE_ID="${1:?usage: 01b_mosdepth_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]}"
LINE_NUM="${2:?usage: 01b_mosdepth_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]}"
# Optional: when given, the CRAM is verified here before mosdepth. The
# download manager verifies before submitting and omits it; the stage watcher
# (01c) passes it, because it dispatches on arrival without reading the file -
# verification then runs on the compute node, in parallel across jobs.
CRAM_MD5="${3:-}"

MOSDEPTH_OUTPUT="${MOSDEPTH_DIR}/${SAMPLE_ID}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
MOSDEPTH_FAST_OUTPUT="${MOSDEPTH_FAST_DIR}/${SAMPLE_ID}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
LOCAL_CRAM="${CRAM_DIR}/${SAMPLE_ID}.cram"
LOCAL_CRAI="${CRAM_DIR}/${SAMPLE_ID}.cram.crai"

# The manager creates this before submitting, but a watcher-dispatched job on
# a fresh tree must not hand apptainer a bind source that does not exist.
mkdir -p "${MOSDEPTH_DIR}"

echo "============================================================"
echo " mosdepth: ${SAMPLE_ID} (manifest line ${LINE_NUM})"
echo " Started: $(date)"
echo "============================================================"

# Idempotent under a requeue: if the work is already done, only tidy up.
needs_processing() {
  if [[ "${COMPARE_FAST_MODE}" == "1" ]]; then
    [[ ! -s "${MOSDEPTH_OUTPUT}" || ! -s "${MOSDEPTH_FAST_OUTPUT}" ]]
  else
    [[ ! -s "${MOSDEPTH_OUTPUT}" ]]
  fi
}
if ! needs_processing; then
  echo "SKIP: mosdepth output already exists: ${MOSDEPTH_OUTPUT}"
  rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
  exit 0
fi

if [[ ! -s "${LOCAL_CRAM}" || ! -s "${LOCAL_CRAI}" ]]; then
  echo "ERROR: CRAM or CRAI missing for ${SAMPLE_ID} - the download manager verifies both before"
  echo "       submitting this job, so something removed them. Re-run 01_download_and_mosdepth.sh."
  exit 1
fi

if [[ -n "${CRAM_MD5}" ]]; then
  echo "Verifying staged CRAM MD5..."
  ACTUAL_MD5=$(md5sum "${LOCAL_CRAM}" | awk '{print $1}')
  if [[ "${ACTUAL_MD5}" != "${CRAM_MD5}" ]]; then
    echo "ERROR: staged CRAM failed MD5 (expected: ${CRAM_MD5}, got: ${ACTUAL_MD5})."
    echo "       The transfer may still be in flight - the file is left in place for the"
    echo "       watcher to re-dispatch once it settles, or for a later manager sweep."
    exit 1
  fi
  echo "  MD5 verified."
fi

# run_mosdepth <sample> <output_dir> [extra mosdepth flags...]
run_mosdepth() {
  local sample="$1" out_dir="$2"
  shift 2
  apptainer exec \
    --bind "${CRAM_DIR}":/crams \
    --bind "${out_dir}":/mosdepth \
    --bind "${REF_DIR}":/ref \
    "${SIF_IMAGE}" \
    mosdepth \
      -n \
      -t "${MOSDEPTH_THREADS}" \
      --by "${MOSDEPTH_BIN_SIZE}" \
      --fasta "/ref/$(basename "${REF_FASTA}")" \
      "$@" \
      "/mosdepth/${sample}.by${MOSDEPTH_BIN_SIZE}" \
      "/crams/${sample}.cram"
}

# timed_mosdepth <sample> <mode> <first_mode> <output_dir> [extra flags...]
# Wall-clock seconds of the run, one headerless row per sample and mode so
# concurrent jobs never share a file. Columns:
#   sample  mode  first_mode  wall_s  mosdepth_version  threads  bin_size  host  epoch
timed_mosdepth() {
  local sample="$1" mode="$2" first_mode="$3" out_dir="$4"
  shift 4
  local t0 t1
  t0=$(date +%s)
  run_mosdepth "${sample}" "${out_dir}" "$@"
  t1=$(date +%s)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample}" "${mode}" "${first_mode}" "$(( t1 - t0 ))" \
    "${MOSDEPTH_VERSION_STR}" "${MOSDEPTH_THREADS}" "${MOSDEPTH_BIN_SIZE}" \
    "$(hostname)" "${t0}" \
    > "${MOSDEPTH_TIMING_DIR}/${sample}.${mode}.tsv"
  echo "  ${mode}: $(( t1 - t0 )) s"
}

if [[ "${COMPARE_FAST_MODE}" == "1" ]]; then
  mkdir -p "${MOSDEPTH_FAST_DIR}" "${MOSDEPTH_TIMING_DIR}"
  # recorded in every timing row, so the eval can say what was measured
  MOSDEPTH_VERSION_STR="$(apptainer exec "${SIF_IMAGE}" mosdepth --version 2>/dev/null | head -1 | awk '{print $NF}')"
  MOSDEPTH_VERSION_STR="${MOSDEPTH_VERSION_STR:-unknown}"

  # Both runs read the same already-verified CRAM, so download time is in
  # neither measurement. The order alternates by manifest line: whichever mode
  # runs second reads a warmer page cache, and alternating cancels that in
  # aggregate. first_mode is recorded so the eval can check it did.
  first_mode=normal
  (( LINE_NUM % 2 == 0 )) && first_mode=fast
  echo "[1/2] Running mosdepth twice for the fast-mode comparison (${first_mode} first; bin size: ${MOSDEPTH_BIN_SIZE} bp, threads: ${MOSDEPTH_THREADS})..."
  if [[ "${first_mode}" == "normal" ]]; then
    timed_mosdepth "${SAMPLE_ID}" normal "${first_mode}" "${MOSDEPTH_DIR}"
    timed_mosdepth "${SAMPLE_ID}" fast "${first_mode}" "${MOSDEPTH_FAST_DIR}" --fast-mode
  else
    timed_mosdepth "${SAMPLE_ID}" fast "${first_mode}" "${MOSDEPTH_FAST_DIR}" --fast-mode
    timed_mosdepth "${SAMPLE_ID}" normal "${first_mode}" "${MOSDEPTH_DIR}"
  fi
  if [[ ! -f "${MOSDEPTH_OUTPUT}" || ! -f "${MOSDEPTH_FAST_OUTPUT}" ]]; then
    echo "ERROR: mosdepth output not found: ${MOSDEPTH_OUTPUT} and/or ${MOSDEPTH_FAST_OUTPUT}"
    exit 1
  fi
  echo "  mosdepth complete: ${MOSDEPTH_OUTPUT}"
  echo "  mosdepth complete: ${MOSDEPTH_FAST_OUTPUT}"
else
  echo "[1/2] Running mosdepth (bin size: ${MOSDEPTH_BIN_SIZE} bp, threads: ${MOSDEPTH_THREADS})..."
  run_mosdepth "${SAMPLE_ID}" "${MOSDEPTH_DIR}"
  if [[ ! -f "${MOSDEPTH_OUTPUT}" ]]; then
    echo "ERROR: mosdepth output not found: ${MOSDEPTH_OUTPUT}"
    exit 1
  fi
  echo "  mosdepth complete: ${MOSDEPTH_OUTPUT}"
fi

echo "[2/2] Cleaning up downloaded CRAM/CRAI..."
rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
echo "  Removed: ${LOCAL_CRAM}"
echo "  Removed: ${LOCAL_CRAI}"

echo ""
echo "============================================================"
echo " ${SAMPLE_ID} complete."
echo " Finished: $(date)"
echo "============================================================"
