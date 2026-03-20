#!/usr/bin/env bash
# =============================================================================
# 01_download_and_mosdepth.sh — Download CRAM, run mosdepth, clean up
# =============================================================================
#
# SLURM array job: each task processes one sample from the manifest.
#   1. Download the CRAM + CRAI from EBI via Aspera (falls back to wget)
#   2. Run mosdepth to compute 1 kb bin coverage
#   3. Download BAS file (per-sample QC stats)
#   4. Remove the downloaded CRAM/CRAI to free disk space
#
# Usage:
#   # After running 00_setup.sh:
#   sbatch 01_download_and_mosdepth.sh
#
#   # Or process a subset (e.g. first 100 samples):
#   sbatch --array=1-100 01_download_and_mosdepth.sh
#
#   # Retry specific failed tasks:
#   sbatch --array=42,99,256 01_download_and_mosdepth.sh
# =============================================================================

#SBATCH --job-name=1kG_mosdepth
#SBATCH --output=logs/mosdepth_%A_%a.out
#SBATCH --error=logs/mosdepth_%A_%a.err
#SBATCH --array=1-3202
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00

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

# Attempt to load the Aspera module if the HPC module system is available.
# This is a no-op when 'module' is not present or 'aspera' is not a module.
module load aspera 2>/dev/null || true

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit this script as a SLURM array job."
  exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest missing or empty: ${MANIFEST}"
  echo "Run setup first: bash 00_setup.sh"
  exit 1
fi

TOTAL_SAMPLES=$(tail -n +2 "${MANIFEST}" | wc -l)
if (( TOTAL_SAMPLES < MIN_MANIFEST_SAMPLES )); then
  echo "ERROR: Manifest has only ${TOTAL_SAMPLES} samples (minimum expected: ${MIN_MANIFEST_SAMPLES})."
  echo "This likely indicates an incomplete manifest/index download."
  echo "Rebuild it with: rm -f \"${MANIFEST}\" \"${INDEX_FILE}\" && bash 00_setup.sh"
  exit 1
fi
if (( TOTAL_SAMPLES != EXPECTED_MANIFEST_SAMPLES )); then
  echo "WARNING: Manifest has ${TOTAL_SAMPLES} samples (expected ${EXPECTED_MANIFEST_SAMPLES})."
fi
if (( SLURM_ARRAY_TASK_ID > TOTAL_SAMPLES )); then
  echo "ERROR: Array task ${SLURM_ARRAY_TASK_ID} exceeds manifest size (${TOTAL_SAMPLES} samples)."
  echo "Submit with a bounded array, e.g.: sbatch --array=1-${TOTAL_SAMPLES} 01_download_and_mosdepth.sh"
  exit 1
fi

# ── Resolve the sample for this array task ───────────────────────────────────
# SLURM_ARRAY_TASK_ID is 1-based; manifest line 1 is the header.
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
LINE=$(sed -n "${LINE_NUM}p" "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
  echo "ERROR: No manifest entry for task ${SLURM_ARRAY_TASK_ID} (line ${LINE_NUM})"
  exit 1
fi

# Parse tab-separated manifest columns in one pass (preserves empty trailing BAS field).
IFS=$'\t' read -r SAMPLE_ID CRAM_FTP_URL CRAI_FTP_URL CRAM_MD5 BAS_FTP_URL _ <<< "${LINE}"

if [[ -z "${SAMPLE_ID}" || -z "${CRAM_FTP_URL}" || -z "${CRAI_FTP_URL}" ]]; then
  echo "ERROR: Manifest entry at line ${LINE_NUM} is missing required columns."
  exit 1
fi
if [[ "${CRAM_FTP_URL}" != "${EXPECTED_FTP_PREFIX}"* || "${CRAI_FTP_URL}" != "${EXPECTED_FTP_PREFIX}"* ]]; then
  echo "ERROR: Manifest entry at line ${LINE_NUM} has unsupported CRAM/CRAI source."
  echo "  CRAM: ${CRAM_FTP_URL}"
  echo "  CRAI: ${CRAI_FTP_URL}"
  exit 1
fi

echo "============================================================"
echo " Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID}"
echo " CRAM: ${CRAM_FTP_URL}"
echo " Started: $(date)"
echo "============================================================"

# ── Skip if mosdepth output already exists ───────────────────────────────────
MOSDEPTH_OUTPUT="${MOSDEPTH_DIR}/${SAMPLE_ID}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
LOCAL_CRAM="${CRAM_DIR}/${SAMPLE_ID}.cram"
LOCAL_CRAI="${CRAM_DIR}/${SAMPLE_ID}.cram.crai"
if [[ -f "${MOSDEPTH_OUTPUT}" ]]; then
  echo "SKIP: mosdepth output already exists: ${MOSDEPTH_OUTPUT}"
  # Clean up any leftover CRAM/CRAI from a previous partial run to free disk space
  rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
  exit 0
fi

# ── Stage 1: Download CRAM + CRAI ───────────────────────────────────────────

download_aspera() {
  # Convert FTP URL to Aspera path and download using system ascp.
  #   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/... → vol1/ftp/...
  # ascp saves with the original remote filename; caller renames as needed.
  local ftp_url="$1"
  local dest="$2"
  local aspera_path
  # Strip the FTP server prefix to get the Aspera-compatible remote path
  # e.g. "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/..." → "vol1/ftp/..."
  aspera_path="${ftp_url//ftp:\/\/ftp.1000genomes.ebi.ac.uk\//}"
  # Use system ascp (install Aspera Connect on your HPC or load it as a module)
  command -v ascp &>/dev/null || return 1
  local remote_basename
  remote_basename="$(basename "${aspera_path}")"
  local downloaded="${CRAM_DIR}/${remote_basename}"
  if ! ascp -i "${ASPERA_SSH_KEY}" \
    -Tr -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" -L- \
    "${EBI_ASPERA_USER}:${aspera_path}" \
    "${CRAM_DIR}/"; then
    # Clean up any partial download left by ascp
    rm -f "${downloaded}" "${dest}"
    return 1
  fi
  # Rename to the desired local filename if it differs from the remote name
  if [[ "${downloaded}" != "${dest}" && -f "${downloaded}" ]]; then
    mv "${downloaded}" "${dest}"
  fi
}

download_wget() {
  local ftp_url="$1"
  local dest="$2"
  if ! wget -q -O "${dest}" "${ftp_url}"; then
    rm -f "${dest}"
    return 1
  fi
}

if [[ ! -f "${LOCAL_CRAM}" ]]; then
  echo "[1/4] Downloading CRAM..."
  if download_aspera "${CRAM_FTP_URL}" "${LOCAL_CRAM}"; then
    echo "  Aspera download complete."
  else
    echo "  Aspera failed; falling back to wget..."
    download_wget "${CRAM_FTP_URL}" "${LOCAL_CRAM}"
    echo "  wget download complete."
  fi
else
  echo "[1/4] CRAM already present: ${LOCAL_CRAM}"
fi

if [[ ! -f "${LOCAL_CRAI}" ]]; then
  echo "  Downloading CRAI..."
  if download_aspera "${CRAI_FTP_URL}" "${LOCAL_CRAI}"; then
    echo "  Aspera download complete."
  else
    echo "  Aspera failed; falling back to wget..."
    download_wget "${CRAI_FTP_URL}" "${LOCAL_CRAI}"
    echo "  wget download complete."
  fi
else
  echo "  CRAI already present: ${LOCAL_CRAI}"
fi

# Optional: verify CRAM MD5
if [[ -n "${CRAM_MD5}" ]]; then
  echo "  Verifying CRAM MD5..."
  ACTUAL_MD5=$(md5sum "${LOCAL_CRAM}" | awk '{print $1}')
  if [[ "${ACTUAL_MD5}" != "${CRAM_MD5}" ]]; then
    echo "  WARNING: MD5 mismatch (expected: ${CRAM_MD5}, got: ${ACTUAL_MD5})"
    echo "  Removing corrupt file and exiting. Re-submit this task to retry."
    rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
    exit 1
  fi
  echo "  MD5 verified."
fi

# ── Stage 2: Run mosdepth ───────────────────────────────────────────────────
echo "[2/4] Running mosdepth (bin size: ${MOSDEPTH_BIN_SIZE} bp, threads: ${MOSDEPTH_THREADS})..."

apptainer exec \
  --bind "${CRAM_DIR}":/crams \
  --bind "${MOSDEPTH_DIR}":/mosdepth \
  --bind "${REF_DIR}":/ref \
  "${SIF_IMAGE}" \
  mosdepth \
    -n \
    -t "${MOSDEPTH_THREADS}" \
    --by "${MOSDEPTH_BIN_SIZE}" \
    --fasta "/ref/$(basename "${REF_FASTA}")" \
    "/mosdepth/${SAMPLE_ID}.by${MOSDEPTH_BIN_SIZE}" \
    "/crams/${SAMPLE_ID}.cram"

# Verify output was created
if [[ ! -f "${MOSDEPTH_OUTPUT}" ]]; then
  echo "ERROR: mosdepth output not found: ${MOSDEPTH_OUTPUT}"
  exit 1
fi
echo "  mosdepth complete: ${MOSDEPTH_OUTPUT}"

# ── Stage 3: Download BAS file (per-sample QC stats, ~KB) ───────────────────
# The BAS file contains Picard-equivalent statistics: mapped reads, duplication
# rate, mean coverage, etc. It is retained permanently (unlike the CRAM).
mkdir -p "${BAS_DIR}"
LOCAL_BAS="${BAS_DIR}/${SAMPLE_ID}.bam.bas"

if [[ -f "${LOCAL_BAS}" ]]; then
  echo "[3/4] BAS file already present: ${LOCAL_BAS}"
elif [[ -z "${BAS_FTP_URL}" ]]; then
  echo "[3/4] No BAS URL in manifest; skipping BAS download."
else
  echo "[3/4] Downloading BAS file (per-sample QC)..."
  if download_aspera "${BAS_FTP_URL}" "${LOCAL_BAS}"; then
    echo "  Aspera download complete."
  else
    echo "  Aspera failed; falling back to wget..."
    download_wget "${BAS_FTP_URL}" "${LOCAL_BAS}"
    echo "  wget download complete."
  fi
fi

# ── Stage 4: Clean up CRAM/CRAI ─────────────────────────────────────────────
echo "[4/4] Cleaning up downloaded CRAM/CRAI..."
rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
echo "  Removed: ${LOCAL_CRAM}"
echo "  Removed: ${LOCAL_CRAI}"

echo ""
echo "============================================================"
echo " Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE_ID}) complete."
echo " Finished: $(date)"
echo "============================================================"
