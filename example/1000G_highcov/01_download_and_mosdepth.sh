#!/usr/bin/env bash
# =============================================================================
# 01_download_and_mosdepth.sh — Download CRAM, run mosdepth, clean up
# =============================================================================
#
# SLURM array job: each task processes one sample from the manifest.
#   1. Download the CRAM + CRAI from EBI via Aspera (falls back to wget)
#   2. Run mosdepth to compute 1 kb bin coverage
#   3. Remove the downloaded CRAM/CRAI to free disk space
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
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

# ── Resolve the sample for this array task ───────────────────────────────────
# SLURM_ARRAY_TASK_ID is 1-based; manifest line 1 is the header.
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
LINE=$(sed -n "${LINE_NUM}p" "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
  echo "ERROR: No manifest entry for task ${SLURM_ARRAY_TASK_ID} (line ${LINE_NUM})"
  exit 1
fi

SAMPLE_ID=$(echo "${LINE}" | cut -f1)
CRAM_FTP_URL=$(echo "${LINE}" | cut -f2)
CRAI_FTP_URL=$(echo "${LINE}" | cut -f3)
CRAM_MD5=$(echo "${LINE}" | cut -f4)

echo "============================================================"
echo " Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID}"
echo " CRAM: ${CRAM_FTP_URL}"
echo " Started: $(date)"
echo "============================================================"

# ── Skip if mosdepth output already exists ───────────────────────────────────
MOSDEPTH_OUTPUT="${MOSDEPTH_DIR}/${SAMPLE_ID}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
if [[ -f "${MOSDEPTH_OUTPUT}" ]]; then
  echo "SKIP: mosdepth output already exists: ${MOSDEPTH_OUTPUT}"
  exit 0
fi

# ── Stage 1: Download CRAM + CRAI ───────────────────────────────────────────
LOCAL_CRAM="${CRAM_DIR}/${SAMPLE_ID}.cram"
LOCAL_CRAI="${CRAM_DIR}/${SAMPLE_ID}.cram.crai"

download_aspera() {
  # Convert FTP URL to Aspera path and download to CRAM_DIR.
  #   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/... → vol1/ftp/...
  # ascp saves with the original remote filename; caller renames as needed.
  local ftp_url="$1"
  local dest="$2"
  local aspera_path
  aspera_path="${ftp_url//ftp:\/\/ftp.1000genomes.ebi.ac.uk\//}"
  apptainer exec \
    --bind "${CRAM_DIR}":/download \
    "${SIF_IMAGE}" \
    ascp -i "${ASPERA_SSH_KEY}" \
      -Tr -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" -L- \
      "${EBI_ASPERA_USER}:${aspera_path}" \
      /download/
  # Rename to the desired local filename if it differs from the remote name
  local remote_basename
  remote_basename="$(basename "${aspera_path}")"
  local downloaded="${CRAM_DIR}/${remote_basename}"
  if [[ "${downloaded}" != "${dest}" && -f "${downloaded}" ]]; then
    mv "${downloaded}" "${dest}"
  fi
}

download_wget() {
  local ftp_url="$1"
  local dest="$2"
  wget -q -O "${dest}" "${ftp_url}"
}

if [[ ! -f "${LOCAL_CRAM}" ]]; then
  echo "[1/3] Downloading CRAM..."
  if download_aspera "${CRAM_FTP_URL}" "${LOCAL_CRAM}"; then
    echo "  Aspera download complete."
  else
    echo "  Aspera failed; falling back to wget..."
    download_wget "${CRAM_FTP_URL}" "${LOCAL_CRAM}"
    echo "  wget download complete."
  fi
else
  echo "[1/3] CRAM already present: ${LOCAL_CRAM}"
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
echo "[2/3] Running mosdepth (bin size: ${MOSDEPTH_BIN_SIZE} bp, threads: ${MOSDEPTH_THREADS})..."

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

# ── Stage 3: Clean up CRAM/CRAI ─────────────────────────────────────────────
echo "[3/3] Cleaning up downloaded CRAM/CRAI..."
rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
echo "  Removed: ${LOCAL_CRAM}"
echo "  Removed: ${LOCAL_CRAI}"

echo ""
echo "============================================================"
echo " Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE_ID}) complete."
echo " Finished: $(date)"
echo "============================================================"
