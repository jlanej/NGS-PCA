#!/usr/bin/env bash
# =============================================================================
# 01_download_and_mosdepth.sh — Download CRAM, run mosdepth, clean up
# =============================================================================
#
# SLURM array job: each task processes one or more samples from the manifest.
#   1. Download the CRAM + CRAI from ENA via Aspera (falls back to wget)
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
#SBATCH --array=1-3202%200
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00

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

# Returns the expected mosdepth output path for a given sample ID.
mosdepth_output_path() {
  echo "${MOSDEPTH_DIR}/${1}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  # ── Pre-submission mode ────────────────────────────────────────────────────
  # When run directly (not as a SLURM array task), scan the manifest to find
  # which samples still need mosdepth output, then submit only those tasks.
  if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Manifest missing or empty: ${MANIFEST}"
    echo "Run setup first: bash 00_setup.sh"
    exit 1
  fi

  echo "Pre-submission check: scanning manifest for samples that need processing..."
  NEEDED_TASKS=()
  PREV_TASK_ID=""
  SAMPLE_IDX=0
  while IFS=$'\t' read -r SAMPLE_ID _; do
    SAMPLE_IDX=$(( SAMPLE_IDX + 1 ))
    MOSDEPTH_OUTPUT="$(mosdepth_output_path "${SAMPLE_ID}")"
    if [[ ! -s "${MOSDEPTH_OUTPUT}" ]]; then
      TASK_ID=$(( (SAMPLE_IDX - 1) / SAMPLES_PER_TASK + 1 ))
      if [[ "${TASK_ID}" != "${PREV_TASK_ID}" ]]; then
        NEEDED_TASKS+=("${TASK_ID}")
        PREV_TASK_ID="${TASK_ID}"
      fi
    fi
  done < <(tail -n +2 "${MANIFEST}")

  TOTAL_SAMPLES=${SAMPLE_IDX}
  TOTAL_TASKS=$(( (TOTAL_SAMPLES + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
  NEEDED_COUNT=${#NEEDED_TASKS[@]}
  ALREADY_DONE=$(( TOTAL_TASKS - NEEDED_COUNT ))
  echo "  Total samples: ${TOTAL_SAMPLES} (${TOTAL_TASKS} tasks at SAMPLES_PER_TASK=${SAMPLES_PER_TASK})"
  echo "  Already done:  ${ALREADY_DONE} tasks"
  echo "  To submit:     ${NEEDED_COUNT} tasks"

  if [[ ${NEEDED_COUNT} -eq 0 ]]; then
    echo "All samples already have mosdepth output. Nothing to submit."
    exit 0
  fi

  # Build a compact SLURM array spec (e.g. "1-5,7,9-12") from the needed task IDs.
  ARRAY_SPEC=""
  RANGE_START=""
  PREV_ID=""
  for TASK_ID in "${NEEDED_TASKS[@]}"; do
    if [[ -z "${PREV_ID}" ]]; then
      RANGE_START="${TASK_ID}"
      PREV_ID="${TASK_ID}"
    elif (( TASK_ID == PREV_ID + 1 )); then
      PREV_ID="${TASK_ID}"
    else
      [[ -n "${ARRAY_SPEC}" ]] && ARRAY_SPEC+=","
      if [[ "${RANGE_START}" == "${PREV_ID}" ]]; then
        ARRAY_SPEC+="${RANGE_START}"
      else
        ARRAY_SPEC+="${RANGE_START}-${PREV_ID}"
      fi
      RANGE_START="${TASK_ID}"
      PREV_ID="${TASK_ID}"
    fi
  done
  # Flush the final range.
  [[ -n "${ARRAY_SPEC}" ]] && ARRAY_SPEC+=","
  if [[ "${RANGE_START}" == "${PREV_ID}" ]]; then
    ARRAY_SPEC+="${RANGE_START}"
  else
    ARRAY_SPEC+="${RANGE_START}-${PREV_ID}"
  fi

  echo "Submitting: sbatch --array=${ARRAY_SPEC}%${MAX_CONCURRENT_TASKS} $(basename "${BASH_SOURCE[0]}")"
  exec sbatch --array="${ARRAY_SPEC}%${MAX_CONCURRENT_TASKS}" "${BASH_SOURCE[0]}"
fi

if (( SAMPLES_PER_TASK < 1 )); then
  echo "ERROR: SAMPLES_PER_TASK must be >= 1 (got: ${SAMPLES_PER_TASK})."
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
  echo "Rebuild it with: rm -f \"${MANIFEST}\" \"${INDEX_FILE_2504}\" \"${INDEX_FILE_698}\" && bash 00_setup.sh"
  exit 1
fi
if (( TOTAL_SAMPLES != EXPECTED_MANIFEST_SAMPLES )); then
  echo "WARNING: Manifest has ${TOTAL_SAMPLES} samples (expected ${EXPECTED_MANIFEST_SAMPLES})."
fi
TOTAL_TASKS=$(( (TOTAL_SAMPLES + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
if (( SLURM_ARRAY_TASK_ID > TOTAL_TASKS )); then
  echo "ERROR: Array task ${SLURM_ARRAY_TASK_ID} exceeds task count (${TOTAL_TASKS}) for ${TOTAL_SAMPLES} samples (SAMPLES_PER_TASK=${SAMPLES_PER_TASK})."
  echo "Submit with: sbatch --array=1-${TOTAL_TASKS}%${MAX_CONCURRENT_TASKS} 01_download_and_mosdepth.sh"
  exit 1
fi

TASK_START=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
TASK_END=$(( TASK_START + SAMPLES_PER_TASK - 1 ))
if (( TASK_END > TOTAL_SAMPLES )); then
  TASK_END=${TOTAL_SAMPLES}
fi

# ── Stage 1: Download CRAM + CRAI ───────────────────────────────────────────

download_aspera() {
  # Download a file via Aspera. Supports both EBI 1000G FTP and ENA SRA FTP URLs.
  #   ftp://ftp.sra.ebi.ac.uk/vol1/run/...         → era-fasp@fasp.sra.ebi.ac.uk:vol1/run/...
  #   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/...  → fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/...
  local ftp_url="$1"
  local dest="$2"
  local aspera_path aspera_user

  command -v ascp &>/dev/null || return 1

  if [[ "${ftp_url}" == *"ftp.sra.ebi.ac.uk"* ]]; then
    aspera_path="${ftp_url//ftp:\/\/ftp.sra.ebi.ac.uk\//}"
    aspera_user="${ENA_ASPERA_USER}"
  elif [[ "${ftp_url}" == *"ftp.1000genomes.ebi.ac.uk"* ]]; then
    aspera_path="${ftp_url//ftp:\/\/ftp.1000genomes.ebi.ac.uk\//}"
    aspera_user="${EBI_ASPERA_USER}"
  else
    return 1
  fi

  local remote_basename
  remote_basename="$(basename "${aspera_path}")"
  local downloaded="${CRAM_DIR}/${remote_basename}"
  if ! ascp -i "${ASPERA_SSH_KEY}" \
    -Tr -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" -L- \
    "${aspera_user}:${aspera_path}" \
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

process_manifest_line() {
  local line_num="$1"
  local line
  line=$(sed -n "${line_num}p" "${MANIFEST}")
  if [[ -z "${line}" ]]; then
    echo "ERROR: No manifest entry for line ${line_num}"
    return 1
  fi

  local SAMPLE_ID CRAM_FTP_URL CRAI_FTP_URL CRAM_MD5
  # Parse tab-separated manifest columns.
  IFS=$'\t' read -r SAMPLE_ID CRAM_FTP_URL CRAI_FTP_URL CRAM_MD5 _ <<< "${line}"
  if [[ -z "${SAMPLE_ID}" || -z "${CRAM_FTP_URL}" || -z "${CRAI_FTP_URL}" ]]; then
    echo "ERROR: Manifest entry at line ${line_num} is missing required columns."
    return 1
  fi
  if [[ "${CRAM_FTP_URL}" != "${EXPECTED_FTP_PREFIX}"* || "${CRAI_FTP_URL}" != "${EXPECTED_FTP_PREFIX}"* ]]; then
    echo "ERROR: Manifest entry at line ${line_num} has unsupported CRAM/CRAI source."
    echo "  CRAM: ${CRAM_FTP_URL}"
    echo "  CRAI: ${CRAI_FTP_URL}"
    return 1
  fi

  echo "============================================================"
  echo " Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID} (manifest line ${line_num})"
  echo " CRAM: ${CRAM_FTP_URL}"
  echo " Started: $(date)"
  echo "============================================================"

  # ── Skip if mosdepth output already exists and is non-empty ──────────────────
  MOSDEPTH_OUTPUT="$(mosdepth_output_path "${SAMPLE_ID}")"
  LOCAL_CRAM="${CRAM_DIR}/${SAMPLE_ID}.cram"
  LOCAL_CRAI="${CRAM_DIR}/${SAMPLE_ID}.cram.crai"
  if [[ -s "${MOSDEPTH_OUTPUT}" ]]; then
    echo "SKIP: mosdepth output already exists: ${MOSDEPTH_OUTPUT}"
    rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
    return 0
  fi

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
      return 1
    fi
    echo "  MD5 verified."
  fi

  # ── Stage 2: Run mosdepth ─────────────────────────────────────────────────
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
    return 1
  fi
  echo "  mosdepth complete: ${MOSDEPTH_OUTPUT}"

  # ── Stage 3: Clean up CRAM/CRAI ───────────────────────────────────────────
  echo "[3/3] Cleaning up downloaded CRAM/CRAI..."
  rm -f "${LOCAL_CRAM}" "${LOCAL_CRAI}"
  echo "  Removed: ${LOCAL_CRAM}"
  echo "  Removed: ${LOCAL_CRAI}"

  echo ""
  echo "============================================================"
  echo " Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE_ID}) complete."
  echo " Finished: $(date)"
  echo "============================================================"
}

echo "Task ${SLURM_ARRAY_TASK_ID} will process manifest sample indices ${TASK_START}-${TASK_END} of ${TOTAL_SAMPLES} (SAMPLES_PER_TASK=${SAMPLES_PER_TASK})."
FAILED_LINES=()
for SAMPLE_IDX in $(seq "${TASK_START}" "${TASK_END}"); do
  LINE_NUM=$(( SAMPLE_IDX + 1 ))
  if ! process_manifest_line "${LINE_NUM}"; then
    FAILED_LINES+=("${LINE_NUM}")
  fi
done

if (( ${#FAILED_LINES[@]} > 0 )); then
  echo "ERROR: ${#FAILED_LINES[@]} manifest line(s) failed in task ${SLURM_ARRAY_TASK_ID}: ${FAILED_LINES[*]}"
  exit 1
fi
