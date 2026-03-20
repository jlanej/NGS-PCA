#!/usr/bin/env bash
# =============================================================================
# 00_setup.sh — One-time setup for the 1000G high-coverage NGS-PCA pipeline
# =============================================================================
#
# This script:
#   1. Creates the working directory structure
#   2. Pulls the NGS-PCA container image (includes mosdepth + ascp)
#   3. Downloads the GRCh38 reference genome (required for CRAM decoding)
#   4. Downloads the official IGSR alignment index and builds a manifest
#
# Prerequisites:
#   - Apptainer (Singularity) available on the system
#   - Internet access (or pre-downloaded files)
#
# Usage:
#   # Edit WORK_DIR in config.sh (or export it), then:
#   bash 00_setup.sh
#
# See also:
#   https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
#   https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data
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

validate_manifest() {
  local manifest_file="$1"
  if [[ ! -s "${manifest_file}" ]]; then
    echo "ERROR: Manifest missing or empty: ${manifest_file}"
    exit 1
  fi

  local header
  header="$(head -n1 "${manifest_file}")"
  if [[ "${header}" != $'SAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\tCRAM_MD5\tBAS_FTP_URL' ]]; then
    echo "ERROR: Manifest header is invalid:"
    echo "  got:      ${header}"
    echo "  expected: SAMPLE_ID<TAB>CRAM_FTP_URL<TAB>CRAI_FTP_URL<TAB>CRAM_MD5<TAB>BAS_FTP_URL"
    exit 1
  fi

  local invalid
  invalid="$(awk -F'\t' '
    NR==1 {next}
    NF<5 || $1=="" || $2=="" || $3=="" || $2 !~ /^ftp:\/\/ftp\.1000genomes\.ebi\.ac\.uk\// || $3 !~ /^ftp:\/\/ftp\.1000genomes\.ebi\.ac\.uk\// {print NR; exit}
  ' "${manifest_file}")"
  if [[ -n "${invalid}" ]]; then
    echo "ERROR: Manifest has invalid required fields/CRAM sources at line ${invalid}: ${manifest_file}"
    echo "  Re-run setup to regenerate it:"
    echo "    rm -f \"${manifest_file}\" && bash 00_setup.sh"
    exit 1
  fi

  local nfound
  nfound=$(tail -n +2 "${manifest_file}" | wc -l)
  if (( nfound < MIN_MANIFEST_SAMPLES )); then
    echo "ERROR: Manifest has only ${nfound} samples (minimum expected: ${MIN_MANIFEST_SAMPLES})."
    echo "  This often means the index download or parse was incomplete."
    echo "  Re-run setup to regenerate it:"
    echo "    rm -f \"${manifest_file}\" \"${INDEX_FILE}\" && bash 00_setup.sh"
    exit 1
  fi
  if (( nfound != EXPECTED_MANIFEST_SAMPLES )); then
    echo "WARNING: Manifest has ${nfound} samples (expected ${EXPECTED_MANIFEST_SAMPLES})."
    echo "  Continuing, but confirm the selected index is up to date."
  fi
}

echo "============================================================"
echo " NGS-PCA 1000G High-Coverage Pipeline — Setup"
echo "============================================================"
echo ""
echo " WORK_DIR:      ${WORK_DIR}"
echo " SIF_IMAGE:     ${SIF_IMAGE}"
echo " REF_FASTA:     ${REF_FASTA}"
echo " INDEX_FILE:    ${INDEX_FILE}"
echo " MANIFEST:      ${MANIFEST}"
echo ""

# ── 1. Create directory structure ────────────────────────────────────────────
echo "[1/5] Creating directory structure..."
mkdir -p "${CRAM_DIR}" "${MOSDEPTH_DIR}" "${BAS_DIR}" "${QC_OUTPUT}" "${NGSPCA_OUTPUT}" "${LOG_DIR}" "${REF_DIR}"
echo "  Done."

# ── 2. Pull the NGS-PCA container image ─────────────────────────────────────
if [[ -f "${SIF_IMAGE}" ]]; then
  echo "[2/5] Container image already exists: ${SIF_IMAGE}  (skipping pull)"
else
  echo "[2/5] Pulling NGS-PCA container image (includes mosdepth)..."
  echo "  apptainer pull ${SIF_IMAGE} docker://ghcr.io/jlanej/ngs-pca:latest"
  apptainer pull "${SIF_IMAGE}" docker://ghcr.io/jlanej/ngs-pca:latest
  echo "  Done."
fi

# Verify bundled tools
echo "  Verifying bundled tools..."
apptainer exec "${SIF_IMAGE}" mosdepth --version || echo "  WARNING: mosdepth not found in image"
echo ""

# ── 3. Download reference genome ────────────────────────────────────────────
if [[ -f "${REF_FASTA}" ]]; then
  echo "[3/5] Reference genome already exists: ${REF_FASTA}  (skipping download)"
else
  echo "[3/5] Downloading GRCh38 reference genome (~3 GB)..."
  echo "  Source: ${REF_URL}"

  # Try Aspera first (fast), fall back to wget
  ASPERA_REF_PATH="vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  if command -v ascp &>/dev/null; then
    echo "  Using Aspera for high-speed download..."
    ascp -i "${ASPERA_SSH_KEY}" \
      -Tr -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" -L- \
      "${EBI_ASPERA_USER}:${ASPERA_REF_PATH}" \
      "${REF_DIR}/" || {
        echo "  Aspera failed; falling back to wget..."
        wget -O "${REF_FASTA}" "${REF_URL}"
      }
  else
    echo "  Using wget..."
    wget -O "${REF_FASTA}" "${REF_URL}"
  fi
  echo "  Done."
fi

# Index the reference if needed
if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "  NOTE: Reference FASTA index (.fai) not found."
  echo "  mosdepth will create it automatically on first use, or run:"
  echo "    samtools faidx ${REF_FASTA}"
fi

# ── 4. Download IGSR index and generate manifest ────────────────────────────
if [[ -f "${MANIFEST}" ]]; then
  echo "[4/5] Manifest already exists: ${MANIFEST}  (skipping generation)"
  validate_manifest "${MANIFEST}"
  NFOUND=$(tail -n +2 "${MANIFEST}" | wc -l)
  echo "  ${NFOUND} samples listed."
else
  echo "[4/5] Downloading official IGSR alignment index and building manifest..."
  echo "  Source: ${INDEX_URL}"

  # Download the index file
  if [[ ! -f "${INDEX_FILE}" ]]; then
    curl -sSL -o "${INDEX_FILE}" "${INDEX_URL}"
  fi

  # Parse the index into a tab-separated manifest:
  #   SAMPLE_ID  CRAM_FTP_URL  CRAI_FTP_URL  CRAM_MD5  BAS_FTP_URL
  #
  # The index has paths like:
  #   ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/.../{SAMPLE}.alt_bwamem_GRCh38DH.*.cram
  # Note the single slash after "ftp:" — we fix it to "ftp://".
  echo -e "SAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\tCRAM_MD5\tBAS_FTP_URL" > "${MANIFEST}"
  grep -v '^#' "${INDEX_FILE}" | while IFS=$'\t' read -r cram cram_md5 crai _ bas _; do
    [[ -z "${cram}" ]] && continue
    # The IGSR index stores URLs with a single slash: "ftp:/ftp.1000genomes..."
    # Fix to standard double-slash: "ftp://ftp.1000genomes..."
    cram_url="${cram/ftp:\/ftp./ftp:\/\/ftp.}"
    crai_url="${crai/ftp:\/ftp./ftp:\/\/ftp.}"
    bas_url="${bas/ftp:\/ftp./ftp:\/\/ftp.}"
    # Extract sample ID from the CRAM filename
    cram_basename="$(basename "${cram}")"
    sample_id="${cram_basename%%.*}"
    echo -e "${sample_id}\t${cram_url}\t${crai_url}\t${cram_md5}\t${bas_url}"
  done >> "${MANIFEST}"

  NFOUND=$(tail -n +2 "${MANIFEST}" | wc -l)
  validate_manifest "${MANIFEST}"
  echo "  Generated manifest with ${NFOUND} samples."
  echo ""
  echo "  Preview (first 5 lines):"
  head -6 "${MANIFEST}" | column -t -s$'\t'
fi

# ── 5. Download IGSR sample panel (population + sex metadata) ────────────────
if [[ -f "${PANEL_FILE}" ]]; then
  echo "[5/5] Sample panel already exists: ${PANEL_FILE}  (skipping download)"
else
  echo "[5/5] Downloading IGSR sample panel (population + sex metadata)..."
  echo "  Source: ${PANEL_URL}"
  wget -q -O "${PANEL_FILE}" "${PANEL_URL}" || \
    curl -sSL -o "${PANEL_FILE}" "${PANEL_URL}"
  NPANEL=$(tail -n +2 "${PANEL_FILE}" | wc -l)
  echo "  Downloaded panel with ${NPANEL} samples."
fi

echo ""
echo "============================================================"
echo " Setup complete."
echo ""
echo " Next steps:"
echo "   1. Review the manifest:  head ${MANIFEST}"
echo "   2. Submit the download + mosdepth array job:"
echo "      sbatch 01_download_and_mosdepth.sh"
echo "   3. After step 2 completes, run NGS-PCA:"
echo "      sbatch 02_run_ngspca.sh"
echo "   4. Collect per-sample QC for batch-effect overlay:"
echo "      bash 03_collect_qc.sh"
echo "============================================================"
