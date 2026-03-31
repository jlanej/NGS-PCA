#!/usr/bin/env bash
# =============================================================================
# fetch_alignment_indexes.sh — Download NYGC 30x alignment index files
# =============================================================================
#
# Downloads the 1000G NYGC 30x alignment index files (which list BAS file URLs
# for all 3,202 high-coverage samples) and saves them into example/1000G_highcov/data/
# so they can be committed to the repository.
#
# These files are used by 03_collect_qc.sh to discover and download per-sample
# BAS alignment statistics.
#
# Usage:
#   bash example/1000G_highcov/fetch_alignment_indexes.sh
#
# After running, commit the downloaded files:
#   git add example/1000G_highcov/data/*.alignment.index
#   git commit -m "Add NYGC 30x alignment index files for BAS lookup"
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
mkdir -p "${DATA_DIR}"

# ── Helper: validate alignment index content ──
validate_alignment_index() {
  local file="$1"
  [[ ! -s "${file}" ]] && return 1
  if head -5 "${file}" | grep -qiE '<html|<h1>|<!DOCTYPE'; then
    echo "  REJECTED: file contains HTML (likely an error page)"
    return 1
  fi
  local data_lines
  data_lines=$(grep -cvE '^#|^$' "${file}" 2>/dev/null || echo "0")
  if [[ "${data_lines}" -eq 0 ]]; then
    echo "  REJECTED: no data lines found"
    return 1
  fi
  if ! awk -F'\t' '/^[^#]/ && $5 ~ /\.bas/ { found=1; exit } END { exit !found }' "${file}" 2>/dev/null; then
    echo "  REJECTED: no .bas URLs in column 5"
    return 1
  fi
  echo "  VALID: ${data_lines} data lines with BAS URLs"
  return 0
}

# ── Helper: try downloading from multiple URLs ──
try_download() {
  local dest="$1"
  shift
  local urls=("$@")

  for url in "${urls[@]}"; do
    echo "  Trying: ${url}"
    if curl -sSL --fail --max-time 120 -o "${dest}" "${url}" 2>/dev/null && validate_alignment_index "${dest}"; then
      return 0
    elif wget -q --timeout=120 -O "${dest}" "${url}" 2>/dev/null && validate_alignment_index "${dest}"; then
      return 0
    else
      rm -f "${dest}"
    fi
  done
  return 1
}

echo "============================================================"
echo " Fetching NYGC 30x alignment index files"
echo "============================================================"
echo ""
echo " Target directory: ${DATA_DIR}"
echo ""

# ── 1. 2504 unrelated samples ──
echo "--- 2504 unrelated samples ---"
DEST_2504="${DATA_DIR}/1000G_2504_high_coverage.GRCh38DH.alignment.index"
if [[ -f "${DEST_2504}" ]] && validate_alignment_index "${DEST_2504}"; then
  echo "  Already present and valid: ${DEST_2504}"
else
  URLS_2504=(
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.GRCh38DH.alignment.index"
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1000G_2504_high_coverage.GRCh38DH.alignment.index"
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.GRCh38DH.alignment.index"
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1000G_2504_high_coverage.GRCh38DH.alignment.index"
    "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/1000G_2504_high_coverage.GRCh38DH.alignment.index"
    "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/alignment_indices/1000G_2504_high_coverage.GRCh38DH.alignment.index"
  )
  if try_download "${DEST_2504}" "${URLS_2504[@]}"; then
    echo "  SUCCESS: downloaded 2504-sample alignment index"
  else
    echo "  FAILED: could not download 2504-sample alignment index from any source"
    echo "  You may need to download it manually and place it at:"
    echo "    ${DEST_2504}"
  fi
fi
echo ""

# ── 2. 698 related samples ──
echo "--- 698 related samples ---"
DEST_698="${DATA_DIR}/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
if [[ -f "${DEST_698}" ]] && validate_alignment_index "${DEST_698}"; then
  echo "  Already present and valid: ${DEST_698}"
else
  URLS_698=(
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
    "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
    "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/alignment_indices/1000G_698_related_high_coverage.GRCh38DH.alignment.index"
  )
  if try_download "${DEST_698}" "${URLS_698[@]}"; then
    echo "  SUCCESS: downloaded 698-sample alignment index"
  else
    echo "  FAILED: could not download 698-sample alignment index from any source"
    echo "  You may need to download it manually and place it at:"
    echo "    ${DEST_698}"
  fi
fi
echo ""

# ── 3. IGSR pilot (24 samples, from GitHub — always accessible) ──
echo "--- IGSR pilot (24 samples, GitHub fallback) ---"
DEST_PILOT="${DATA_DIR}/1000genomes.high_coverage.GRCh38DH.alignment.index"
if [[ -f "${DEST_PILOT}" ]] && validate_alignment_index "${DEST_PILOT}"; then
  echo "  Already present and valid: ${DEST_PILOT}"
else
  PILOT_URL="https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.high_coverage.GRCh38DH.alignment.index"
  echo "  Trying: ${PILOT_URL}"
  if curl -sSL --fail --max-time 60 -o "${DEST_PILOT}" "${PILOT_URL}" 2>/dev/null && validate_alignment_index "${DEST_PILOT}"; then
    echo "  SUCCESS: downloaded IGSR pilot alignment index (24 samples)"
  else
    echo "  FAILED: could not download pilot alignment index"
    rm -f "${DEST_PILOT}"
  fi
fi
echo ""

# ── Summary ──
echo "============================================================"
echo " Summary"
echo "============================================================"
for f in "${DATA_DIR}"/*.alignment.index; do
  [[ -f "${f}" ]] || continue
  lines=$(grep -cvE '^#|^$' "${f}" 2>/dev/null || echo "0")
  echo "  $(basename "${f}"): ${lines} samples"
done
echo ""
echo " To commit these files to the repository:"
echo "   git add example/1000G_highcov/data/*.alignment.index"
echo "   git commit -m 'Add alignment index files for BAS lookup'"
echo "============================================================"
