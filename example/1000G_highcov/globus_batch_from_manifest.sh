#!/usr/bin/env bash
# =============================================================================
# globus_batch_from_manifest.sh — emit a Globus batch file for the cohort
# =============================================================================
#
# ENA exposes ftp.sra.ebi.ac.uk through the "EMBL-EBI Public Data" Globus
# collection under /vol1, so a manifest URL becomes a Globus source path by
# stripping the host. Each emitted line renames the file on arrival to the
# <SAMPLE>.cram / <SAMPLE>.cram.crai names the pipeline expects, so a finished
# transfer is exactly what 01_download_and_mosdepth.sh treats as "already
# present": it will verify every file's MD5 and go straight to mosdepth.
#
# Only samples still needing work are emitted - mosdepth output missing, and
# no CRAM already on local disk - so the batch shrinks as the cohort finishes
# and a re-generated file is always safe to submit.
#
# Usage:
#   bash globus_batch_from_manifest.sh > globus_batch.txt
#   globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
#     --batch globus_batch.txt --label "1kG highcov CRAMs" --notify failed,inactive
#
# Destination paths are absolute (under CRAM_DIR); find your site collection
# with: globus endpoint search "<your site>"
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SCRIPT_DIR}/config.sh"
elif [[ -f "$(pwd)/config.sh" ]]; then
  CONFIG_FILE="$(pwd)/config.sh"
else
  echo "ERROR: Could not find config.sh." >&2
  exit 1
fi
source "${CONFIG_FILE}"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest missing or empty: ${MANIFEST} - run 00_setup.sh first." >&2
  exit 1
fi

mosdepth_output_path() {
  echo "${MOSDEPTH_DIR}/${1}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
}
mosdepth_fast_output_path() {
  echo "${MOSDEPTH_FAST_DIR}/${1}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
}
needs_processing() {
  local sample="$1"
  if [[ "${COMPARE_FAST_MODE}" == "1" ]]; then
    [[ ! -s "$(mosdepth_output_path "${sample}")" || ! -s "$(mosdepth_fast_output_path "${sample}")" ]]
  else
    [[ ! -s "$(mosdepth_output_path "${sample}")" ]]
  fi
}

EMITTED=0
SKIPPED=0
LINE_NUM=1
while IFS=$'\t' read -r SAMPLE_ID CRAM_URL CRAI_URL _; do
  LINE_NUM=$(( LINE_NUM + 1 ))   # manifest line, counting its header
  [[ -z "${SAMPLE_ID}" ]] && continue
  # Whitespace inside a field - a stray space in a sample name, a carriage
  # return from a damaged row - would split a batch line into extra tokens
  # and break the pipeline's own file naming later. Surface it, don't ship it.
  if [[ "${SAMPLE_ID}" =~ [[:space:]] || "${CRAM_URL}" =~ [[:space:]] || "${CRAI_URL}" =~ [[:space:]] ]]; then
    echo "WARNING: skipping manifest line ${LINE_NUM}: whitespace inside a field (sample '${SAMPLE_ID}') - fix that manifest row" >&2
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  if ! needs_processing "${SAMPLE_ID}" || [[ -f "${CRAM_DIR}/${SAMPLE_ID}.cram" ]]; then
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  if [[ "${CRAM_URL}" != ftp://ftp.sra.ebi.ac.uk/* ]]; then
    echo "WARNING: skipping ${SAMPLE_ID} - not an ENA URL: ${CRAM_URL}" >&2
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  # Quoted because the batch parser splits lines shlex-style
  printf '"%s" "%s"\n' "${CRAM_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE_ID}.cram"
  printf '"%s" "%s"\n' "${CRAI_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE_ID}.cram.crai"
  EMITTED=$(( EMITTED + 1 ))
done < <(tail -n +2 "${MANIFEST}")

echo "Emitted ${EMITTED} samples ($(( EMITTED * 2 )) files) to transfer; skipped ${SKIPPED} already done or staged." >&2
if (( EMITTED == 0 )); then
  echo "Nothing to transfer." >&2
fi
