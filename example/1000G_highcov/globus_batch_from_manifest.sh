#!/usr/bin/env bash
# =============================================================================
# globus_batch_from_manifest.sh — emit Globus batch files for the cohort
# =============================================================================
#
# ENA exposes ftp.sra.ebi.ac.uk through the "EMBL-EBI Public Data" Globus
# collection under /vol1, so a manifest URL becomes a Globus source path by
# stripping the host. Each emitted line renames the file on arrival to the
# <SAMPLE>.cram / <SAMPLE>.cram.crai names the pipeline expects, quoted
# because the batch parser splits lines shlex-style.
#
# Two batch files are written, indexes and CRAMs apart, because Globus decides
# file order within a task itself - a single mixed task was observed running
# the CRAMs first, leaving no completed pair for hours. Submitted as two
# tasks, the CRAI task finishes in minutes, and from then on every CRAM that
# completes completes its pair, which is what lets 01c_dispatch_staged.sh
# start mosdepth while the bulk of the transfer is still running.
#
# Every sample still needing mosdepth output is emitted, whether or not files
# for it already sit in CRAM_DIR - always submit with --sync-level checksum,
# which skips already-complete files nearly for free and re-transfers partial
# or corrupt ones. That makes a re-run of this script plus a re-submission the
# whole restart story after a cancelled or failed transfer.
#
# Usage:
#   bash globus_batch_from_manifest.sh [output_prefix]     # default: globus_batch
#   globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
#     --batch globus_batch_crai.txt --sync-level checksum --label "1kG crai" \
#     --notify failed,inactive
#   globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
#     --batch globus_batch_cram.txt --sync-level checksum --label "1kG cram" \
#     --notify failed,inactive
#
# Find your site collection with: globus endpoint search "<your site>"
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

OUT_PREFIX="${1:-globus_batch}"
CRAI_BATCH="${OUT_PREFIX}_crai.txt"
CRAM_BATCH="${OUT_PREFIX}_cram.txt"

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
: > "${CRAI_BATCH}"
: > "${CRAM_BATCH}"
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
  if ! needs_processing "${SAMPLE_ID}"; then
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  if [[ "${CRAM_URL}" != ftp://ftp.sra.ebi.ac.uk/* ]]; then
    echo "WARNING: skipping ${SAMPLE_ID} - not an ENA URL: ${CRAM_URL}" >&2
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  printf '"%s" "%s"\n' "${CRAI_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE_ID}.cram.crai" >> "${CRAI_BATCH}"
  printf '"%s" "%s"\n' "${CRAM_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE_ID}.cram" >> "${CRAM_BATCH}"
  EMITTED=$(( EMITTED + 1 ))
done < <(tail -n +2 "${MANIFEST}")

echo "Emitted ${EMITTED} samples: ${CRAI_BATCH} (indexes) and ${CRAM_BATCH} (CRAMs); skipped ${SKIPPED} done or damaged." >&2
echo "Submit the CRAI batch first, both with --sync-level checksum (see the header of this script)." >&2
if (( EMITTED == 0 )); then
  echo "Nothing to transfer." >&2
fi
