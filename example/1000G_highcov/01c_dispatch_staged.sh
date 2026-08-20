#!/usr/bin/env bash
# =============================================================================
# 01c_dispatch_staged.sh — dispatch mosdepth as staged CRAM pairs arrive
# =============================================================================
#
# Run this alongside a bulk transfer (Globus) writing into $CRAM_DIR: it polls
# for samples whose CRAM and CRAI have both arrived and settled - mtimes
# unchanged across one poll interval - and submits each one's mosdepth job
# (01b) the moment it has, passing the manifest MD5 so verification happens on
# the compute node, in parallel across jobs, rather than serially here.
#
# It NEVER deletes or downloads anything, so it cannot fight the transfer.
# A file dispatched while secretly still in flight fails its MD5 in the job,
# loudly, with the file left in place; the watcher sees the job die, waits for
# the pair to settle again, and re-dispatches - up to WATCH_DISPATCH_ATTEMPTS
# times, after which the sample is parked with a warning for the final sweep.
#
# Do not run the download manager (01) at the same time as the transfer or
# this watcher: its MD5-mismatch path deletes and re-downloads, which is
# exactly wrong for a file another mover is mid-write on. Sequence is:
# transfer + watcher together, then one `bash 01_download_and_mosdepth.sh`
# sweep for anything the transfer missed.
#
# Exits when every sample is done, in flight, or parked - or after
# WATCH_IDLE_EXIT seconds in which nothing in CRAM_DIR changed at all
# (0 = wait forever). Growing files count as activity: a transfer that
# delivers all the CRAMs before their CRAIs completes no pairs for hours,
# and must not read as a stall.
#
# Usage:
#   bash 01c_dispatch_staged.sh              # submits itself as one job
#   WATCHER_LOCAL=1 bash 01c_dispatch_staged.sh   # run in place (tmux advised)
# =============================================================================

#SBATCH --job-name=1kG_stagewatch
#SBATCH --output=logs/stagewatch_%j.out
#SBATCH --error=logs/stagewatch_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=96:00:00

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

# CRAM_DIR too: the bulk transfer needs its destination to exist.
mkdir -p "${LOG_DIR}" "${DOWNLOAD_STATE_DIR}" "${CRAM_DIR}"
WATCH_STATE="${DOWNLOAD_STATE_DIR}/watch"
mkdir -p "${WATCH_STATE}"
MOSDEPTH_JOB_SCRIPT="${CONFIG_DIR}/01b_mosdepth_sample.sh"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest missing or empty: ${MANIFEST} - run 00_setup.sh first."
  exit 1
fi

if [[ -z "${SLURM_JOB_ID:-}" && "${WATCHER_LOCAL:-0}" != "1" ]]; then
  echo "Submitting the stage watcher (poll every ${WATCH_INTERVAL} s)..."
  echo "  Watcher log: ${LOG_DIR}/stagewatch_<jobid>.out"
  exec sbatch \
    --output="${LOG_DIR}/stagewatch_%j.out" \
    --error="${LOG_DIR}/stagewatch_%j.err" \
    "${BASH_SOURCE[0]}"
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

mtime_of() {
  stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null
}

# submit_mosdepth <sample> <line_num> <cram_md5>
submit_mosdepth() {
  local sample="$1" line_num="$2" cram_md5="$3"
  local try jobid
  for try in 1 2 3; do
    if jobid=$(sbatch --parsable \
        --job-name="1kG_mosdepth_${sample}" \
        --output="${LOG_DIR}/mosdepth_${sample}_%j.out" \
        --error="${LOG_DIR}/mosdepth_${sample}_%j.err" \
        --cpus-per-task="${MOSDEPTH_THREADS}" \
        --mem="${MOSDEPTH_MEM}" \
        --time="${MOSDEPTH_TIME}" \
        "${MOSDEPTH_JOB_SCRIPT}" "${sample}" "${line_num}" "${cram_md5}"); then
      echo "[$(date '+%H:%M:%S')] dispatched ${sample} (job ${jobid})"
      return 0
    fi
    sleep 30
  done
  return 1
}

# Manifest into arrays once; the loop re-reads nothing but the filesystem.
SAMPLES=()
LINE_NUMS=()
MD5S=()
LINE_NUM=1
while IFS=$'\t' read -r SAMPLE_ID _ _ CRAM_MD5 _; do
  LINE_NUM=$(( LINE_NUM + 1 ))
  [[ -z "${SAMPLE_ID}" ]] && continue
  SAMPLES+=("${SAMPLE_ID}")
  LINE_NUMS+=("${LINE_NUM}")
  MD5S+=("${CRAM_MD5}")
done < <(tail -n +2 "${MANIFEST}")

echo "Stage watcher: ${#SAMPLES[@]} samples in the manifest; polling ${CRAM_DIR} every ${WATCH_INTERVAL} s."
LAST_ACTIVITY=$(date +%s)
IN_FLIGHT_FILE="${WATCH_STATE}/in_flight.$$"
# Liveness is any change in CRAM_DIR, not just completed pairs: a transfer
# can spend hours delivering CRAMs before their CRAIs, and growing files must
# hold the idle clock open even though nothing is dispatchable yet.
ACTIVITY_STAMP="${WATCH_STATE}/activity.stamp"
touch "${ACTIVITY_STAMP}"

while :; do
  if [[ -n "$(find "${CRAM_DIR}" -maxdepth 1 -newer "${ACTIVITY_STAMP}" -print 2>/dev/null | head -1)" ]]; then
    touch "${ACTIVITY_STAMP}"
    LAST_ACTIVITY=$(date +%s)
  fi

  : > "${IN_FLIGHT_FILE}"
  if command -v squeue &>/dev/null; then
    squeue -h -u "${USER:-$(id -un)}" -o "%j" 2>/dev/null \
      | sed -n 's/^1kG_mosdepth_//p' > "${IN_FLIGHT_FILE}" || true
  fi

  WAITING=0
  DONE=0
  PARKED=0
  IN_FLIGHT=0
  for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[${idx}]}"
    if ! needs_processing "${sample}"; then
      DONE=$(( DONE + 1 ))
      continue
    fi
    if [[ -f "${WATCH_STATE}/${sample}.parked" ]]; then
      PARKED=$(( PARKED + 1 ))
      continue
    fi
    if grep -qxF "${sample}" "${IN_FLIGHT_FILE}"; then
      IN_FLIGHT=$(( IN_FLIGHT + 1 ))
      continue
    fi
    # A dispatch marker with no queued job and no output means the job died
    # (a premature dispatch failing its MD5 does this by design); clear it and
    # let the sample earn another dispatch by settling again.
    rm -f "${WATCH_STATE}/${sample}.dispatched"

    cram="${CRAM_DIR}/${sample}.cram"
    crai="${CRAM_DIR}/${sample}.cram.crai"
    if [[ ! -s "${cram}" || ! -s "${crai}" ]]; then
      WAITING=$(( WAITING + 1 ))
      continue
    fi
    mtimes="$(mtime_of "${cram}") $(mtime_of "${crai}")"
    state_file="${WATCH_STATE}/${sample}.mtimes"
    prev="$(cat "${state_file}" 2>/dev/null || true)"
    echo "${mtimes}" > "${state_file}"
    if [[ "${prev}" != "${mtimes}" ]]; then
      # first sighting, or still being written
      WAITING=$(( WAITING + 1 ))
      continue
    fi

    # An attempt counts toward parking only when the pair is unchanged since
    # the last failed dispatch: a file still filling in fails legitimately,
    # and new bytes reset the count. Parking is for content that stopped
    # changing and still cannot produce a completed run.
    attempts_file="${WATCH_STATE}/${sample}.attempts"
    prev_count=0
    prev_mtimes=""
    if [[ -f "${attempts_file}" ]]; then
      IFS=' ' read -r prev_count prev_mtimes < "${attempts_file}" || true
    fi
    if [[ "${prev_mtimes}" == "${mtimes}" ]]; then
      attempts=$(( prev_count + 1 ))
    else
      attempts=1
    fi
    echo "${attempts} ${mtimes}" > "${attempts_file}"
    if (( attempts > WATCH_DISPATCH_ATTEMPTS )); then
      echo "WARNING: parking ${sample} - dispatched ${WATCH_DISPATCH_ATTEMPTS} times without a"
      echo "         completed mosdepth run (its MD5 keeps failing or its jobs keep dying);"
      echo "         the final 01_download_and_mosdepth.sh sweep will handle it."
      touch "${WATCH_STATE}/${sample}.parked"
      PARKED=$(( PARKED + 1 ))
      continue
    fi
    if submit_mosdepth "${sample}" "${LINE_NUMS[${idx}]}" "${MD5S[${idx}]}"; then
      touch "${WATCH_STATE}/${sample}.dispatched"
      IN_FLIGHT=$(( IN_FLIGHT + 1 ))
      LAST_ACTIVITY=$(date +%s)
    else
      echo "WARNING: sbatch failed for ${sample}; will retry next round."
      WAITING=$(( WAITING + 1 ))
    fi
  done

  if (( WAITING == 0 && IN_FLIGHT == 0 )); then
    break
  fi
  if (( WATCH_IDLE_EXIT > 0 && WAITING > 0 && IN_FLIGHT == 0 )) \
     && (( $(date +%s) - LAST_ACTIVITY > WATCH_IDLE_EXIT )); then
    echo "Nothing in ${CRAM_DIR} has changed for ${WATCH_IDLE_EXIT} s with ${WAITING} sample(s)"
    echo "still missing - the transfer is finished or stalled. Run bash 01_download_and_mosdepth.sh to sweep."
    break
  fi
  sleep "${WATCH_INTERVAL}"
done
rm -f "${IN_FLIGHT_FILE}"

echo ""
echo "============================================================"
echo " Stage watcher finished: ${DONE} done, ${IN_FLIGHT} still in flight, ${WAITING} never arrived, ${PARKED} parked."
if (( PARKED > 0 )); then
  echo " Parked samples:"
  for parked_file in "${WATCH_STATE}"/*.parked; do
    [[ -f "${parked_file}" ]] && echo "   $(basename "${parked_file}" .parked)"
  done
fi
echo " Finish with: bash 01_download_and_mosdepth.sh   # verifies and sweeps the remainder"
echo "============================================================"
(( PARKED == 0 )) || exit 1
