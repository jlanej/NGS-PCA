#!/usr/bin/env bash
# =============================================================================
# 01_download_and_mosdepth.sh — Download manager: fetch CRAMs, spawn mosdepth
# =============================================================================
#
# One long-running job downloads the cohort a few samples at a time -
# DOWNLOAD_SLOTS concurrent transfers, each with the full Aspera->HTTPS
# fallback chain and MD5 verification - and the moment a sample verifies, its
# mosdepth run is submitted as an independent SLURM job (01b_mosdepth_sample.sh).
# The WAN is the scarce, failure-prone resource, so it gets few gentle streams;
# mosdepth is abundant compute, so its concurrency is whatever the scheduler
# grants and never holds a download slot.
#
# The sweep is idempotent: samples with mosdepth output are skipped, samples
# whose mosdepth job is already queued or running are skipped, and everything
# else - including a CRAM already on disk from a dead run - is verified and
# dispatched. Re-running this script continues where the last sweep stopped.
#
# Usage:
#   bash 01_download_and_mosdepth.sh          # submits the manager as one job
#   DOWNLOADER_LOCAL=1 bash 01_download_and_mosdepth.sh   # run in place
#                                             # (login/DTN node, tmux advised)
#   DOWNLOAD_LIMIT=10 bash 01_download_and_mosdepth.sh    # smoke-test subset
# =============================================================================

#SBATCH --job-name=1kG_download
#SBATCH --output=logs/download_manager_%j.out
#SBATCH --error=logs/download_manager_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
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

# Attempt to load the Aspera module if the HPC module system is available.
# This is a no-op when 'module' is not present or 'aspera' is not a module.
module load aspera 2>/dev/null || true

mkdir -p "${LOG_DIR}" "${CRAM_DIR}" "${MOSDEPTH_DIR}" "${DOWNLOAD_STATE_DIR}"

# Inside a SLURM job the batch script runs from the spool directory, so 01b is
# addressed through the repository directory config.sh resolved itself in.
MOSDEPTH_JOB_SCRIPT="${CONFIG_DIR}/01b_mosdepth_sample.sh"

# Persistent distrust flag: a completed Aspera transfer that failed its MD5 is
# systematic corruption, not bad luck (measured at 506 of 506 on one cluster),
# so one bad payload switches the whole run - and later sweeps - to HTTPS.
# Delete the file to give Aspera another chance.
ASPERA_DISTRUST_FLAG="${DOWNLOAD_STATE_DIR}/aspera_disabled"

# Returns the expected mosdepth output path for a given sample ID.
mosdepth_output_path() {
  echo "${MOSDEPTH_DIR}/${1}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
}

# Fast-mode twin of the above. Same file name, so NGS-PCA derives identical
# sample identifiers from either tree and the two PC tables join by sample.
mosdepth_fast_output_path() {
  echo "${MOSDEPTH_FAST_DIR}/${1}.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz"
}

# Whether a sample still needs downloading or mosdepth. With the fast-mode
# comparison on, a sample missing either output is redone in full so that
# every timed pair comes from one node and one downloaded CRAM.
needs_processing() {
  local sample="$1"
  if [[ "${COMPARE_FAST_MODE}" == "1" ]]; then
    [[ ! -s "$(mosdepth_output_path "${sample}")" || ! -s "$(mosdepth_fast_output_path "${sample}")" ]]
  else
    [[ ! -s "$(mosdepth_output_path "${sample}")" ]]
  fi
}

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest missing or empty: ${MANIFEST}"
  echo "Run setup first: bash 00_setup.sh"
  exit 1
fi
TOTAL_SAMPLES=$(tail -n +2 "${MANIFEST}" | wc -l | tr -d "[:space:]")
if (( TOTAL_SAMPLES < MIN_MANIFEST_SAMPLES )); then
  echo "ERROR: Manifest has only ${TOTAL_SAMPLES} samples (minimum expected: ${MIN_MANIFEST_SAMPLES})."
  echo "This likely indicates an incomplete manifest/index download."
  echo "Rebuild it with: rm -f \"${MANIFEST}\" \"${INDEX_FILE_2504}\" \"${INDEX_FILE_698}\" && bash 00_setup.sh"
  exit 1
fi
if (( TOTAL_SAMPLES != EXPECTED_MANIFEST_SAMPLES )); then
  echo "WARNING: Manifest has ${TOTAL_SAMPLES} samples (expected ${EXPECTED_MANIFEST_SAMPLES})."
fi

if [[ -z "${SLURM_JOB_ID:-}" && "${DOWNLOADER_LOCAL:-0}" != "1" ]]; then
  # ── Submission mode ────────────────────────────────────────────────────────
  # Provision Aspera on the submit host - the build downloads ~68 MB and
  # writes one shared .sif, better done once here than inside the job.
  if [[ "${USE_ASPERA}" == "1" && ! -x "${ASPERA_HOME}/bin/ascp" ]]; then
    echo "Provisioning Aspera (one-time)..."
    if ! bash "${SCRIPT_DIR}/install_aspera.sh"; then
      echo "  Continuing without Aspera; downloads will use parallel HTTPS."
    fi
  fi
  # aria2c the same way when the host lacks it: a small bespoke image, built
  # once on the submit host. Non-fatal - parallel curl remains the fallback.
  if ! command -v aria2c &>/dev/null && command -v apptainer &>/dev/null \
     && [[ ! -s "${ARIA2_SIF}" ]]; then
    echo "Provisioning aria2 image (one-time)..."
    if ! apptainer build --fakeroot "${ARIA2_SIF}" "${SCRIPT_DIR}/aria2.def"; then
      rm -f "${ARIA2_SIF}"
      echo "  Continuing without aria2c; HTTPS downloads will use parallel curl."
    fi
  fi
  echo "Submitting the download manager (${DOWNLOAD_SLOTS} concurrent downloads)..."
  echo "  Manager log: ${LOG_DIR}/download_manager_<jobid>.out"
  echo "  Per-sample download logs: ${LOG_DIR}/download_<sample>.log"
  exec sbatch \
    --output="${LOG_DIR}/download_manager_%j.out" \
    --error="${LOG_DIR}/download_manager_%j.err" \
    "${BASH_SOURCE[0]}"
fi

# ── Stage 1 transports ───────────────────────────────────────────────────────

download_aspera() {
  # Download a file via Aspera. Supports both EBI 1000G FTP and ENA SRA FTP URLs.
  #   ftp://ftp.sra.ebi.ac.uk/vol1/run/...         → era-fasp@fasp.sra.ebi.ac.uk:vol1/run/...
  #   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/...  → fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/...
  local ftp_url="$1"
  local dest="$2"
  local aspera_path aspera_user

  [[ "${USE_ASPERA}" == "1" ]] || return 1
  # File-backed so every download slot, and every later sweep, sees it.
  [[ -f "${ASPERA_DISTRUST_FLAG}" ]] && return 1
  [[ "${ASPERA_DISABLED:-0}" == "1" ]] && return 1
  # ASPERA_BIN is resolved in config.sh: work-dir install if install_aspera.sh
  # has been run, else whatever ascp is on PATH.
  command -v "${ASPERA_BIN}" &>/dev/null || return 1
  [[ -s "${ASPERA_SSH_KEY}" ]] || return 1

  if [[ "${ftp_url}" == *"ftp.sra.ebi.ac.uk"* ]]; then
    aspera_path="${ftp_url//ftp:\/\/ftp.sra.ebi.ac.uk\//}"
    aspera_user="${ENA_ASPERA_USER}"
  elif [[ "${ftp_url}" == *"ftp.1000genomes.ebi.ac.uk"* ]]; then
    aspera_path="${ftp_url//ftp:\/\/ftp.1000genomes.ebi.ac.uk\//}"
    aspera_user="${EBI_ASPERA_USER}"
  else
    return 1
  fi

  local remote_basename attempt
  remote_basename="$(basename "${aspera_path}")"
  local downloaded="${CRAM_DIR}/${remote_basename}"
  # -k 2 resumes a partial from an earlier attempt after checksumming the
  # blocks already on disk. Only ascp may resume its own partials: FASP writes
  # blocks out of order, so a byte-range transport continuing the file would
  # seal holes into it, to be caught much later by the MD5 check. That is why
  # the partial is kept between attempts here and deleted before falling back.
  for (( attempt = 1; attempt <= ASPERA_RETRIES; attempt++ )); do
    if "${ASPERA_BIN}" -i "${ASPERA_SSH_KEY}" \
      -Tr -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" -L- -k 2 \
      "${aspera_user}:${aspera_path}" \
      "${CRAM_DIR}/"; then
      # Rename to the desired local filename if it differs from the remote name
      if [[ "${downloaded}" != "${dest}" && -f "${downloaded}" ]]; then
        mv "${downloaded}" "${dest}"
      fi
      return 0
    fi
    if (( attempt < ASPERA_RETRIES )); then
      echo "  Aspera attempt ${attempt}/${ASPERA_RETRIES} failed; retrying with resume..."
      sleep 15
    fi
  done
  # Clean up the partial before the byte-range fallbacks
  rm -f "${downloaded}" "${dest}"
  return 1
}

# Rewrite an ENA ftp:// manifest URL to its HTTPS equivalent. ENA serves the
# same tree over HTTPS with range support, which is what makes the parallel
# download paths below possible.
https_url_for() {
  local ftp_url="$1"
  echo "${ftp_url/ftp:\/\/ftp.sra.ebi.ac.uk/${ENA_HTTPS_BASE}}"
}

# The AWS Open Data mirror of the same files, or empty when the URL does not
# map or the mirror is disabled. Verified byte-identical to ENA's copies
# (content lengths match exactly, and every file passes the manifest MD5), it
# serves with S3 throughput and free egress - the preferred source, with ENA
# as fallback. The layout drops ENA's run-prefix grouping level and splits the
# release batches:
#   vol1/run/ERR324/ERR3239480/NA12718.final.cram (batch 2504)
#     -> 1000G_2504_high_coverage/data/ERR3239480/NA12718.final.cram
#   batch 698 -> 1000G_2504_high_coverage/additional_698_related/data/...
s3_url_for() {
  local ftp_url="$1" batch="$2"
  [[ -n "${S3_HTTPS_BASE}" ]] || { echo ""; return 0; }
  if [[ "${ftp_url}" != ftp://ftp.sra.ebi.ac.uk/vol1/run/*/*/* ]]; then
    echo ""
    return 0
  fi
  local rest="${ftp_url#ftp://ftp.sra.ebi.ac.uk/vol1/run/}"
  local run_and_file="${rest#*/}"
  case "${batch}" in
    2504) echo "${S3_HTTPS_BASE}/1000G_2504_high_coverage/data/${run_and_file}" ;;
    698) echo "${S3_HTTPS_BASE}/1000G_2504_high_coverage/additional_698_related/data/${run_and_file}" ;;
    *) echo "" ;;
  esac
}

# Content-Length of a remote file, or empty if it cannot be determined.
# Uses tolower() rather than gawk's IGNORECASE so this also works under the
# BSD awk and mawk found on some systems.
remote_size() {
  curl -sS -I -L --max-time 60 "$1" 2>/dev/null \
    | awk '{ if (tolower($0) ~ /^content-length:/) v = $2 }
           END { gsub(/\r/, "", v); print v }'
}

local_size() {
  stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null
}

# aria2c as a command array: the host's when present, else the bespoke
# aria2.sif built at submission (see aria2.def - kept out of the analysis
# image on purpose). The same-path bind makes --dir/--out resolve identically
# inside the container; apptainer shares the host network.
ARIA2C=()
ARIA2C_HOW=""
resolve_aria2c() {
  if command -v aria2c &>/dev/null; then
    ARIA2C=(aria2c)
    ARIA2C_HOW="host"
    return 0
  fi
  if command -v apptainer &>/dev/null && [[ -s "${ARIA2_SIF}" ]] \
     && apptainer exec "${ARIA2_SIF}" aria2c --version &>/dev/null; then
    ARIA2C=(apptainer exec --bind "${CRAM_DIR}" "${ARIA2_SIF}" aria2c)
    ARIA2C_HOW="via ${ARIA2_SIF} (apptainer)"
    return 0
  fi
  return 1
}

download_aria2() {
  local url="$1"
  local dest="$2"

  (( ${#ARIA2C[@]} )) || return 1

  if ! "${ARIA2C[@]}" \
    --quiet=true \
    --continue=true \
    --max-connection-per-server="${DOWNLOAD_CONNECTIONS}" \
    --split="${DOWNLOAD_CONNECTIONS}" \
    --min-split-size=1M \
    --max-tries=5 \
    --retry-wait=10 \
    --auto-file-renaming=false \
    --allow-overwrite=true \
    --dir "$(dirname "${dest}")" \
    --out "$(basename "${dest}")" \
    "${url}"; then
    rm -f "${dest}" "${dest}.aria2"
    return 1
  fi
  rm -f "${dest}.aria2"
}

# Parallel byte-range download using only curl, for hosts without aria2c.
download_curl_parallel() {
  local url="$1"
  local dest="$2"
  local total chunk parts_dir rc=0 i pid

  command -v curl &>/dev/null || return 1

  total="$(remote_size "${url}")"
  [[ "${total}" =~ ^[0-9]+$ ]] || return 1
  (( total > 0 )) || return 1

  local n="${DOWNLOAD_CONNECTIONS}"
  (( n < 1 )) && n=1
  chunk=$(( (total + n - 1) / n ))
  parts_dir="${dest}.parts.$$"
  rm -rf "${parts_dir}"
  mkdir -p "${parts_dir}" || return 1

  local pids=()
  for (( i = 0; i < n; i++ )); do
    local start=$(( i * chunk ))
    local end=$(( start + chunk - 1 ))
    (( start >= total )) && break
    (( end >= total )) && end=$(( total - 1 ))
    # Zero-padded part names so the assembly glob sorts numerically.
    curl -sS -L --fail --retry 5 --retry-delay 10 \
      -r "${start}-${end}" \
      -o "${parts_dir}/part.$(printf '%04d' "${i}")" \
      "${url}" &
    pids+=("$!")
  done

  for pid in "${pids[@]}"; do
    wait "${pid}" || rc=1
  done
  if (( rc != 0 )); then
    rm -rf "${parts_dir}"
    return 1
  fi

  if ! cat "${parts_dir}"/part.* > "${dest}"; then
    rm -rf "${parts_dir}" ; rm -f "${dest}"
    return 1
  fi
  rm -rf "${parts_dir}"

  # Guard against a short/torn assembly before mosdepth ever sees the file.
  local got
  got="$(local_size "${dest}")"
  if [[ "${got}" != "${total}" ]]; then
    echo "  ERROR: size mismatch after parallel download (expected ${total}, got ${got:-0})"
    rm -f "${dest}"
    return 1
  fi
}

download_wget() {
  local ftp_url="$1"
  local dest="$2"
  if ! wget -q -c -O "${dest}" "${ftp_url}"; then
    rm -f "${dest}"
    return 1
  fi
}

# Try each source (S3 mirror first, then ENA) with each transport in
# descending order of speed. Returns non-zero only when everything failed.
download_file() {
  local url="$1" dest="$2" label="$3" batch="${4:-}"

  if download_aspera "${url}" "${dest}"; then
    echo "  ${label}: Aspera download complete."
    return 0
  fi

  local sources=() source_names=() https
  https="$(s3_url_for "${url}" "${batch}")"
  if [[ -n "${https}" ]]; then
    sources+=("${https}")
    source_names+=("S3")
  fi
  sources+=("$(https_url_for "${url}")")
  source_names+=("ENA")

  local i
  for i in "${!sources[@]}"; do
    if download_aria2 "${sources[${i}]}" "${dest}"; then
      echo "  ${label}: aria2c download complete (${source_names[${i}]}, ${DOWNLOAD_CONNECTIONS} streams)."
      return 0
    fi
    if download_curl_parallel "${sources[${i}]}" "${dest}"; then
      echo "  ${label}: parallel curl download complete (${source_names[${i}]}, ${DOWNLOAD_CONNECTIONS} ranges)."
      return 0
    fi
  done
  if download_wget "${url}" "${dest}"; then
    echo "  ${label}: wget download complete (single stream)."
    return 0
  fi
  echo "  ERROR: all download methods failed for ${url}"
  return 1
}

# ── One download slot: fetch, verify, hand off to mosdepth ──────────────────

# download_sample <sample> <line_num> <cram_url> <crai_url> <md5> [batch]
download_sample() {
  local sample="$1" line_num="$2" cram_url="$3" crai_url="$4" cram_md5="$5" batch="${6:-}"
  local local_cram="${CRAM_DIR}/${sample}.cram"
  local local_crai="${CRAM_DIR}/${sample}.cram.crai"

  echo "============================================================"
  echo " Sample: ${sample} (manifest line ${line_num})"
  echo " CRAM: ${cram_url}"
  echo " Started: $(date)"
  echo "============================================================"

  # A file that fails MD5 - a transfer corrupted in flight, or a partial left
  # by a dead run and picked up as "already present" - is deleted and
  # downloaded once more before the sample is given up on.
  local download_attempt actual_md5
  for download_attempt in 1 2; do
    if [[ ! -f "${local_cram}" ]]; then
      echo "[1/2] Downloading CRAM..."
      if ! download_file "${cram_url}" "${local_cram}" "CRAM" "${batch}"; then
        return 1
      fi
    else
      echo "[1/2] CRAM already present: ${local_cram}"
    fi

    if [[ ! -f "${local_crai}" ]]; then
      echo "  Downloading CRAI..."
      if ! download_file "${crai_url}" "${local_crai}" "CRAI" "${batch}"; then
        rm -f "${local_cram}"
        return 1
      fi
    else
      echo "  CRAI already present: ${local_crai}"
    fi

    # Verify CRAM MD5 when the manifest carries one
    if [[ -n "${cram_md5}" ]]; then
      echo "  Verifying CRAM MD5..."
      actual_md5=$(md5sum "${local_cram}" | awk '{print $1}')
      if [[ "${actual_md5}" != "${cram_md5}" ]]; then
        echo "  WARNING: MD5 mismatch (expected: ${cram_md5}, got: ${actual_md5})"
        rm -f "${local_cram}" "${local_crai}"
        if (( download_attempt == 1 )); then
          touch "${ASPERA_DISTRUST_FLAG}"
          echo "  Removing corrupt file and re-downloading over HTTPS (Aspera distrusted from here on)."
          continue
        fi
        echo "  ERROR: MD5 mismatch again after a fresh download."
        return 1
      fi
      echo "  MD5 verified."
    fi
    break
  done

  echo "[2/2] Submitting mosdepth job..."
  local try jobid
  for try in 1 2 3; do
    if jobid=$(sbatch --parsable \
        --job-name="1kG_mosdepth_${sample}" \
        --output="${LOG_DIR}/mosdepth_${sample}_%j.out" \
        --error="${LOG_DIR}/mosdepth_${sample}_%j.err" \
        --cpus-per-task="${MOSDEPTH_THREADS}" \
        --mem="${MOSDEPTH_MEM}" \
        --time="${MOSDEPTH_TIME}" \
        "${MOSDEPTH_JOB_SCRIPT}" "${sample}" "${line_num}"); then
      echo "  mosdepth job ${jobid} submitted for ${sample}."
      return 0
    fi
    echo "  sbatch failed (attempt ${try}/3); retrying in 30 s..."
    sleep 30
  done
  echo "  ERROR: could not submit the mosdepth job for ${sample}; the verified CRAM is kept"
  echo "         and the next sweep will dispatch it without re-downloading."
  return 1
}

# ── Manager ──────────────────────────────────────────────────────────────────

RESULT_DIR="${DOWNLOAD_STATE_DIR}/results"
rm -rf "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"

# One scheduler snapshot: samples whose mosdepth job is already queued or
# running are not re-dispatched, so a manager restart cannot double-download.
# A file rather than an associative array keeps this bash-3 clean.
IN_FLIGHT_FILE="${DOWNLOAD_STATE_DIR}/in_flight.$$"
: > "${IN_FLIGHT_FILE}"
if command -v squeue &>/dev/null; then
  squeue -h -u "${USER:-$(id -un)}" -o "%j" 2>/dev/null \
    | sed -n 's/^1kG_mosdepth_//p' > "${IN_FLIGHT_FILE}" || true
fi

NEEDED_LINES=()
NEEDED_SAMPLES=()
SKIPPED_DONE=0
SKIPPED_IN_FLIGHT=0
LINE_NUM=1
while IFS=$'\t' read -r SAMPLE_ID CRAM_URL CRAI_URL CRAM_MD5 _; do
  LINE_NUM=$(( LINE_NUM + 1 ))   # manifest line, counting its header
  [[ -z "${SAMPLE_ID}" ]] && continue
  if ! needs_processing "${SAMPLE_ID}"; then
    SKIPPED_DONE=$(( SKIPPED_DONE + 1 ))
    continue
  fi
  if grep -qxF "${SAMPLE_ID}" "${IN_FLIGHT_FILE}"; then
    SKIPPED_IN_FLIGHT=$(( SKIPPED_IN_FLIGHT + 1 ))
    continue
  fi
  if [[ "${CRAM_URL}" != "${EXPECTED_FTP_PREFIX}"* || "${CRAI_URL}" != "${EXPECTED_FTP_PREFIX}"* ]]; then
    echo "ERROR: Manifest entry at line ${LINE_NUM} has unsupported CRAM/CRAI source: ${CRAM_URL}"
    echo "fail" > "${RESULT_DIR}/${SAMPLE_ID}"
    continue
  fi
  NEEDED_LINES+=("${LINE_NUM}")
  NEEDED_SAMPLES+=("${SAMPLE_ID}")
done < <(tail -n +2 "${MANIFEST}")
rm -f "${IN_FLIGHT_FILE}"

echo "Download manager: ${TOTAL_SAMPLES} samples in the manifest."
# Which mosdepth work each dispatched job will do is decided here by
# environment, and getting it wrong is only discovered after the CRAMs are
# gone - so it is announced, loudly, every run.
if [[ "${COMPARE_FAST_MODE}" == "1" ]]; then
  echo "  Fast-mode comparison: ON - each sample runs mosdepth twice (normal + fast, timed)."
else
  echo "  Fast-mode comparison: off - single mosdepth run per sample. If this cohort is for"
  echo "  the comparison, stop and set COMPARE_FAST_MODE=1 (export it, or pin it in config.sh)."
fi
resolve_aria2c || true
if (( ${#ARIA2C[@]} )); then
  echo "  aria2c: ${ARIA2C_HOW}"
else
  echo "  aria2c: not available - HTTPS transfers use parallel curl (${DOWNLOAD_CONNECTIONS} ranges)"
fi
echo "  Already done:        ${SKIPPED_DONE}"
echo "  mosdepth in flight:  ${SKIPPED_IN_FLIGHT}"
echo "  To download:         ${#NEEDED_LINES[@]} (${DOWNLOAD_SLOTS} at a time)"
if (( DOWNLOAD_LIMIT > 0 )); then
  echo "  DOWNLOAD_LIMIT=${DOWNLOAD_LIMIT}: stopping after that many dispatches."
fi

if (( ${#NEEDED_LINES[@]} == 0 )); then
  echo "Nothing to download."
  exit 0
fi

DISPATCHED=0
for idx in "${!NEEDED_LINES[@]}"; do
  if (( DOWNLOAD_LIMIT > 0 && DISPATCHED >= DOWNLOAD_LIMIT )); then
    break
  fi

  line_num="${NEEDED_LINES[${idx}]}"
  sample="${NEEDED_SAMPLES[${idx}]}"

  # Keep at most DOWNLOAD_SLOTS transfers running. Polling rather than
  # wait -n: transfers run for minutes, and bash older than 4.3 has no wait -n.
  while (( $(jobs -rp | wc -l) >= DOWNLOAD_SLOTS )); do
    sleep 2
  done

  # Disk backpressure: a stalled mosdepth queue must not fill scratch. Only
  # completed jobs delete their CRAMs, so the count of *.cram files bounds
  # local usage at roughly MAX_LOCAL_CRAMS x 16 GB. It gates only samples that
  # would DOWNLOAD: one whose CRAM is already on disk - staged by Globus, or
  # left by a dead run - consumes no new space, and dispatching it frees some,
  # so pausing those would deadlock a pre-staged cohort against itself.
  if [[ ! -f "${CRAM_DIR}/${sample}.cram" ]]; then
    while (( $(find "${CRAM_DIR}" -maxdepth 1 -name "*.cram" 2>/dev/null | wc -l) >= MAX_LOCAL_CRAMS )); do
      echo "Backpressure: ${MAX_LOCAL_CRAMS} CRAMs on disk await mosdepth; pausing downloads for 60 s..."
      sleep 60
    done
  fi
  line=$(sed -n "${line_num}p" "${MANIFEST}")
  IFS=$'\t' read -r _ CRAM_URL CRAI_URL CRAM_MD5 RELEASE_BATCH _ <<< "${line}"

  echo "[$(date '+%H:%M:%S')] slot -> ${sample} (log: ${LOG_DIR}/download_${sample}.log)"
  (
    if download_sample "${sample}" "${line_num}" "${CRAM_URL}" "${CRAI_URL}" "${CRAM_MD5}" "${RELEASE_BATCH}"; then
      echo "ok" > "${RESULT_DIR}/${sample}"
    else
      echo "fail" > "${RESULT_DIR}/${sample}"
    fi
  ) > "${LOG_DIR}/download_${sample}.log" 2>&1 &
  DISPATCHED=$(( DISPATCHED + 1 ))
done
wait || true

OK_COUNT=$(grep -lx "ok" "${RESULT_DIR}"/* 2>/dev/null | wc -l | tr -d "[:space:]") || true
FAIL_COUNT=$(grep -lx "fail" "${RESULT_DIR}"/* 2>/dev/null | wc -l | tr -d "[:space:]") || true
echo ""
echo "============================================================"
echo " Download sweep finished: ${OK_COUNT} dispatched to mosdepth, ${FAIL_COUNT} failed."
if (( FAIL_COUNT > 0 )); then
  echo " Failed samples (see ${LOG_DIR}/download_<sample>.log):"
  grep -lx "fail" "${RESULT_DIR}"/* | while IFS= read -r f; do echo "   $(basename "${f}")"; done
  echo " Re-run this script to retry them."
fi
if [[ -f "${ASPERA_DISTRUST_FLAG}" ]]; then
  echo " NOTE: Aspera is distrusted (a completed transfer failed MD5); downloads are HTTPS-only."
  echo "       Remove ${ASPERA_DISTRUST_FLAG} to try Aspera again."
fi
echo " mosdepth jobs run independently; watch them with: squeue -u \$USER -n 1kG_download; squeue -u \$USER | grep 1kG_mosdepth"
echo "============================================================"
(( FAIL_COUNT == 0 )) || exit 1
