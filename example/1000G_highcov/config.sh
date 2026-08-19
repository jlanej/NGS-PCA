#!/usr/bin/env bash
# =============================================================================
# config.sh — Shared configuration for the 1000G high-coverage NGS-PCA pipeline
# =============================================================================
#
# Source this file from every pipeline script:
#   source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
#
# Override any variable by exporting it before sourcing, e.g.:
#   export WORK_DIR=/scratch/user/1000G_highcov
#   source config.sh
# =============================================================================

set -euo pipefail

# Directory containing this config file (used for repository-relative defaults)
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Working directories (adjust to your HPC environment) ─────────────────────
WORK_DIR="${WORK_DIR:-/scratch/${USER}/1000G_highcov}"
CRAM_DIR="${CRAM_DIR:-${WORK_DIR}/crams}"
MOSDEPTH_DIR="${MOSDEPTH_DIR:-${WORK_DIR}/mosdepth_output}"
QC_OUTPUT="${QC_OUTPUT:-${WORK_DIR}/qc_output}"
NGSPCA_OUTPUT="${NGSPCA_OUTPUT:-${WORK_DIR}/ngspca_output}"
LOG_DIR="${LOG_DIR:-${WORK_DIR}/logs}"

# ── Container image ──────────────────────────────────────────────────────────
SIF_IMAGE="${SIF_IMAGE:-${WORK_DIR}/ngs-pca.sif}"

# ── Reference genome (GRCh38, required for CRAM decoding by mosdepth) ───────
REF_DIR="${REF_DIR:-${WORK_DIR}/reference}"
REF_FASTA="${REF_FASTA:-${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa}"
REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# ── 1000G high-coverage data (3,202 NYGC 30x WGS samples) ──────────────────
# NYGC 30x sequence indexes — one for the 2,504 unrelated samples, one for the
# 698 related samples.  Together they cover all 3,202 high-coverage samples.
#
# IMPORTANT: The GitHub-hosted file at igsr/1000Genomes_data_indexes (master)
# named "1000genomes.high_coverage.GRCh38DH.alignment.index" contains only the
# original 2015 PCR-free pilot data (24 samples, one per population).
# The full NYGC 30x indexes live on the EBI FTP under 1000G_2504_high_coverage.
#
# The sequence.index files list CRAMs on the ENA FTP (ftp.sra.ebi.ac.uk).
# Columns: ENA_FILE_PATH  MD5SUM  RUN_ID  ...  SAMPLE_NAME  POPULATION  ...
#
# See: https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
# Ref: Byrska-Bishop et al. (2022) Cell 185(18):3426-3440
INDEX_URL_2504="${INDEX_URL_2504:-ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index}"
INDEX_URL_698="${INDEX_URL_698:-ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index}"
INDEX_FILE_2504="${INDEX_FILE_2504:-${WORK_DIR}/1000G_2504_high_coverage.sequence.index}"
INDEX_FILE_698="${INDEX_FILE_698:-${WORK_DIR}/1000G_698_related_high_coverage.sequence.index}"

MANIFEST="${MANIFEST:-${WORK_DIR}/manifest.tsv}"
EXPECTED_MANIFEST_SAMPLES="${EXPECTED_MANIFEST_SAMPLES:-3202}"
MIN_MANIFEST_SAMPLES="${MIN_MANIFEST_SAMPLES:-3000}"
# CRAM URLs in the NYGC indexes point to the ENA FTP (ftp.sra.ebi.ac.uk)
EXPECTED_FTP_PREFIX="${EXPECTED_FTP_PREFIX:-ftp://ftp.sra.ebi.ac.uk/}"

# IGSR sample panel — population, reported sex, family relationships (publicly available)
# NOTE: The PED file does NOT contain a superpopulation column; superpopulation
#       is derived from the population code in 03_collect_qc.sh.
PANEL_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped"
PANEL_FILE="${PANEL_FILE:-${WORK_DIR}/igsr_sample_panel.ped}"

# EBI/ENA Aspera settings (high-speed FASP transfers — requires system ascp)
#
# Two things changed at EMBL-EBI and BOTH must be handled, so a site
# 'module load aspera' will not work on its own:
#   1. fasp.sra.ebi.ac.uk now runs OpenSSH 8.7 with all SHA-1 key exchange
#      removed and encrypt-then-MAC MACs only. ascp 3.9.x's libssh2 cannot
#      negotiate that, and reports the failure misleadingly as
#      "failed to authenticate" — so changing the -i key does not help.
#   2. The anonymous key asperaweb_id_dsa.openssh is no longer accepted; EBI
#      replaced it with an RSA key for its public accounts.
# 01_download_and_mosdepth.sh provisions both automatically on the submit host
# via install_aspera.sh; see aspera.def for the full detail.
# See: https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data
ASPERA_HOME="${ASPERA_HOME:-${WORK_DIR}/aspera}"
ASPERA_SIF="${ASPERA_SIF:-${ASPERA_HOME}/aspera.sif}"
# Falls back to any ascp on PATH, so a working site install is still used.
if [[ -z "${ASPERA_BIN:-}" ]]; then
  if [[ -x "${ASPERA_HOME}/bin/ascp" ]]; then
    ASPERA_BIN="${ASPERA_HOME}/bin/ascp"
  else
    ASPERA_BIN="ascp"
  fi
fi
ASPERA_SSH_KEY="${ASPERA_SSH_KEY:-${ASPERA_HOME}/etc/ebi_public.key}"
# Per-transfer FASP target rate. The aggregate asked of EBI is this times
# DOWNLOAD_SLOTS (6 x 100m ≈ 600 Mbit/s), and it has to fit both the site
# link and what EBI will serve one site: a retired design asked for 60 Gbit/s
# across 200 tasks and collapsed into packet loss and session timeouts, with
# sessions crawling at ~13 Mbit/s before dying. Raise one factor only after a
# sweep finishes without timeouts.
ASPERA_BANDWIDTH="${ASPERA_BANDWIDTH:-100m}"
# ascp attempts per file before falling back to HTTPS. Between attempts the
# partial stays on disk and -k 2 resumes it, so a session that dies 7 GB into
# a CRAM does not start over.
ASPERA_RETRIES="${ASPERA_RETRIES:-3}"
ASPERA_PORT=33001
# Set USE_ASPERA=0 to skip Aspera entirely: no image is built, and downloads
# go straight to the parallel-HTTPS paths below. Also the right setting when
# MD5 mismatches cluster on Aspera transfers (see the README's triage): one
# cluster measured 506 of 506 Aspera payloads corrupt while HTTPS ran clean,
# with ascp exiting success each time. 01 stops trusting Aspera within a task
# after its first mismatch, but a wave known to corrupt should not pay one
# wasted transfer per sample to rediscover it.
USE_ASPERA="${USE_ASPERA:-1}"

# ── Non-Aspera download tuning ──────────────────────────────────────────────
# ENA serves the same paths over HTTPS with Accept-Ranges: bytes, so a single
# CRAM can be pulled as N concurrent range requests. Measured off-site, 8
# streams ran ~3x a single stream. aria2c is preferred when present; otherwise
# the pipeline falls back to parallel curl ranges, then to single-stream wget.
# What ENA sees from the site at once is DOWNLOAD_SLOTS x this - 6 x 16 = 96
# connections, polite - and per-file speed is what sets the cohort's pace, so
# this stays high now that few files transfer at a time. (The retired 200-task
# array at 16 connections each presented 3,200 and got throttled.)
DOWNLOAD_CONNECTIONS="${DOWNLOAD_CONNECTIONS:-16}"
# Base used to rewrite ftp:// manifest URLs for the HTTP download paths.
ENA_HTTPS_BASE="${ENA_HTTPS_BASE:-https://ftp.sra.ebi.ac.uk}"
# NYGC CRAMs are on the ENA FTP; reference genome is on the 1000G FTP.
ENA_ASPERA_USER="era-fasp@fasp.sra.ebi.ac.uk"
EBI_ASPERA_USER="fasp-g1k@fasp.1000genomes.ebi.ac.uk"
ENA_FTP_BASE="ftp://ftp.sra.ebi.ac.uk"
EBI_FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk"

# ── Download manager controls ────────────────────────────────────────────────
# 01_download_and_mosdepth.sh is one long-running job holding this many
# concurrent transfers; each verified sample immediately becomes its own
# mosdepth job (01b), whose concurrency the scheduler decides. The aggregate
# asked of EBI is DOWNLOAD_SLOTS x per-file rate - few fast transfers beat
# many starved ones.
DOWNLOAD_SLOTS="${DOWNLOAD_SLOTS:-6}"
# Pause downloading when this many CRAMs sit on disk awaiting mosdepth, so a
# stalled queue cannot fill scratch (~16 GB each).
MAX_LOCAL_CRAMS="${MAX_LOCAL_CRAMS:-60}"
# Stop after this many dispatches; 0 means the whole manifest. For smoke tests.
DOWNLOAD_LIMIT="${DOWNLOAD_LIMIT:-0}"
# Manager state: per-sweep results, and the persistent Aspera distrust flag a
# corrupt payload sets (delete the flag file to try Aspera again).
DOWNLOAD_STATE_DIR="${DOWNLOAD_STATE_DIR:-${WORK_DIR}/download_state}"
# Resources for each spawned mosdepth job.
MOSDEPTH_MEM="${MOSDEPTH_MEM:-8G}"
MOSDEPTH_TIME="${MOSDEPTH_TIME:-06:00:00}"

# ── Stage watcher (01c_dispatch_staged.sh) ──────────────────────────────────
# Poll interval while watching a bulk transfer land in CRAM_DIR. A pair is
# dispatched once both files' mtimes are unchanged across one interval.
WATCH_INTERVAL="${WATCH_INTERVAL:-60}"
# Exit after this many seconds with samples still missing and nothing in
# flight - the transfer is finished or stalled either way (0 = wait forever).
WATCH_IDLE_EXIT="${WATCH_IDLE_EXIT:-7200}"
# Dispatches per sample before it is parked for the final sweep.
WATCH_DISPATCH_ATTEMPTS="${WATCH_DISPATCH_ATTEMPTS:-3}"

# ── mosdepth parameters ─────────────────────────────────────────────────────
MOSDEPTH_BIN_SIZE="${MOSDEPTH_BIN_SIZE:-1000}"
MOSDEPTH_THREADS="${MOSDEPTH_THREADS:-2}"

# ── mosdepth fast-mode comparison (optional) ────────────────────────────────
# COMPARE_FAST_MODE=1 makes 01_download_and_mosdepth.sh run mosdepth twice per
# sample - once as configured, once with --fast-mode - into two output trees,
# recording the wall time of each run. Download time is in neither measurement:
# both runs share the one CRAM the task just downloaded, and their order
# alternates by manifest line so page-cache warming cannot favour one mode in
# aggregate. A sample is processed when either tree lacks its output, so every
# timed pair comes from one node and one download. See the README's
# "Fast-mode comparison" section; 04_fast_mode_eval.sh evaluates the results.
COMPARE_FAST_MODE="${COMPARE_FAST_MODE:-0}"
MOSDEPTH_FAST_DIR="${MOSDEPTH_FAST_DIR:-${WORK_DIR}/mosdepth_output_fast}"
NGSPCA_FAST_OUTPUT="${NGSPCA_FAST_OUTPUT:-${WORK_DIR}/ngspca_output_fast}"
QC_FAST_OUTPUT="${QC_FAST_OUTPUT:-${WORK_DIR}/qc_output_fast}"
MOSDEPTH_TIMING_DIR="${MOSDEPTH_TIMING_DIR:-${QC_OUTPUT}/mosdepth_timing}"

# ── NGS-PCA parameters ──────────────────────────────────────────────────────
NUM_PC="${NUM_PC:-200}"
ITERS="${ITERS:-10}"
OVERSAMPLE="${OVERSAMPLE:-200}"
RANDOM_SEED="${RANDOM_SEED:-42}"
NGSPCA_THREADS="${NGSPCA_THREADS:-32}"
SAMPLE_EVERY="${SAMPLE_EVERY:-0}"
BED_EXCLUDE_REPO_DEFAULT="${CONFIG_DIR}/../../resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz"
BED_EXCLUDE_CONTAINER_DEFAULT="/app/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz"
if [[ -f "${BED_EXCLUDE_REPO_DEFAULT}" ]]; then
  BED_EXCLUDE_DEFAULT="${BED_EXCLUDE_REPO_DEFAULT}"
else
  BED_EXCLUDE_DEFAULT="${BED_EXCLUDE_CONTAINER_DEFAULT}"
fi
BED_EXCLUDE="${BED_EXCLUDE:-${BED_EXCLUDE_DEFAULT}}"
# Java heap size — set slightly below SBATCH --mem to leave OS headroom.
# For 3,202 samples × 200 PCs with SBATCH --mem=256G, 240g is appropriate.
XMX="${XMX:-240g}"
