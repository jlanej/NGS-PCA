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

# ── Working directories (adjust to your HPC environment) ─────────────────────
WORK_DIR="${WORK_DIR:-/scratch/${USER}/1000G_highcov}"
CRAM_DIR="${CRAM_DIR:-${WORK_DIR}/crams}"
MOSDEPTH_DIR="${MOSDEPTH_DIR:-${WORK_DIR}/mosdepth_output}"
BAS_DIR="${BAS_DIR:-${WORK_DIR}/bas_files}"   # per-sample Picard-equivalent QC
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

# IGSR sample panel — population, superpopulation, reported sex (publicly available)
PANEL_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped"
PANEL_FILE="${PANEL_FILE:-${WORK_DIR}/igsr_sample_panel.ped}"

# EBI/ENA Aspera settings (high-speed FASP transfers — requires system ascp)
# Install Aspera Connect from https://www.ibm.com/products/aspera/downloads
# or load it as an HPC module: module load aspera-connect
# See: https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data
ASPERA_SSH_KEY="${ASPERA_SSH_KEY:-${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh}"
ASPERA_BANDWIDTH="${ASPERA_BANDWIDTH:-300m}"
ASPERA_PORT=33001
# NYGC CRAMs are on the ENA FTP; reference genome is on the 1000G FTP.
ENA_ASPERA_USER="era-fasp@fasp.sra.ebi.ac.uk"
EBI_ASPERA_USER="fasp-g1k@fasp.1000genomes.ebi.ac.uk"
ENA_FTP_BASE="ftp://ftp.sra.ebi.ac.uk"
EBI_FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk"

# ── Download array batching/concurrency controls ────────────────────────────
# Number of manifest entries processed sequentially by each SLURM array task.
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-1}"
# Recommended max number of concurrently running array tasks when submitting:
#   sbatch --array=1-N%${MAX_CONCURRENT_TASKS} 01_download_and_mosdepth.sh
MAX_CONCURRENT_TASKS="${MAX_CONCURRENT_TASKS:-200}"

# ── mosdepth parameters ─────────────────────────────────────────────────────
MOSDEPTH_BIN_SIZE="${MOSDEPTH_BIN_SIZE:-1000}"
MOSDEPTH_THREADS="${MOSDEPTH_THREADS:-2}"

# ── NGS-PCA parameters ──────────────────────────────────────────────────────
NUM_PC="${NUM_PC:-200}"
ITERS="${ITERS:-10}"
OVERSAMPLE="${OVERSAMPLE:-200}"
RANDOM_SEED="${RANDOM_SEED:-42}"
NGSPCA_THREADS="${NGSPCA_THREADS:-32}"
SAMPLE_EVERY="${SAMPLE_EVERY:-0}"
BED_EXCLUDE="${BED_EXCLUDE:-/app/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz}"
# Java heap size — set slightly below SBATCH --mem to leave OS headroom.
# For 3,202 samples × 200 PCs with SBATCH --mem=256G, 240g is appropriate.
XMX="${XMX:-240g}"
