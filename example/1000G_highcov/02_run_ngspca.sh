#!/usr/bin/env bash
# =============================================================================
# 02_run_ngspca.sh — Run NGS-PCA on all 3,202 1000G high-coverage samples
# =============================================================================
#
# Single SLURM job: computes ~200 PCs from the mosdepth output produced by
# 01_download_and_mosdepth.sh.
#
# Resource requirements depend on sample count:
#   - 3,202 samples × ~200 PCs: ~256 GB RAM, ~32 threads, ~6–12 h walltime
#   - Larger cohorts (e.g. 142k): ~1800 GB RAM, 120 threads, ~60 h
# Adjust #SBATCH directives below to match your cluster.
#
# Usage:
#   # After all array tasks from step 01 have completed:
#   sbatch 02_run_ngspca.sh
# =============================================================================

#SBATCH --job-name=1kG_ngspca
#SBATCH --output=logs/ngspca_%j.out
#SBATCH --error=logs/ngspca_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
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

# Ensure log and output directories exist
mkdir -p "${LOG_DIR}" "${NGSPCA_OUTPUT}"

echo "============================================================"
echo " NGS-PCA — 1000G High-Coverage (${NUM_PC} PCs)"
echo "============================================================"
echo ""
echo " Input:       ${MOSDEPTH_DIR}/"
echo " Output:      ${NGSPCA_OUTPUT}/"
echo " Parameters:"
echo "   numPC        = ${NUM_PC}"
echo "   iters        = ${ITERS}"
echo "   oversample   = ${OVERSAMPLE}"
echo "   randomSeed   = ${RANDOM_SEED}"
echo "   threads      = ${NGSPCA_THREADS}"
echo "   sampleEvery  = ${SAMPLE_EVERY}"
echo "   bedExclude   = ${BED_EXCLUDE}"
echo ""
echo " Started: $(date)"
echo ""

# ── Pre-flight check: verify mosdepth output ────────────────────────────────
echo "Checking mosdepth outputs..."
EXPECTED=$(tail -n +2 "${MANIFEST}" | wc -l)
FOUND=$(find "${MOSDEPTH_DIR}" -name "*.by${MOSDEPTH_BIN_SIZE}.regions.bed.gz" | wc -l)
echo "  Expected: ${EXPECTED} samples"
echo "  Found:    ${FOUND} mosdepth result files"
if [[ "${FOUND}" -eq 0 ]]; then
  echo "ERROR: No mosdepth output found. Run 01_download_and_mosdepth.sh first."
  exit 1
fi
if [[ "${FOUND}" -lt "${EXPECTED}" ]]; then
  echo "WARNING: ${FOUND}/${EXPECTED} samples have mosdepth output."
  echo "  Proceeding with available samples. Re-run step 01 for missing samples."
fi
echo ""

# ── Run NGS-PCA ──────────────────────────────────────────────────────────────
echo "Running NGS-PCA..."
apptainer run \
  --bind "${MOSDEPTH_DIR}":/mosdepth \
  --bind "${NGSPCA_OUTPUT}":/output \
  "${SIF_IMAGE}" \
  -input /mosdepth/ \
  -outputDir /output/ \
  -numPC "${NUM_PC}" \
  -sampleEvery "${SAMPLE_EVERY}" \
  -threads "${NGSPCA_THREADS}" \
  -iters "${ITERS}" \
  -oversample "${OVERSAMPLE}" \
  -randomSeed "${RANDOM_SEED}" \
  -distribution UNIFORM \
  -bedExclude "${BED_EXCLUDE}"

echo ""
echo "============================================================"
echo " NGS-PCA complete."
echo ""
echo " Output files:"
for f in svd.pcs.txt svd.loadings.txt svd.singularvalues.txt svd.bins.txt svd.samples.txt; do
  if [[ -f "${NGSPCA_OUTPUT}/${f}" ]]; then
    SIZE=$(du -h "${NGSPCA_OUTPUT}/${f}" | cut -f1)
    ROWS=$(wc -l < "${NGSPCA_OUTPUT}/${f}")
    echo "   ${f}  (${SIZE}, ${ROWS} lines)"
  else
    echo "   ${f}  MISSING"
  fi
done
echo ""
echo " Finished: $(date)"
echo "============================================================"
