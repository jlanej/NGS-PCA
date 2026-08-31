#!/usr/bin/env bash
# =============================================================================
# 04b_fast_mode_report.sh — the fast-mode evaluation as one shareable HTML
# =============================================================================
#
# Re-runs 04_fast_mode_eval.py inside a small python+matplotlib container -
# so the figures exist even where the host python has no matplotlib - then
# assembles tables, figures, methods and generated interpretation into a
# self-contained fast_mode_report.html (print to PDF from any browser).
#
# When ${NGSPCA_SEED_CONTROL_OUTPUT}/svd.pcs.txt exists (the normal tree
# recomputed under a different -randomSeed), the subspace analysis also runs
# on that pair and the report gains a calibration section: fast-vs-normal
# differences judged against the estimator's own seed-to-seed noise.
#
# The container is provisioned like aria2's: pull the CI-published :eval tag
# (unprivileged, works on any node), fall back to building eval_report.def
# (needs fakeroot), and failing both, run with the host python3 - the report
# still builds, carrying tables without figures.
#
# Usage:
#   bash 04b_fast_mode_report.sh            # after both PCA runs and the QC
#   bash 04b_fast_mode_report.sh --pcs 50   # extra args reach the eval
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

# Self-defaulted so a config.sh copy predating the seed control still works
SEED_CONTROL_SEED="${SEED_CONTROL_SEED:-43}"
NGSPCA_SEED_CONTROL_OUTPUT="${NGSPCA_SEED_CONTROL_OUTPUT:-${WORK_DIR}/ngspca_output_seed${SEED_CONTROL_SEED}}"

EVAL_OUT="${QC_OUTPUT}/fast_mode_eval"
mkdir -p "${EVAL_OUT}"

provision_eval() {
  [[ -s "${EVAL_SIF}" ]] && return 0
  command -v apptainer &>/dev/null || return 1
  echo "Provisioning eval image (one-time)..."
  if apptainer pull "${EVAL_SIF}" "${EVAL_IMAGE_URI}" &>/dev/null && [[ -s "${EVAL_SIF}" ]]; then
    echo "  Pulled ${EVAL_IMAGE_URI}."
    return 0
  fi
  rm -f "${EVAL_SIF}"
  if [[ -f "${CONFIG_DIR}/eval_report.def" ]] \
     && apptainer build --fakeroot "${EVAL_SIF}" "${CONFIG_DIR}/eval_report.def" \
     && [[ -s "${EVAL_SIF}" ]]; then
    return 0
  fi
  rm -f "${EVAL_SIF}"
  return 1
}

# One python command array, containered when possible, host otherwise
PY=(python3)
if provision_eval; then
  BINDS=()
  for dir in "${NGSPCA_OUTPUT}" "${NGSPCA_FAST_OUTPUT}" "${NGSPCA_SEED_CONTROL_OUTPUT}" \
             "${QC_OUTPUT}" "${QC_FAST_OUTPUT}" \
             "${MOSDEPTH_TIMING_DIR}" "${CONFIG_DIR}"; do
    [[ -d "${dir}" ]] && BINDS+=(--bind "${dir}")
  done
  PY=(apptainer exec ${BINDS[@]+"${BINDS[@]}"} "${EVAL_SIF}" python3)
  echo "Running the evaluation in ${EVAL_SIF} (figures included)."
else
  echo "No eval container available; using the host python3 - figures appear only if it has matplotlib."
fi

"${PY[@]}" "${CONFIG_DIR}/04_fast_mode_eval.py" \
  --normal "${NGSPCA_OUTPUT}" \
  --fast "${NGSPCA_FAST_OUTPUT}" \
  --timing "${MOSDEPTH_TIMING_DIR}" \
  --qc-normal "${QC_OUTPUT}/sample_qc.tsv" \
  --qc-fast "${QC_FAST_OUTPUT}/sample_qc.tsv" \
  --out "${EVAL_OUT}" \
  "$@"

# Set-level concordance - principal angles, containment, cross-correlation -
# needs numpy, which rides along with matplotlib in the eval image. Skipped
# with a note otherwise; the report renders its absence honestly.
if "${PY[@]}" -c "import numpy" 2>/dev/null; then
  "${PY[@]}" "${CONFIG_DIR}/04c_pc_set_concordance.py" \
    --normal "${NGSPCA_OUTPUT}/svd.pcs.txt" \
    --fast "${NGSPCA_FAST_OUTPUT}/svd.pcs.txt" \
    --out-dir "${EVAL_OUT}"

  # Seed-control calibration: the SAME normal tree under a different
  # -randomSeed, run through the identical analysis. Differences there are
  # the estimator's own truncation noise - the yardstick the report holds the
  # fast-mode differences against.
  SEED_PCS="${NGSPCA_SEED_CONTROL_OUTPUT}/svd.pcs.txt"
  if [[ -s "${SEED_PCS}" ]]; then
    mkdir -p "${EVAL_OUT}/seed_control"
    "${PY[@]}" "${CONFIG_DIR}/04c_pc_set_concordance.py" \
      --normal "${NGSPCA_OUTPUT}/svd.pcs.txt" \
      --fast "${SEED_PCS}" \
      --label-b seed-control \
      --out-dir "${EVAL_OUT}/seed_control"
  else
    echo "Seed-control PCs not found (${SEED_PCS}) - the report will carry the"
    echo "fast-mode comparison without its calibration yardstick. Produce them with:"
    echo "  NGSPCA_OUTPUT=\"\${WORK_DIR}/ngspca_output_seed${SEED_CONTROL_SEED}\" RANDOM_SEED=${SEED_CONTROL_SEED} sbatch 02_run_ngspca.sh"
  fi
else
  echo "numpy not available - skipping the PC set-concordance section."
fi

"${PY[@]}" "${CONFIG_DIR}/04b_fast_mode_report.py" --eval-dir "${EVAL_OUT}"

echo ""
echo "Report: ${EVAL_OUT}/fast_mode_report.html"
echo "  Self-contained - copy it anywhere, open in a browser, print to PDF for a fixed copy."
