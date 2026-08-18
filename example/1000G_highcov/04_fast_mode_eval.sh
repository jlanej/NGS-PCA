#!/usr/bin/env bash
# =============================================================================
# 04_fast_mode_eval.sh — evaluate the mosdepth fast-mode comparison
# =============================================================================
#
# Runs 04_fast_mode_eval.py against the paths in config.sh, once both NGS-PCA
# runs of the comparison have finished (see the README's "Fast-mode
# comparison" section). Runs on the host: it needs python3, and produces the
# summary figure only where matplotlib is importable - the tables and report
# are written either way.
#
# Usage:
#   bash 04_fast_mode_eval.sh            # defaults from config.sh
#   bash 04_fast_mode_eval.sh --pcs 50   # extra args reach the python script
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

exec python3 "${SCRIPT_DIR}/04_fast_mode_eval.py" \
  --normal "${NGSPCA_OUTPUT}" \
  --fast "${NGSPCA_FAST_OUTPUT}" \
  --timing "${MOSDEPTH_TIMING_DIR}" \
  --out "${QC_OUTPUT}/fast_mode_eval" \
  "$@"
