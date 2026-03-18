#!/usr/bin/env bash
# =============================================================================
# 03_collect_qc.sh — Aggregate per-sample QC into a single overlay-ready table
# =============================================================================
#
# Collects three sources of publicly available or pipeline-generated QC:
#
#   1. mosdepth summary (.mosdepth.summary.txt) — generated during step 01.
#      Provides per-chromosome mean coverage.  From this we derive:
#        - Mean autosomal coverage (chr1–22)
#        - X and Y coverage ratios → inferred genetic sex
#
#   2. BAS files (.bam.bas) — per-sample Picard-equivalent stats downloaded
#      from EBI during step 01 alongside the CRAM.  Provides:
#        - % mapped reads
#        - Duplication rate
#        - Total sequenced bases
#
#   3. IGSR sample panel (integrated_call_samples_v3.*.ALL.ped) — downloaded
#      during setup.  Provides:
#        - Population (e.g. GBR, YRI, CHB)
#        - Superpopulation (AFR, AMR, EAS, EUR, SAS)
#        - Reported sex
#
# Output:
#   $QC_OUTPUT/sample_qc.tsv — one row per sample, ready to join with
#   svd.pcs.txt on SAMPLE_ID for PCA overlays.
#
# Usage:
#   bash 03_collect_qc.sh
#
#   # Or submit as a short SLURM job:
#   sbatch --mem=4G --time=01:00:00 --wrap="bash 03_collect_qc.sh"
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${QC_OUTPUT}"

OUT_TSV="${QC_OUTPUT}/sample_qc.tsv"
MOSDEPTH_SUMMARY_SUFFIX=".mosdepth.summary.txt"

echo "============================================================"
echo " Collecting per-sample QC"
echo "============================================================"
echo ""
echo " mosdepth summaries: ${MOSDEPTH_DIR}/*${MOSDEPTH_SUMMARY_SUFFIX}"
echo " BAS files:          ${BAS_DIR}/*.bam.bas"
echo " Sample panel:       ${PANEL_FILE}"
echo " Output:             ${OUT_TSV}"
echo ""

# ── Helper: parse mosdepth summary ──────────────────────────────────────────
# Summary file columns (tab-separated): chrom length bases mean min max
# We compute mean autosomal coverage from chr1–chr22, and X/Y coverage ratios.
parse_mosdepth_summary() {
  local summary_file="$1"
  awk 'BEGIN { FS="\t" }
    NR == 1 { next }            # skip header
    $1 ~ /^chr[0-9]+$/ {       # autosomal chromosomes (chr1–chr22)
      sum += $4 * $2; len += $2
    }
    $1 == "chrX" { x_mean = $4 }
    $1 == "chrY" { y_mean = $4 }
    END {
      auto_cov = (len > 0) ? sum / len : "NA"
      x_ratio  = (auto_cov > 0 && x_mean != "") ? x_mean / auto_cov : "NA"
      y_ratio  = (auto_cov > 0 && y_mean != "") ? y_mean / auto_cov : "NA"
      # Infer sex: males have X~0.5× autosome and detectable Y; females X~1.0×
      if (x_ratio != "NA" && y_ratio != "NA") {
        inferred_sex = (y_ratio > 0.1 && x_ratio < 0.75) ? "M" : "F"
      } else {
        inferred_sex = "NA"
      }
      printf "%.4f\t%.4f\t%.4f\t%s\n", auto_cov, x_ratio, y_ratio, inferred_sex
    }
  ' "${summary_file}"
}

# ── Helper: parse BAS file ───────────────────────────────────────────────────
# BAS files use a tab-separated format from the 1000G alignment pipeline.
# Column names vary slightly; we use a flexible key=value or positional parse.
# We extract: percent mapped, duplication rate, total bases.
parse_bas_file() {
  local bas_file="$1"
  awk 'BEGIN { FS="\t" }
    NR == 1 {
      # Find column indices by header name (case-insensitive)
      for (i = 1; i <= NF; i++) {
        h = tolower($i)
        if (h == "mapped_bases" || h == "mapped_bases_cigar") mapped_col = i
        if (h == "total_bases")                                total_col  = i
        if (h == "duplicate_mapped_bases" || h == "dup_mapped_bases") dup_col = i
      }
      next
    }
    NR == 2 {
      total  = (total_col  > 0) ? $total_col  : "NA"
      mapped = (mapped_col > 0) ? $mapped_col : "NA"
      dup    = (dup_col    > 0) ? $dup_col    : "NA"
      pct_mapped = (total != "NA" && total > 0 && mapped != "NA") \
                   ? mapped / total * 100 : "NA"
      pct_dup    = (mapped != "NA" && mapped > 0 && dup != "NA") \
                   ? dup / mapped * 100 : "NA"
      printf "%.4f\t%.4f\t%s\n",
        (pct_mapped != "NA") ? pct_mapped : 0,
        (pct_dup    != "NA") ? pct_dup    : 0,
        total
    }
  ' "${bas_file}"
}

# ── Build QC table ───────────────────────────────────────────────────────────
echo "Building QC table..."

HEADER="SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX\tPCT_MAPPED\tPCT_DUPLICATE\tTOTAL_BASES"
echo -e "${HEADER}" > "${OUT_TSV}"

# Iterate over mosdepth summary files (one per sample)
FOUND=0
MISSING_BAS=0
for summary in "${MOSDEPTH_DIR}"/*"${MOSDEPTH_SUMMARY_SUFFIX}"; do
  [[ -f "${summary}" ]] || continue

  # Derive sample ID from filename (strip suffix: .mosdepth.summary.txt)
  fname="$(basename "${summary}")"
  sample_id="${fname%%.*}"

  # Parse mosdepth summary
  mos_fields=$(parse_mosdepth_summary "${summary}")

  # Parse BAS file if available
  bas_file="${BAS_DIR}/${sample_id}.bam.bas"
  if [[ -f "${bas_file}" ]]; then
    bas_fields=$(parse_bas_file "${bas_file}")
  else
    bas_fields="NA\tNA\tNA"
    MISSING_BAS=$(( MISSING_BAS + 1 ))
  fi

  echo -e "${sample_id}\t${mos_fields}\t${bas_fields}" >> "${OUT_TSV}"
  FOUND=$(( FOUND + 1 ))
done

echo "  Processed ${FOUND} samples (${MISSING_BAS} missing BAS files)."

# ── Join with IGSR sample panel ──────────────────────────────────────────────
# The PED file columns: family_id sample_id paternal maternal sex phenotype pop superpop
# We extract: sample_id, pop, superpop, reported_sex (1=male,2=female)
echo "Joining with sample panel metadata..."

if [[ ! -f "${PANEL_FILE}" ]]; then
  echo "  WARNING: Sample panel not found at ${PANEL_FILE}."
  echo "  Run 00_setup.sh to download it, or set PANEL_FILE."
  echo "  Skipping population/sex annotation."
else
  PANEL_TSV="${QC_OUTPUT}/.panel_extracted.tmp"
  # Extract relevant columns from PED file (tab-separated, skip header if present)
  awk 'NR > 1 {
    sex_label = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    print $2 "\t" $7 "\t" $8 "\t" sex_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}" 2>/dev/null || \
  awk '{
    sex_label = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    print $2 "\t" $7 "\t" $8 "\t" sex_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}"

  # Join QC table with panel on SAMPLE_ID (column 1)
  MERGED_TSV="${QC_OUTPUT}/.merged.tmp"
  echo -e "SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX\tPCT_MAPPED\tPCT_DUPLICATE\tTOTAL_BASES\tPOPULATION\tSUPERPOPULATION\tREPORTED_SEX" \
    > "${MERGED_TSV}"

  # Use awk join: load panel into associative array, then annotate QC rows
  awk 'BEGIN { OFS="\t" }
    NR == FNR {
      pop[$1] = $2; super[$1] = $3; sex[$1] = $4
      next
    }
    FNR == 1 { next }   # skip QC header
    {
      sid = $1
      p  = (sid in pop)   ? pop[sid]   : "NA"
      s  = (sid in super) ? super[sid] : "NA"
      sx = (sid in sex)   ? sex[sid]   : "NA"
      print $0, p, s, sx
    }
  ' "${PANEL_TSV}" "${OUT_TSV}" >> "${MERGED_TSV}"

  mv "${MERGED_TSV}" "${OUT_TSV}"
  rm -f "${PANEL_TSV}"
  echo "  Population metadata joined."
fi

# ── Summary ──────────────────────────────────────────────────────────────────
NROWS=$(tail -n +2 "${OUT_TSV}" | wc -l)
echo ""
echo "============================================================"
echo " QC table complete: ${OUT_TSV}"
echo " Rows: ${NROWS} samples"
echo ""
echo " Columns:"
head -1 "${OUT_TSV}" | tr '\t' '\n' | nl
echo ""
echo " Preview (first 3 data rows):"
head -4 "${OUT_TSV}" | column -t -s$'\t'
echo ""
echo " To join with NGS-PCA results:"
printf "   join -t $'\\t' -1 1 -2 1 \\\n"
echo "     <(sort ${OUT_TSV}) \\"
echo "     <(sort ${NGSPCA_OUTPUT}/svd.pcs.txt) \\"
echo "     > ${QC_OUTPUT}/pcs_with_qc.tsv"
echo "============================================================"
