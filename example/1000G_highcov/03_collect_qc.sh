#!/usr/bin/env bash
# =============================================================================
# 03_collect_qc.sh — Aggregate per-sample QC into a single overlay-ready table
# =============================================================================
#
# Collects publicly available or pipeline-generated QC from multiple sources:
#
#   1. mosdepth summary (.mosdepth.summary.txt) — generated during step 01.
#      Provides per-chromosome mean coverage.  From this we derive:
#        - Mean autosomal coverage (chr1–22)
#        - X and Y coverage ratios → inferred genetic sex
#        - Mitochondrial (chrM) coverage ratio — useful for sample QC
#
#   1b. mosdepth global distribution (.mosdepth.global.dist.txt) — generated
#       during step 01.  Cumulative coverage depth distribution.  From this:
#        - Median genome coverage
#        - % genome at ≥ 10× depth
#        - % genome at ≥ 20× depth
#
#   1c. mosdepth coverage summary (mosdepth_coverage_summary.tsv) — generated
#       by 03a_mosdepth_coverage_summary.sh.  Per-bin autosomal coverage
#       distribution statistics:
#        - Median bin coverage, SD, MAD, IQR
#        - HQ (high-quality) variants of each stat — bins not overlapping
#          the exclusion BED used in step 2
#        - mtDNA CN estimated from HQ autosomal median and chrM coverage
#
#   2. IGSR sample panel (integrated_call_samples_v3.*.ALL.ped) — downloaded
#      during setup.  PED file has 12 columns; we extract:
#        - Population ($7, e.g. GBR, YRI, CHB)
#        - Superpopulation — derived from population via lookup (AFR/AMR/EAS/EUR/SAS)
#        - Reported sex ($5, 1=M / 2=F)
#        - Family role ($8, PED Relationship: father, mother, child, unrel …)
#        - Relatedness — derived: both parental IDs "0" → unrelated, else related
#
#   3. Sample manifest (manifest.tsv) — built during setup from the two NYGC
#      sequence.index files.  Provides batch-level annotations:
#        - Release batch (2504 vs 698)
#        - Sequencing center
#        - Study ID
#        - Instrument model
#        - Library name
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

mkdir -p "${QC_OUTPUT}"

OUT_TSV="${QC_OUTPUT}/sample_qc.tsv"
MOSDEPTH_SUMMARY_SUFFIX=".mosdepth.summary.txt"
MOSDEPTH_GLOBAL_DIST_SUFFIX=".mosdepth.global.dist.txt"
COV_SUMMARY_TSV="${QC_OUTPUT}/mosdepth_coverage_summary.tsv"

echo "============================================================"
echo " Collecting per-sample QC"
echo "============================================================"
echo ""
echo " mosdepth summaries: ${MOSDEPTH_DIR}/*${MOSDEPTH_SUMMARY_SUFFIX}"
echo " Coverage summary:   ${COV_SUMMARY_TSV}"
echo " Sample panel:       ${PANEL_FILE}"
echo " Manifest:           ${MANIFEST}"
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
    $1 == "chrX"  { x_mean = $4 }
    $1 == "chrY"  { y_mean = $4 }
    $1 == "chrM"  { m_mean = $4 }
    END {
      auto_cov = (len > 0) ? sum / len : "NA"
      x_ratio  = (auto_cov > 0 && x_mean != "") ? x_mean / auto_cov : "NA"
      y_ratio  = (auto_cov > 0 && y_mean != "") ? y_mean / auto_cov : "NA"
      m_ratio  = (auto_cov > 0 && m_mean != "") ? m_mean / auto_cov : "NA"
      # Infer sex: males have X~0.5× autosome and detectable Y; females X~1.0×
      if (x_ratio != "NA" && y_ratio != "NA") {
        inferred_sex = (y_ratio > 0.1 && x_ratio < 0.75) ? "M" : "F"
      } else {
        inferred_sex = "NA"
      }
      printf "%.4f\t%.4f\t%.4f\t%s\t%s\n", auto_cov, x_ratio, y_ratio, inferred_sex, \
        (m_ratio != "NA") ? sprintf("%.4f", m_ratio) : "NA"
    }
  ' "${summary_file}"
}

# ── Helper: parse mosdepth global coverage distribution ─────────────────────
# The mosdepth .global.dist.txt file records the cumulative coverage
# distribution for the whole genome.  Format (tab-separated):
#   chrom/total   depth   fraction_at_>=_depth
# Example:
#   total   0   1.00
#   total   1   0.99
#   ...
# From this we derive:
#   - MEDIAN_GENOME_COV: depth at which cumulative fraction crosses 0.50
#   - PCT_GENOME_COV_10X: fraction of genome at >= 10× depth
#   - PCT_GENOME_COV_20X: fraction of genome at >= 20× depth
# Output: 3 tab-separated fields
parse_mosdepth_global_dist() {
  local dist_file="$1"
  awk 'BEGIN { FS="\t"; median="NA"; pct10x="NA"; pct20x="NA"; last_depth_above_half=-1 }
    $1 == "total" {
      depth = $2 + 0
      frac  = $3 + 0
      if (depth == 10) pct10x = sprintf("%.4f", frac * 100)
      if (depth == 20) pct20x = sprintf("%.4f", frac * 100)
      if (frac >= 0.5) last_depth_above_half = depth
    }
    END {
      if (last_depth_above_half >= 0) median = last_depth_above_half
      printf "%s\t%s\t%s\n", median, pct10x, pct20x
    }
  ' "${dist_file}"
}

# ── Build QC table ───────────────────────────────────────────────────────────
echo "Building QC table..."

DIST_NA="NA\tNA\tNA"
COV_SUMMARY_NA="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

HEADER="SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX\tMITO_COV_RATIO"
HEADER="${HEADER}\tMEDIAN_GENOME_COV\tPCT_GENOME_COV_10X\tPCT_GENOME_COV_20X"
HEADER="${HEADER}\tSD_COV\tMAD_COV\tIQR_COV\tMEDIAN_BIN_COV"
HEADER="${HEADER}\tHQ_MEDIAN_COV\tHQ_SD_COV\tHQ_MAD_COV\tHQ_IQR_COV\tMTDNA_CN"
echo -e "${HEADER}" > "${OUT_TSV}"

# Load coverage summary into an associative-style lookup (awk join later)
# Format: SAMPLE_ID  MEAN_COV  MEDIAN_COV  SD_COV  MAD_COV  IQR_COV
#         HQ_MEAN_COV  HQ_MEDIAN_COV  HQ_SD_COV  HQ_MAD_COV  HQ_IQR_COV
if [[ -f "${COV_SUMMARY_TSV}" ]]; then
  echo "  Loading coverage summary from: ${COV_SUMMARY_TSV}"
  COV_SUMMARY_LOADED=1
else
  echo "  WARNING: Coverage summary not found at ${COV_SUMMARY_TSV}."
  echo "  Run 03a_mosdepth_coverage_summary.sh first for SD/MAD/IQR/HQ metrics."
  echo "  Those columns will be NA."
  COV_SUMMARY_LOADED=0
fi

# Iterate over mosdepth summary files (one per sample)
FOUND=0
for summary in "${MOSDEPTH_DIR}"/*"${MOSDEPTH_SUMMARY_SUFFIX}"; do
  [[ -f "${summary}" ]] || continue

  # Derive sample ID from filename (strip suffix: .mosdepth.summary.txt)
  fname="$(basename "${summary}")"
  sample_id="${fname%%.*}"

  # Parse mosdepth summary
  mos_fields=$(parse_mosdepth_summary "${summary}")

  # Parse mosdepth global coverage distribution (if available)
  # The dist file has the same prefix as the summary but with .global.dist.txt suffix.
  dist_file="${summary/${MOSDEPTH_SUMMARY_SUFFIX}/${MOSDEPTH_GLOBAL_DIST_SUFFIX}}"
  if [[ -f "${dist_file}" ]]; then
    dist_fields=$(parse_mosdepth_global_dist "${dist_file}")
  else
    dist_fields="${DIST_NA}"
  fi

  # Look up coverage summary stats (SD, MAD, IQR, MEDIAN_BIN_COV, HQ stats)
  # Coverage summary columns: $1=SAMPLE_ID $2=MEAN_COV $3=MEDIAN_COV
  #   $4=SD_COV $5=MAD_COV $6=IQR_COV $7=HQ_MEAN_COV $8=HQ_MEDIAN_COV
  #   $9=HQ_SD_COV $10=HQ_MAD_COV $11=HQ_IQR_COV
  if [[ "${COV_SUMMARY_LOADED}" -eq 1 ]]; then
    cov_fields=$(awk -F'\t' -v sid="${sample_id}" '
      NR > 1 && $1 == sid {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $4, $5, $6, $3, $8, $9, $10, $11
        found = 1; exit
      }
      END { if (!found) printf "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" }
    ' "${COV_SUMMARY_TSV}")
  else
    cov_fields="${COV_SUMMARY_NA}"
  fi

  # Compute mtDNA CN = 2 × chrM_mean / HQ_autosomal_median
  # Combined mos_fields + cov_fields layout (13 tab-separated fields):
  #   $1=auto_cov  $2=x_ratio  $3=y_ratio  $4=inferred_sex  $5=m_ratio
  #   $6=SD_COV  $7=MAD_COV  $8=IQR_COV  $9=MEDIAN_BIN_COV  $10=HQ_MEDIAN_COV
  #   $11=HQ_SD_COV  $12=HQ_MAD_COV  $13=HQ_IQR_COV
  mtdna_cn=$(printf "%s\t%s" "${mos_fields}" "${cov_fields}" | awk -F'\t' '{
    auto_cov = $1 + 0; m_ratio = $5; hq_med = $10
    if (m_ratio != "NA" && auto_cov > 0 && hq_med != "NA" && hq_med + 0 > 0)
      printf "%.4f\n", 2 * m_ratio * auto_cov / (hq_med + 0)
    else
      print "NA"
  }')

  echo -e "${sample_id}\t${mos_fields}\t${dist_fields}\t${cov_fields}\t${mtdna_cn}" >> "${OUT_TSV}"
  FOUND=$(( FOUND + 1 ))
done

echo "  Processed ${FOUND} samples."

# ── Join with IGSR sample panel ──────────────────────────────────────────────
# The PED file (integrated_call_samples_v3.*.ALL.ped) has 12 tab-separated
# columns with header:
#   1  Family ID       6  Phenotype      11  Third Order
#   2  Individual ID   7  Population     12  Other Comments
#   3  Paternal ID     8  Relationship
#   4  Maternal ID     9  Siblings
#   5  Gender          10 Second Order
#
# NOTE: The PED file does NOT contain a SuperPopulation column.  We derive it
#       from the Population code ($7) using the standard 1000G Phase 3 mapping:
#         AFR: ACB ASW ESN GWD LWK MSL YRI
#         AMR: CLM MXL PEL PUR
#         EAS: CDX CHB CHS JPT KHV
#         EUR: CEU FIN GBR IBS TSI
#         SAS: BEB GIH ITU PJL STU
#
# Two relatedness annotations are produced — they are complementary:
#   FAMILY_ROLE  — the PED Relationship field ($8).  Describes each
#                  individual's declared role in the family structure
#                  (e.g. father, mother, child, unrel).
#   RELATEDNESS  — derived from paternal ($3) and maternal ($4) IDs.
#                  "unrelated" when both are "0" or empty (founders /
#                  unrelated individuals), "related" otherwise (children
#                  with known parents in the pedigree).
#
# These two fields are NOT fully concordant.  For example a "father" has
# RELATEDNESS = "unrelated" because founder parents are coded as 0/0 in the
# PED, while his child has RELATEDNESS = "related" because the parental IDs
# are filled in.  Conversely a small number of "unrel" individuals have non-
# zero parental IDs in the PED, giving RELATEDNESS = "related".
echo "Joining with sample panel metadata..."

if [[ ! -f "${PANEL_FILE}" ]]; then
  echo "  WARNING: Sample panel not found at ${PANEL_FILE}."
  echo "  Run 00_setup.sh to download it, or set PANEL_FILE."
  echo "  Skipping population/sex/relatedness annotation."
else
  PANEL_TSV="${QC_OUTPUT}/.panel_extracted.tmp"
  # Extract relevant columns from PED file (tab-separated, skip header if present)
  # Derive superpopulation from population code, and derive relatedness from
  # paternal_id ($3) and maternal_id ($4): both "0" or empty → unrelated.
  awk 'BEGIN {
    # 1000G Phase 3 population-to-superpopulation mapping (26 populations)
    sp["ACB"]="AFR"; sp["ASW"]="AFR"; sp["ESN"]="AFR"; sp["GWD"]="AFR"
    sp["LWK"]="AFR"; sp["MSL"]="AFR"; sp["YRI"]="AFR"
    sp["CLM"]="AMR"; sp["MXL"]="AMR"; sp["PEL"]="AMR"; sp["PUR"]="AMR"
    sp["CDX"]="EAS"; sp["CHB"]="EAS"; sp["CHS"]="EAS"; sp["JPT"]="EAS"; sp["KHV"]="EAS"
    sp["CEU"]="EUR"; sp["FIN"]="EUR"; sp["GBR"]="EUR"; sp["IBS"]="EUR"; sp["TSI"]="EUR"
    sp["BEB"]="SAS"; sp["GIH"]="SAS"; sp["ITU"]="SAS"; sp["PJL"]="SAS"; sp["STU"]="SAS"
  }
  NR > 1 {
    sex_label  = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    rel_label  = (($3 == "0" || $3 == "") && ($4 == "0" || $4 == "")) ? "unrelated" : "related"
    superpop   = ($7 in sp) ? sp[$7] : "NA"
    print $2 "\t" $7 "\t" superpop "\t" sex_label "\t" $8 "\t" rel_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}" 2>/dev/null || \
  # Fallback: some PED files lack a header line — process all rows
  awk 'BEGIN {
    sp["ACB"]="AFR"; sp["ASW"]="AFR"; sp["ESN"]="AFR"; sp["GWD"]="AFR"
    sp["LWK"]="AFR"; sp["MSL"]="AFR"; sp["YRI"]="AFR"
    sp["CLM"]="AMR"; sp["MXL"]="AMR"; sp["PEL"]="AMR"; sp["PUR"]="AMR"
    sp["CDX"]="EAS"; sp["CHB"]="EAS"; sp["CHS"]="EAS"; sp["JPT"]="EAS"; sp["KHV"]="EAS"
    sp["CEU"]="EUR"; sp["FIN"]="EUR"; sp["GBR"]="EUR"; sp["IBS"]="EUR"; sp["TSI"]="EUR"
    sp["BEB"]="SAS"; sp["GIH"]="SAS"; sp["ITU"]="SAS"; sp["PJL"]="SAS"; sp["STU"]="SAS"
  }
  {
    sex_label  = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    rel_label  = (($3 == "0" || $3 == "") && ($4 == "0" || $4 == "")) ? "unrelated" : "related"
    superpop   = ($7 in sp) ? sp[$7] : "NA"
    print $2 "\t" $7 "\t" superpop "\t" sex_label "\t" $8 "\t" rel_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}"

  # Join QC table with panel on SAMPLE_ID (column 1)
  MERGED_TSV="${QC_OUTPUT}/.merged.tmp"
  # Build merged header: current QC columns + panel columns
  CURRENT_HEADER="$(head -1 "${OUT_TSV}")"
  echo -e "${CURRENT_HEADER}\tPOPULATION\tSUPERPOPULATION\tREPORTED_SEX\tFAMILY_ROLE\tRELATEDNESS" \
    > "${MERGED_TSV}"

  # Use awk join: load panel into associative array, then annotate QC rows
  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == FNR {
      pop[$1] = $2; super[$1] = $3; sex[$1] = $4; role[$1] = $5; rel[$1] = $6
      next
    }
    FNR == 1 { next }   # skip QC header
    {
      sid = $1
      p  = (sid in pop)   ? pop[sid]   : "NA"
      s  = (sid in super) ? super[sid] : "NA"
      sx = (sid in sex)   ? sex[sid]   : "NA"
      ro = (sid in role)  ? role[sid]  : "NA"
      r  = (sid in rel)   ? rel[sid]   : "NA"
      print $0, p, s, sx, ro, r
    }
  ' "${PANEL_TSV}" "${OUT_TSV}" >> "${MERGED_TSV}"

  mv "${MERGED_TSV}" "${OUT_TSV}"
  rm -f "${PANEL_TSV}"
  echo "  Population/sex/relatedness metadata joined."
fi

# ── Join with manifest batch annotations ─────────────────────────────────────
# The manifest (built by 00_setup.sh) contains batch-level metadata parsed from
# the sequence.index files: RELEASE_BATCH, CENTER_NAME, STUDY_ID,
# INSTRUMENT_MODEL, LIBRARY_NAME.
echo "Joining with manifest batch annotations..."

if [[ ! -f "${MANIFEST}" ]]; then
  echo "  WARNING: Manifest not found at ${MANIFEST}."
  echo "  Run 00_setup.sh to generate it, or set MANIFEST."
  echo "  Skipping batch annotations."
else
  MANIFEST_TSV="${QC_OUTPUT}/.manifest_extracted.tmp"
  # Extract batch columns from manifest (skip header)
  # Manifest columns: SAMPLE_ID(1) CRAM_FTP_URL(2) CRAI_FTP_URL(3) CRAM_MD5(4)
  #   RELEASE_BATCH(5) CENTER_NAME(6) STUDY_ID(7) INSTRUMENT_MODEL(8) LIBRARY_NAME(9)
  awk -F'\t' 'NR > 1 && $1 != "" {
    print $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9
  }' OFS='\t' "${MANIFEST}" > "${MANIFEST_TSV}"

  # Determine current header and append new column names
  CURRENT_HEADER="$(head -1 "${OUT_TSV}")"
  MERGED_TSV="${QC_OUTPUT}/.merged_batch.tmp"
  echo -e "${CURRENT_HEADER}\tRELEASE_BATCH\tCENTER_NAME\tSTUDY_ID\tINSTRUMENT_MODEL\tLIBRARY_NAME" \
    > "${MERGED_TSV}"

  # Use awk join: load manifest batch data, then annotate QC rows
  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == FNR {
      batch[$1]  = $2; center[$1] = $3; study[$1] = $4
      inst[$1]   = $5; lib[$1]    = $6
      next
    }
    FNR == 1 { next }   # skip QC header
    {
      sid = $1
      b  = (sid in batch)  ? batch[sid]  : "NA"
      c  = (sid in center) ? center[sid] : "NA"
      st = (sid in study)  ? study[sid]  : "NA"
      i  = (sid in inst)   ? inst[sid]   : "NA"
      l  = (sid in lib)    ? lib[sid]    : "NA"
      print $0, b, c, st, i, l
    }
  ' "${MANIFEST_TSV}" "${OUT_TSV}" >> "${MERGED_TSV}"

  mv "${MERGED_TSV}" "${OUT_TSV}"
  rm -f "${MANIFEST_TSV}"
  echo "  Batch annotations joined."
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
echo " To join with NGS-PCA results (use a literal tab for -t):"
echo "   join -t \$'\\t' -1 1 -2 1 \\"
echo "     <(sort ${OUT_TSV}) \\"
echo "     <(sort ${NGSPCA_OUTPUT}/svd.pcs.txt) \\"
echo "     > ${QC_OUTPUT}/pcs_with_qc.tsv"
echo "============================================================"
