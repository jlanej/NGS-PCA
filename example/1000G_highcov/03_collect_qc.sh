#!/usr/bin/env bash
# =============================================================================
# 03_collect_qc.sh — Aggregate per-sample QC into a single overlay-ready table
# =============================================================================
#
# Collects four sources of publicly available or pipeline-generated QC:
#
#   1. mosdepth summary (.mosdepth.summary.txt) — generated during step 01.
#      Provides per-chromosome mean coverage.  From this we derive:
#        - Mean autosomal coverage (chr1–22)
#        - X and Y coverage ratios → inferred genetic sex
#
#   2. BAS files (.bam.bas) — per-readgroup alignment statistics from the 1000G
#      alignment pipeline.  Automatically downloaded from the NYGC alignment
#      index files on the EBI FTP.  Provides all standard BAS QC metrics
#      (see https://www.internationalgenome.org/category/bas/):
#        - Total bases, mapped bases, total reads, mapped reads
#        - Mapped reads paired in sequencing, mapped reads properly paired
#        - % mismatched bases, average quality of mapped bases
#        - Mean/median insert size and variability (SD, MAD)
#        - Duplicate reads and duplicate bases
#        - Derived: % mapped reads, % duplicate reads
#
#   3. IGSR sample panel (integrated_call_samples_v3.*.ALL.ped) — downloaded
#      during setup.  Provides:
#        - Population (e.g. GBR, YRI, CHB)
#        - Superpopulation (AFR, AMR, EAS, EUR, SAS)
#        - Reported sex
#        - Relatedness (unrelated vs related, from paternal/maternal IDs)
#
#   4. Sample manifest (manifest.tsv) — built during setup from the two NYGC
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

mkdir -p "${QC_OUTPUT}" "${BAS_DIR}"

OUT_TSV="${QC_OUTPUT}/sample_qc.tsv"
MOSDEPTH_SUMMARY_SUFFIX=".mosdepth.summary.txt"

echo "============================================================"
echo " Collecting per-sample QC"
echo "============================================================"
echo ""
echo " mosdepth summaries: ${MOSDEPTH_DIR}/*${MOSDEPTH_SUMMARY_SUFFIX}"
echo " BAS files:          ${BAS_DIR}/*.bam.bas"
echo " Sample panel:       ${PANEL_FILE}"
echo " Manifest:           ${MANIFEST}"
echo " Output:             ${OUT_TSV}"
echo ""

# ── Download BAS files from alignment index ─────────────────────────────────
# The NYGC 30x alignment index files list BAS file FTP paths in column 5.
# Format: CRAM  CRAM_MD5  CRAI  CRAI_MD5  BAS  BAS_MD5
# We download each BAS file and name it ${SAMPLE_ID}.bam.bas in $BAS_DIR.
download_bas_files() {
  echo "Checking for BAS files to download..."

  local alignment_indices=()
  for idx_url_var in ALIGNMENT_INDEX_URL_2504 ALIGNMENT_INDEX_URL_698; do
    local idx_url="${!idx_url_var}"
    [[ -z "${idx_url}" ]] && continue
    local idx_file_var="${idx_url_var/URL/FILE}"
    local idx_file="${!idx_file_var}"
    [[ -z "${idx_file}" ]] && continue

    if [[ ! -s "${idx_file}" ]]; then
      echo "  Downloading alignment index: ${idx_url}"
      if curl -sSL -o "${idx_file}" "${idx_url}" 2>/dev/null || \
         wget -q -O "${idx_file}" "${idx_url}" 2>/dev/null; then
        echo "    Downloaded: ${idx_file}"
      else
        echo "    WARNING: Could not download alignment index from ${idx_url}"
        echo "    BAS files for this batch will be skipped."
        continue
      fi
    else
      echo "  Alignment index already present: ${idx_file}"
    fi
    [[ -s "${idx_file}" ]] && alignment_indices+=("${idx_file}")
  done

  if [[ ${#alignment_indices[@]} -eq 0 ]]; then
    echo "  No alignment indexes available. BAS columns will be NA."
    return
  fi

  # Build a lookup: SAMPLE_ID → BAS_FTP_URL from the alignment index
  local bas_url_file="${QC_OUTPUT}/.bas_urls.tmp"
  : > "${bas_url_file}"
  for idx_file in "${alignment_indices[@]}"; do
    awk -F'\t' '
      /^#/ { next }
      $5 != "" {
        # Extract sample ID from the BAS filename (column 5)
        # Typical: ftp:/ftp.../data/POP/SAMPLE/.../SAMPLE.*.bam.bas
        n = split($5, parts, "/")
        fname = parts[n]
        # Strip everything after the first "." to get sample ID
        split(fname, namep, ".")
        sample = namep[1]
        # Fix single-slash ftp:/ to ftp:// (common in IGSR alignment indexes)
        bas_url = $5
        sub(/^ftp:\/ftp\./, "ftp://ftp.", bas_url)
        if (sample != "" && bas_url != "")
          print sample "\t" bas_url
      }
    ' "${idx_file}"
  done >> "${bas_url_file}"

  local total_bas
  total_bas=$(wc -l < "${bas_url_file}")
  echo "  Found ${total_bas} BAS file URLs in alignment index(es)."

  # Download BAS files that are not already present
  local downloaded=0
  local skipped=0
  local failed=0
  while IFS=$'\t' read -r sample_id bas_url; do
    local dest="${BAS_DIR}/${sample_id}.bam.bas"
    if [[ -s "${dest}" ]]; then
      skipped=$(( skipped + 1 ))
      continue
    fi
    if wget -q -O "${dest}" "${bas_url}" 2>/dev/null || \
       curl -sSL -o "${dest}" "${bas_url}" 2>/dev/null; then
      if [[ -s "${dest}" ]]; then
        downloaded=$(( downloaded + 1 ))
      else
        rm -f "${dest}"
        failed=$(( failed + 1 ))
      fi
    else
      rm -f "${dest}"
      failed=$(( failed + 1 ))
    fi
  done < "${bas_url_file}"

  rm -f "${bas_url_file}"
  echo "  BAS download complete: ${downloaded} downloaded, ${skipped} already present, ${failed} failed."
}

download_bas_files

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
# See: https://www.internationalgenome.org/category/bas/
#
# Standard 21-column format (one row per readgroup):
#   1: bam_filename           8: #_total_bases           15: average_quality_of_mapped_bases
#   2: md5                    9: #_mapped_bases           16: mean_insert_size
#   3: study                 10: #_total_reads            17: insert_size_sd
#   4: sample                11: #_mapped_reads           18: median_insert_size
#   5: platform              12: #_mapped_reads_paired_in_sequencing
#   6: library               13: #_mapped_reads_properly_paired
#   7: readgroup             14: %_of_mismatched_bases    19: insert_size_median_absolute_deviation
#                                                         20: #_duplicate_reads
#                                                         21: #_duplicate_bases
#
# When multiple readgroups exist, count columns (8-13, 20-21) are summed.
# Percentage/average columns (14-15) are recomputed from aggregated counts.
# Insert size columns (16-19) are averaged across readgroups (weighted).
#
# Output: 16 tab-separated fields:
#   TOTAL_BASES  MAPPED_BASES  TOTAL_READS  MAPPED_READS
#   MAPPED_READS_PAIRED  MAPPED_READS_PROPERLY_PAIRED
#   PCT_MISMATCHED_BASES  AVG_QUALITY_MAPPED_BASES
#   MEAN_INSERT_SIZE  INSERT_SIZE_SD  MEDIAN_INSERT_SIZE  INSERT_SIZE_MAD
#   DUPLICATE_READS  DUPLICATE_BASES  PCT_MAPPED  PCT_DUPLICATE
parse_bas_file() {
  local bas_file="$1"
  awk 'BEGIN { FS="\t"; OFS="\t"
    # Positional defaults (standard 21-column BAS format)
    total_bases_col=8; mapped_bases_col=9; total_reads_col=10; mapped_reads_col=11
    paired_col=12; proper_col=13; mismatch_col=14; avgqual_col=15
    ins_mean_col=16; ins_sd_col=17; ins_med_col=18; ins_mad_col=19
    dup_reads_col=20; dup_bases_col=21
    has_header=0
  }
  NR == 1 {
    # Detect header row: if first field contains letters and is not a file path
    if ($1 ~ /[a-zA-Z]/ && $1 !~ /\.bam$/ && $1 !~ /\.cram$/ && $1 !~ /^[\/.]/) {
      has_header = 1
      # Find column indices by header name (case-insensitive)
      for (i = 1; i <= NF; i++) {
        h = tolower($i)
        gsub(/^#/, "", h)  # strip leading #
        if (h == "_total_bases" || h == "total_bases")           total_bases_col = i
        if (h == "_mapped_bases" || h == "mapped_bases")         mapped_bases_col = i
        if (h == "_total_reads" || h == "total_reads")           total_reads_col = i
        if (h == "_mapped_reads" || h == "mapped_reads")         mapped_reads_col = i
        if (h == "_mapped_reads_paired_in_sequencing" || h == "mapped_reads_paired_in_sequencing")   paired_col = i
        if (h == "_mapped_reads_properly_paired" || h == "mapped_reads_properly_paired")             proper_col = i
        if (h == "%_of_mismatched_bases" || h == "of_mismatched_bases" || h == "mismatched_bases")   mismatch_col = i
        if (h == "average_quality_of_mapped_bases" || h == "avg_quality_mapped_bases")               avgqual_col = i
        if (h == "mean_insert_size")                             ins_mean_col = i
        if (h == "insert_size_sd")                               ins_sd_col = i
        if (h == "median_insert_size")                           ins_med_col = i
        if (h == "insert_size_median_absolute_deviation" || h == "insert_size_mad") ins_mad_col = i
        if (h == "_duplicate_reads" || h == "duplicate_reads")   dup_reads_col = i
        if (h == "_duplicate_bases" || h == "duplicate_bases")   dup_bases_col = i
      }
      next
    }
  }
  {
    nrg++
    # Sum count columns across readgroups
    s_total_bases   += ($total_bases_col   + 0)
    s_mapped_bases  += ($mapped_bases_col  + 0)
    s_total_reads   += ($total_reads_col   + 0)
    s_mapped_reads  += ($mapped_reads_col  + 0)
    s_paired        += ($paired_col        + 0)
    s_proper        += ($proper_col        + 0)
    s_dup_reads     += ($dup_reads_col     + 0)
    s_dup_bases     += ($dup_bases_col     + 0)
    # For mismatch %, average quality, and insert size stats: accumulate for averaging
    if ($mismatch_col + 0 > 0 || $mismatch_col == "0")  { s_mismatch += ($mismatch_col + 0); n_mismatch++ }
    if ($avgqual_col + 0 > 0 || $avgqual_col == "0")    { s_avgqual  += ($avgqual_col  + 0); n_avgqual++  }
    if ($ins_mean_col + 0 > 0 || $ins_mean_col == "0")  { s_ins_mean += ($ins_mean_col + 0); n_ins_mean++ }
    if ($ins_sd_col + 0 > 0 || $ins_sd_col == "0")      { s_ins_sd   += ($ins_sd_col   + 0); n_ins_sd++   }
    if ($ins_med_col + 0 > 0 || $ins_med_col == "0")    { s_ins_med  += ($ins_med_col  + 0); n_ins_med++  }
    if ($ins_mad_col + 0 > 0 || $ins_mad_col == "0")    { s_ins_mad  += ($ins_mad_col  + 0); n_ins_mad++  }
  }
  END {
    if (nrg == 0) {
      # No data rows — output all NAs (16 fields)
      printf "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
      exit
    }
    # Derived percentages
    pct_mapped = "NA"
    if (s_total_bases > 0 && s_mapped_bases >= 0) {
      pct_mapped = sprintf("%.4f", s_mapped_bases / s_total_bases * 100)
    }
    pct_dup = "NA"
    if (s_mapped_bases > 0 && s_dup_bases >= 0) {
      pct_dup = sprintf("%.4f", s_dup_bases / s_mapped_bases * 100)
    }
    # Averaged metrics
    mismatch  = (n_mismatch > 0) ? sprintf("%.2f", s_mismatch / n_mismatch) : "NA"
    avgqual   = (n_avgqual  > 0) ? sprintf("%.2f", s_avgqual  / n_avgqual)  : "NA"
    ins_mean  = (n_ins_mean > 0) ? sprintf("%.0f", s_ins_mean / n_ins_mean) : "NA"
    ins_sd    = (n_ins_sd   > 0) ? sprintf("%.2f", s_ins_sd   / n_ins_sd)   : "NA"
    ins_med   = (n_ins_med  > 0) ? sprintf("%.0f", s_ins_med  / n_ins_med) : "NA"
    ins_mad   = (n_ins_mad  > 0) ? sprintf("%.0f", s_ins_mad  / n_ins_mad) : "NA"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
      s_total_bases, s_mapped_bases, s_total_reads, s_mapped_reads, \
      s_paired, s_proper, \
      mismatch, avgqual, \
      ins_mean, ins_sd, ins_med, ins_mad, \
      s_dup_reads, s_dup_bases, \
      pct_mapped, pct_dup
  }
  ' "${bas_file}"
}

# ── Build QC table ───────────────────────────────────────────────────────────
echo "Building QC table..."

# BAS-derived columns (16 fields from parse_bas_file):
#   TOTAL_BASES  MAPPED_BASES  TOTAL_READS  MAPPED_READS
#   MAPPED_READS_PAIRED  MAPPED_READS_PROPERLY_PAIRED
#   PCT_MISMATCHED_BASES  AVG_QUALITY_MAPPED_BASES
#   MEAN_INSERT_SIZE  INSERT_SIZE_SD  MEDIAN_INSERT_SIZE  INSERT_SIZE_MAD
#   DUPLICATE_READS  DUPLICATE_BASES  PCT_MAPPED  PCT_DUPLICATE
BAS_NA="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

HEADER="SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX"
HEADER="${HEADER}\tTOTAL_BASES\tMAPPED_BASES\tTOTAL_READS\tMAPPED_READS"
HEADER="${HEADER}\tMAPPED_READS_PAIRED\tMAPPED_READS_PROPERLY_PAIRED"
HEADER="${HEADER}\tPCT_MISMATCHED_BASES\tAVG_QUALITY_MAPPED_BASES"
HEADER="${HEADER}\tMEAN_INSERT_SIZE\tINSERT_SIZE_SD\tMEDIAN_INSERT_SIZE\tINSERT_SIZE_MAD"
HEADER="${HEADER}\tDUPLICATE_READS\tDUPLICATE_BASES\tPCT_MAPPED\tPCT_DUPLICATE"
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
    bas_fields="${BAS_NA}"
    MISSING_BAS=$(( MISSING_BAS + 1 ))
  fi

  echo -e "${sample_id}\t${mos_fields}\t${bas_fields}" >> "${OUT_TSV}"
  FOUND=$(( FOUND + 1 ))
done

echo "  Processed ${FOUND} samples (${MISSING_BAS} missing BAS files)."

# ── Join with IGSR sample panel ──────────────────────────────────────────────
# The PED file columns: family_id sample_id paternal maternal sex phenotype pop superpop
# We extract: sample_id, pop, superpop, reported_sex (1=male,2=female),
# and relatedness (paternal_id and maternal_id both "0" → unrelated).
echo "Joining with sample panel metadata..."

if [[ ! -f "${PANEL_FILE}" ]]; then
  echo "  WARNING: Sample panel not found at ${PANEL_FILE}."
  echo "  Run 00_setup.sh to download it, or set PANEL_FILE."
  echo "  Skipping population/sex/relatedness annotation."
else
  PANEL_TSV="${QC_OUTPUT}/.panel_extracted.tmp"
  # Extract relevant columns from PED file (tab-separated, skip header if present)
  # Also derive relatedness: if paternal_id ($3) and maternal_id ($4) are both
  # "0" or empty → unrelated (standard PED convention for unknown parents).
  awk 'NR > 1 {
    sex_label = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    rel_label = (($3 == "0" || $3 == "") && ($4 == "0" || $4 == "")) ? "unrelated" : "related"
    print $2 "\t" $7 "\t" $8 "\t" sex_label "\t" rel_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}" 2>/dev/null || \
  # Fallback: some PED files lack a header line — process all rows
  awk '{
    sex_label = ($5 == 1) ? "M" : ($5 == 2) ? "F" : "NA"
    rel_label = (($3 == "0" || $3 == "") && ($4 == "0" || $4 == "")) ? "unrelated" : "related"
    print $2 "\t" $7 "\t" $8 "\t" sex_label "\t" rel_label
  }' OFS='\t' "${PANEL_FILE}" > "${PANEL_TSV}"

  # Join QC table with panel on SAMPLE_ID (column 1)
  MERGED_TSV="${QC_OUTPUT}/.merged.tmp"
  # Build merged header: current QC columns + panel columns
  CURRENT_HEADER="$(head -1 "${OUT_TSV}")"
  echo -e "${CURRENT_HEADER}\tPOPULATION\tSUPERPOPULATION\tREPORTED_SEX\tRELATEDNESS" \
    > "${MERGED_TSV}"

  # Use awk join: load panel into associative array, then annotate QC rows
  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == FNR {
      pop[$1] = $2; super[$1] = $3; sex[$1] = $4; rel[$1] = $5
      next
    }
    FNR == 1 { next }   # skip QC header
    {
      sid = $1
      p  = (sid in pop)   ? pop[sid]   : "NA"
      s  = (sid in super) ? super[sid] : "NA"
      sx = (sid in sex)   ? sex[sid]   : "NA"
      r  = (sid in rel)   ? rel[sid]   : "NA"
      print $0, p, s, sx, r
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
