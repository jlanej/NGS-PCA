#!/usr/bin/env bash
# generateExcludeBed.sh — Build the NGS-PCA exclusion BED for GRCh38.
#
# This script produces a merged BED file that combines:
#   1. Low-mappability regions (from GEM mappability tracks, kmer-specific)
#   2. Structural-variant blacklist (10x Genomics sv_blacklist.bed)
#   3. Database of Genomic Variants (DGV)
#   4. Genomic SuperDups (UCSC)
#
# Prerequisites (must be on PATH):
#   awk, cut, sort, bedtools, gzip
#
# Required input files (download instructions below each):
#   GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_<kmer>.bed
#       — mappability BED produced by GEM tools (one file per kmer length)
#       — https://github.com/gemtools/gemtools
#   sv_blacklist.bed
#       — wget http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed
#   GRCh38_hg38_variants_2016-08-31.txt
#       — wget "http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2016-08-31.txt"
#   genomicSuperDups.hg38.bed
#       — Download from UCSC Table Browser > All Tables > genomicSuperDups > get BED
#
# Usage:
#   bash generateExcludeBed.sh
#
# Output (for each kmer in the loop):
#   ngs_pca_exclude.sv_blacklist.map.kmer.<kmer>.<threshold>.dgv.gsd.sorted.merge.bed.gz

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Mappability threshold: regions with score < threshold are excluded (1.0 = non-unique)
threshold=1.0

# Directory where finished BED files will be copied
outDir="${HOME}/git/NGS-PCA/resources/GRCh38"

# ---------------------------------------------------------------------------
# Build exclusion BEDs for each kmer length
# ---------------------------------------------------------------------------

for kmer in 100 125 150; do

    echo "=== Processing kmer=${kmer} threshold=${threshold} ==="

    mappa_input="GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_${kmer}.bed"
    mappa_thresh="GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_${kmer}.${threshold}.bed"
    mappa_merged="GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_${kmer}.${threshold}.merge.bed"

    prefix="ngs_pca_exclude.sv_blacklist.map.kmer.${kmer}.${threshold}"

    # ------------------------------------------------------------------
    # Step 1: Low-mappability regions
    # ------------------------------------------------------------------
    echo "  [1/4] Extracting low-mappability regions..."
    awk -v threshold="$threshold" '$5 < threshold {print $1"\t"$2"\t"$3}' \
        "$mappa_input" > "$mappa_thresh"
    bedtools merge -i "$mappa_thresh" > "$mappa_merged"

    # ------------------------------------------------------------------
    # Step 2: Merge with SV blacklist
    # ------------------------------------------------------------------
    # sv_blacklist.bed: wget http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed
    echo "  [2/4] Merging SV blacklist..."
    cut -f 1-3 sv_blacklist.bed > "${prefix}.bed"
    cat "$mappa_merged" >> "${prefix}.bed"
    sort -k1,1 -k2,2n "${prefix}.bed" \
        | bedtools merge > "${prefix}.sorted.merge.bed"

    # ------------------------------------------------------------------
    # Step 3: Add DGV structural variants
    # ------------------------------------------------------------------
    # DGV: wget "http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2016-08-31.txt"
    echo "  [3/4] Adding DGV variants..."
    awk '{
        line = "chr"$2"\t"$3"\t"$4
        if (line !~ /start/) print line
      }' GRCh38_hg38_variants_2016-08-31.txt > "${prefix}.dgv.bed"
    cat "${prefix}.sorted.merge.bed" >> "${prefix}.dgv.bed"
    sort -k1,1 -k2,2n "${prefix}.dgv.bed" \
        | bedtools merge > "${prefix}.dgv.sorted.merge.bed"

    # ------------------------------------------------------------------
    # Step 4: Add genomic SuperDups
    # ------------------------------------------------------------------
    # genomicSuperDups: UCSC Table Browser > All Tables > genomicSuperDups > get BED
    echo "  [4/4] Adding genomicSuperDups..."
    cp "${prefix}.dgv.sorted.merge.bed" "${prefix}.dgv.gsd.bed"
    cut -f 1-3 genomicSuperDups.hg38.bed >> "${prefix}.dgv.gsd.bed"
    sort -k1,1 -k2,2n "${prefix}.dgv.gsd.bed" \
        | bedtools merge > "${prefix}.dgv.gsd.sorted.merge.bed"

    # ------------------------------------------------------------------
    # Copy and compress the final BED
    # ------------------------------------------------------------------
    mkdir -p "$outDir"
    gzip -c "${prefix}.dgv.gsd.sorted.merge.bed" \
        > "${outDir}/${prefix}.dgv.gsd.sorted.merge.bed.gz"
    echo "  Written: ${outDir}/${prefix}.dgv.gsd.sorted.merge.bed.gz"

done

echo "Done."
