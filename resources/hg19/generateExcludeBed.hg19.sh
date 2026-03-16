#!/usr/bin/env bash
# generateExcludeBed.hg19.sh — Build the NGS-PCA exclusion BED for hg19/GRCh37.
#
# This script produces a merged exclusion BED that combines:
#   1. Low-mappability regions (from a GEM mappability track, kmer 50)
#   2. Structural-variant blacklist (10x Genomics sv_blacklist.bed for hg19)
#   3. Database of Genomic Variants (DGV)
#   4. Genomic SuperDups (UCSC)
#
# Prerequisites (must be on PATH):
#   awk, cut, sort, bedtools, bgzip, tabix
#
# Required input files (download instructions below each):
#   hg19_canonical_50.bed
#       — Mappability BED produced by GEM tools for hg19, kmer 50
#         (https://github.com/gemtools/gemtools)
#   sv_blacklist.bed
#       — wget "http://cf.10xgenomics.com/supp/genome/hg19/sv_blacklist.bed"
#   GRCh37_hg19_variants_2016-05-15.txt
#       — wget "http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt"
#   genomicSuperDups.hg19.bed
#       — Download from UCSC Table Browser > All Tables > genomicSuperDups > get BED
#
# Usage:
#   bash generateExcludeBed.hg19.sh
#
# Output:
#   ngs_pca_exclude.hg19.sv_blacklist.map.kmer.<kmer>.<threshold>.dgv.gsd.sorted.merge.bed.gz
#   ngs_pca_exclude.hg19.sv_blacklist.map.kmer.<kmer>.<threshold>.dgv.gsd.sorted.merge.bed.gz.tbi

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Mappability threshold: regions with score < threshold are excluded (1.0 = non-unique)
threshold=1
kmer=50

# Directory where finished BED files will be copied
outDir="${HOME}/git/NGS-PCA/resources/hg19"

# ---------------------------------------------------------------------------
# Build exclusion BED
# ---------------------------------------------------------------------------

prefix="ngs_pca_exclude.hg19.sv_blacklist.map.kmer.${kmer}.${threshold}"

# ------------------------------------------------------------------
# Step 1: Low-mappability regions
# ------------------------------------------------------------------
echo "[1/4] Extracting low-mappability regions..."
awk -v threshold="$threshold" '$5 < threshold {print $1"\t"$2"\t"$3}' \
    "hg19_canonical_${kmer}.bed" \
    | bedtools merge \
    | bgzip > "hg19_canonical_${kmer}.${threshold}.merge.bed.gz"
tabix "hg19_canonical_${kmer}.${threshold}.merge.bed.gz"

# ------------------------------------------------------------------
# Step 2: Merge with SV blacklist
# ------------------------------------------------------------------
# sv_blacklist.bed: wget "http://cf.10xgenomics.com/supp/genome/hg19/sv_blacklist.bed"
echo "[2/4] Merging SV blacklist..."
cut -f 1-3 sv_blacklist.bed > "${prefix}.bed"
zcat "hg19_canonical_${kmer}.${threshold}.merge.bed.gz" >> "${prefix}.bed"
bedtools sort -i "${prefix}.bed" \
    | bedtools merge > "${prefix}.sorted.merge.bed"

# ------------------------------------------------------------------
# Step 3: Add DGV structural variants
# ------------------------------------------------------------------
# DGV: wget "http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt"
echo "[3/4] Adding DGV variants..."
awk '{
    line = "chr"$2"\t"$3"\t"$4
    if (line !~ /start/) print line
  }' GRCh37_hg19_variants_2016-05-15.txt > "${prefix}.dgv.bed"
cat "${prefix}.sorted.merge.bed" >> "${prefix}.dgv.bed"
bedtools sort -i "${prefix}.dgv.bed" \
    | bedtools merge > "${prefix}.dgv.sorted.merge.bed"

# ------------------------------------------------------------------
# Step 4: Add genomic SuperDups
# ------------------------------------------------------------------
# genomicSuperDups: UCSC Table Browser > All Tables > genomicSuperDups > get BED
echo "[4/4] Adding genomicSuperDups..."
cp "${prefix}.dgv.sorted.merge.bed" "${prefix}.dgv.gsd.bed"
cut -f 1-3 genomicSuperDups.hg19.bed >> "${prefix}.dgv.gsd.bed"
bedtools sort -i "${prefix}.dgv.gsd.bed" \
    | bedtools merge \
    | bgzip > "${prefix}.dgv.gsd.sorted.merge.bed.gz"

# ------------------------------------------------------------------
# Copy to repository resources directory and index
# ------------------------------------------------------------------
mkdir -p "$outDir"
cp "${prefix}.dgv.gsd.sorted.merge.bed.gz" \
    "${outDir}/${prefix}.dgv.gsd.sorted.merge.bed.gz"
tabix "${outDir}/${prefix}.dgv.gsd.sorted.merge.bed.gz"

echo "Done. Output: ${outDir}/${prefix}.dgv.gsd.sorted.merge.bed.gz"



