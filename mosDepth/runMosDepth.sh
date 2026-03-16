#!/usr/bin/env bash
# runMosDepth.sh — Generate per-sample mosdepth commands and optionally run them.
#
# Usage:
#   bash runMosDepth.sh
#
# Prerequisites:
#   mosdepth (https://github.com/brentp/mosdepth)
#   GNU parallel (optional, for parallel execution)
#
# Installation (conda/bioconda — Python 3):
#   conda config --add channels defaults
#   conda config --add channels conda-forge
#   conda config --add channels bioconda
#   conda install mosdepth
#
# See also:
#   https://docs.anaconda.com/anaconda/install/linux
#   https://github.com/brentp/mosdepth

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration — edit these variables before running
# ---------------------------------------------------------------------------

# Directory where mosdepth output files will be written
mosDepthResultsDir="${HOME}/mosdepthOutput/"

# Reference genome FASTA (GRCh38 recommended)
ref="${HOME}/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Directory containing input CRAM (or BAM) files
inDir="${HOME}/TOPMedCrams/"

# File extension of input files ("cram" or "bam")
ext="cram"

# File that will hold one mosdepth command per line
cmdFile="${mosDepthResultsDir}mosdepthCommands.txt"

# ---------------------------------------------------------------------------
# Step 1: Create output directory and command file
# ---------------------------------------------------------------------------

mkdir -p "$mosDepthResultsDir"

# Truncate (or create) the command file
: > "$cmdFile"

# ---------------------------------------------------------------------------
# Step 2: Generate one mosdepth command per input file
# ---------------------------------------------------------------------------
# Options used:
#   -n          skip per-base depth output (saves disk space and time)
#   -t 2        use 2 threads per sample
#   --by 1000   compute coverage in 1 000 bp bins
#   --fasta     required for CRAM decoding

for file in "$inDir"*."$ext"; do
    [ -e "$file" ] || { echo "No ${ext} files found in ${inDir}"; exit 1; }
    out=$(basename "$file" ".$ext")
    echo "mosdepth -n -t 2 --by 1000 --fasta $ref ${mosDepthResultsDir}${out}.mos $file" >> "$cmdFile"
done

echo "Commands written to: $cmdFile"

# ---------------------------------------------------------------------------
# Step 3: Run the commands (uncomment one of the options below)
# ---------------------------------------------------------------------------

# Option A — sequential (simple, no dependencies)
# bash "$cmdFile"

# Option B — parallel with GNU parallel (recommended for large cohorts)
#   "parallel --jobs 12 < $cmdFile" uses 12 samples × 2 threads = 24 threads total
# parallel --jobs 12 < "$cmdFile"

