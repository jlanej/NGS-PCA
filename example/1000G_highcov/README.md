# 1000 Genomes 30x High-Coverage: End-to-End NGS-PCA Example

Fully reproducible pipeline for computing ~200 coverage-based principal components from the complete set of **3,202** [1000 Genomes 30x high-coverage WGS samples](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) (NYGC, GRCh38).

> **Cite:** If you use these results in a publication, please cite the NGS-PCA repository and the 1000 Genomes high-coverage data:
>
> - Byrska-Bishop, M. et al. (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell*, 185(18), 3426–3440.e19. https://doi.org/10.1016/j.cell.2022.08.004

---

## Overview

The pipeline has four stages, each implemented as a standalone script that can be submitted to a SLURM-based HPC scheduler:

> **SLURM note:** Submit jobs from this directory (`example/1000G_highcov`) so the scripts can source `config.sh` via `$SLURM_SUBMIT_DIR`.

| Stage | Script | What it does | SLURM type |
|-------|--------|-------------|------------|
| **0** | `00_setup.sh` | Pull container image, download reference genome, build sample manifest, download sample panel | Interactive / login node |
| **1** | `01_download_and_mosdepth.sh` | For each sample: download CRAM via Aspera/wget → run mosdepth → remove CRAM | Array job (3,202 tasks) |
| **2** | `02_run_ngspca.sh` | Run NGS-PCA on all mosdepth results → ~200 PCs | Single large-memory job |
| **3a** | `03a_mosdepth_coverage_summary.sh` | Compute autosomal coverage stats (mean, median, SD, MAD, IQR) and HQ statistics (non-excluded bins) from mosdepth output | Parallelized (all cores) |
| **3** | `03_collect_qc.sh` | Aggregate per-sample QC into one table for batch-effect overlay | Interactive or short job |

All three stages use the **same container image** (`ghcr.io/jlanej/ngs-pca:latest`), which bundles:

- **NGS-PCA** — randomized SVD for coverage PCA
- **mosdepth** v0.3.9 — fast BAM/CRAM depth calculation
- Pre-built **exclusion BED files** for GRCh38

No additional software installation is required beyond [Apptainer](https://apptainer.org/) (Singularity). IBM Aspera Connect (`ascp`) is an optional system-level tool for faster downloads — the pipeline automatically falls back to `wget` if `ascp` is not available.

---

## Prerequisites

| Requirement | Notes |
|---|---|
| **Apptainer ≥ 1.0** | Most HPC clusters provide this via `module load apptainer` or `module load singularity` |
| **SLURM** | Adjust `#SBATCH` directives if using PBS/SGE/LSF |
| **Internet access** | Needed on the node running `00_setup.sh` and the compute nodes running `01_download_and_mosdepth.sh` |
| **Disk space** | ~50 GB for mosdepth outputs (kept), ~3 GB for reference, ~25 GB peak per CRAM (cleaned up) |
| **Memory** | Step 01: 4 GB/task; Step 02: ~256 GB for 3,202 samples |

---

## Quick start

```bash
# Clone the repository
git clone https://github.com/jlanej/NGS-PCA.git
cd NGS-PCA/example/1000G_highcov

# 1. Edit config.sh — set WORK_DIR to a scratch/project directory
export WORK_DIR=/scratch/$USER/1000G_highcov

# 2. Run one-time setup (login node — ~30 min with Aspera, longer via FTP)
bash 00_setup.sh

# 3. Submit the download + mosdepth array job
sbatch 01_download_and_mosdepth.sh

# 4. Once all 3,202 tasks finish, submit the NGS-PCA job
sbatch 02_run_ngspca.sh

# 5. Compute coverage summary statistics (parallelized across all cores)
bash 03a_mosdepth_coverage_summary.sh

# 6. Collect QC metrics for batch-effect overlay
bash 03_collect_qc.sh
```

---

## Detailed walkthrough

### Step 0: Setup

```bash
bash 00_setup.sh
```

This script performs four tasks:

1. **Creates the directory tree** under `$WORK_DIR`:

   ```
   $WORK_DIR/
   ├── ngs-pca.sif           # Container image
   ├── reference/             # GRCh38 FASTA + index
   ├── manifest.tsv           # Sample manifest (auto-generated)
   ├── crams/                 # Temporary CRAM storage (cleaned per sample)
   ├── mosdepth_output/       # Persistent mosdepth results
   ├── ngspca_output/         # Final PCA results
   └── logs/                  # SLURM job logs
   ```

2. **Pulls the container image** from GitHub Container Registry:

   ```bash
   apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest
   ```

3. **Downloads the GRCh38 reference genome** (~3 GB) from the EBI 1000G FTP. mosdepth needs this to decode CRAM files. The download uses Aspera when available, with wget as fallback.

4. **Downloads the NYGC 30x sequence indexes** (2,504 unrelated + 698 related) from the EBI FTP and builds a unified manifest (`manifest.tsv`) listing every sample ID, CRAM FTP URL, CRAI URL, CRAM MD5, and batch-level metadata parsed from the sequence.index files. The indexes are part of the NYGC 30x data collection:

   - 2,504 samples: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index`
   - 698 related: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index`

   > **Note:** The GitHub-hosted file `1000genomes.high_coverage.GRCh38DH.alignment.index` at [igsr/1000Genomes_data_indexes](https://github.com/igsr/1000Genomes_data_indexes) contains only the 2015 PCR-free pilot data (24 samples, one per population) — **not** the full 3,202-sample NYGC 30x dataset.

   CRAMs are hosted on the ENA FTP (`ftp.sra.ebi.ac.uk`). The manifest includes download fields plus batch annotations:

   ```
   SAMPLE_ID  CRAM_FTP_URL  CRAI_FTP_URL  CRAM_MD5  RELEASE_BATCH  CENTER_NAME  STUDY_ID    INSTRUMENT_MODEL       LIBRARY_NAME
   NA12718    ftp://...     ftp://...crai  923ca...  2504           NYGC         PRJEB31736  Illumina NovaSeq 6000  LP600...
   HG00096    ftp://...     ftp://...crai  355346..  2504           NYGC         PRJEB31736  Illumina NovaSeq 6000  LP600...
   ...
   ```

### Step 1: Download + mosdepth (array job)

```bash
bash 01_download_and_mosdepth.sh
```

This script first scans the manifest and checks which samples already have mosdepth output. It then submits a SLURM array job covering only the tasks that still need to run, skipping samples that are already complete. This avoids wasting queue slots on already-finished work.

By default each task processes one sample, but you can batch multiple samples per task with `SAMPLES_PER_TASK`:

```bash
# Batch 4 samples per task + throttle to 20 concurrent tasks
export SAMPLES_PER_TASK=4
bash 01_download_and_mosdepth.sh
```

Each sample in the task is still processed sequentially as:

1. **Download** the CRAM and CRAI from the ENA using **Aspera** (`ascp`) for high-speed transfer, falling back to `wget` if Aspera fails. Aspera typically achieves 10–100× faster transfers than FTP.

   ```
   ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
     -Tr -Q -l 300m -P33001 -L- \
     era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR323/ERR3239480/NA12718.final.cram \
     /download/
   ```

   > **Network requirements for Aspera:** TCP port 33001 (outgoing) and UDP port 33001 (both directions) must be open. See [IGSR download FAQ](https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/) for details.

2. **Verify** the downloaded CRAM's MD5 checksum against the value in the NYGC sequence index.

3. **Run mosdepth** in 1 kb bins:

   ```
   mosdepth -n -t 2 --by 1000 --fasta GRCh38.fa output_prefix input.cram
   ```

4. **Remove** the CRAM and CRAI to free disk space. Only the mosdepth output (`*.regions.bed.gz`, ~15 MB per sample) is retained.

**To process a subset** (e.g. for testing):

```bash
# First 10 samples only
sbatch --array=1-10 01_download_and_mosdepth.sh

# Retry specific failed tasks
sbatch --array=42,99,256 01_download_and_mosdepth.sh
```

**Monitor progress:**

```bash
# Check how many samples are complete
ls $WORK_DIR/mosdepth_output/*.regions.bed.gz | wc -l

# Check SLURM job status
squeue -u $USER -n 1kG_mosdepth

# View logs for a specific task
cat $WORK_DIR/logs/mosdepth_<JOBID>_<TASKID>.out
```

### Step 2: Run NGS-PCA

```bash
sbatch 02_run_ngspca.sh
```

Once all mosdepth results are available, this single job computes ~200 PCs using the randomized SVD algorithm:

| Parameter | Value | Rationale |
|---|---|---|
| `-numPC` | 200 | ~6% of 3,202 samples — captures population structure and batch effects |
| `-iters` | 10 | Power iterations — sufficient for large cohorts (10K+) |
| `-oversample` | 200 | Oversampling parameter — improves approximation quality |
| `-randomSeed` | 42 | Fixed seed for reproducibility |
| `-threads` | 32 | Parallel loading of mosdepth BED files |
| `-sampleEvery` | 0 | Use all genomic bins (no downsampling) |
| `-distribution` | UNIFORM | Random matrix distribution (see [docs/random-matrix-distribution.md](../../docs/random-matrix-distribution.md)) |
| `-bedExclude` | bundled GRCh38 BED | Excludes SV blacklists, low-mappability, segmental duplications |

**Resource estimates for 3,202 samples:**

| Resource | Estimate |
|---|---|
| Memory | ~256 GB |
| CPUs | 32 |
| Walltime | 6–12 hours |
| Disk (output) | ~500 MB |

### Output files

All output is written to `$WORK_DIR/ngspca_output/`:

| File | Contents |
|---|---|
| `svd.pcs.txt` | 3,202 × 200 sample-by-PC matrix |
| `svd.loadings.txt` | Bin-by-loading matrix |
| `svd.singularvalues.txt` | 200 singular values |
| `svd.bins.txt` | Genomic bins retained after filtering |
| `svd.samples.txt` | Sample identifiers (in row order of `svd.pcs.txt`) |

### Step 3: Collect QC metrics for batch-effect overlay

```bash
# First, compute per-sample autosomal coverage summary statistics
# (parallelized across all available CPU cores)
bash 03a_mosdepth_coverage_summary.sh

# Then collect all QC metrics into a single table
bash 03_collect_qc.sh
```

The coverage summary script (`03a_mosdepth_coverage_summary.sh`) reads the mosdepth `*.regions.bed.gz` files generated in step 1, extracts autosomal chromosome bins (chr1–22), and computes per-sample statistics: mean, median, standard deviation, MAD (median absolute deviation), and IQR (interquartile range). It auto-detects the mosdepth output directory and parallelizes across all available CPU cores using `xargs -P`. Progress is reported to stderr as samples complete.

If `bedtools` is available and `BED_EXCLUDE` points to a valid file, the script also computes **HQ (high-quality)** variants of each statistic using only autosomal bins that do **not** overlap the exclusion BED. The HQ autosomal median is then used by `03_collect_qc.sh` to estimate mitochondrial DNA copy number (mtDNA CN = 2 × chrM_mean / HQ_median).

The QC collection script (`03_collect_qc.sh`) then aggregates all QC sources into a single table — `$WORK_DIR/qc_output/sample_qc.tsv` — that can be directly overlaid on PCA plots to demonstrate which batch effects are captured by each PC.

#### Available QC metrics

| Metric | Source | How obtained |
|---|---|---|
| **Mean autosomal coverage** | mosdepth `.mosdepth.summary.txt` | Free — already computed in step 1 |
| **X coverage ratio** (X/autosomal) | mosdepth `.mosdepth.summary.txt` | Free — already computed in step 1 |
| **Y coverage ratio** (Y/autosomal) | mosdepth `.mosdepth.summary.txt` | Free — already computed in step 1 |
| **Inferred sex** (M/F from coverage) | Derived from X/Y ratios | Free — derived automatically |
| **Mitochondrial coverage ratio** (chrM/autosomal) | mosdepth `.mosdepth.summary.txt` | Free — already computed in step 1 |
| **Median genome coverage** | mosdepth `.mosdepth.global.dist.txt` | Free — already computed in step 1 |
| **% genome ≥ 10× depth** | mosdepth `.mosdepth.global.dist.txt` | Free — already computed in step 1 |
| **% genome ≥ 20× depth** | mosdepth `.mosdepth.global.dist.txt` | Free — already computed in step 1 |
| **Coverage SD** (autosomal bin SD) | `03a_mosdepth_coverage_summary.sh` | Computed from mosdepth regions |
| **Coverage MAD** (autosomal bin MAD) | `03a_mosdepth_coverage_summary.sh` | Computed from mosdepth regions |
| **Coverage IQR** (autosomal bin IQR) | `03a_mosdepth_coverage_summary.sh` | Computed from mosdepth regions |
| **Median bin coverage** (autosomal) | `03a_mosdepth_coverage_summary.sh` | Computed from mosdepth regions |
| **HQ median coverage** (non-excluded bins) | `03a_mosdepth_coverage_summary.sh` | Requires `bedtools` + `BED_EXCLUDE` |
| **HQ SD / MAD / IQR** | `03a_mosdepth_coverage_summary.sh` | Requires `bedtools` + `BED_EXCLUDE` |
| **mtDNA CN** (mitochondrial copy number) | Derived: 2 × chrM_mean / HQ_median | Computed in `03_collect_qc.sh` |
| **Population** (e.g. GBR, YRI) | IGSR sample panel | Downloaded once during setup |
| **Superpopulation** (AFR/AMR/EAS/EUR/SAS) | IGSR sample panel | Downloaded once during setup |
| **Reported sex** | IGSR sample panel | Downloaded once during setup |
| **Relatedness** (unrelated/related) | IGSR sample panel (PED) | Derived from paternal/maternal IDs |
| **Release batch** (2504 or 698) | Manifest (`manifest.tsv`) | Tagged by source sequence.index file |
| **Sequencing center** | Manifest (sequence.index col 6) | Parsed during setup — expected `NYGC` for all |
| **Study ID** | Manifest (sequence.index col 4) | Parsed during setup — study accession |
| **Instrument model** | Manifest (sequence.index col 14) | Parsed during setup (e.g. `Illumina NovaSeq 6000`) |
| **Library name** | Manifest (sequence.index col 15) | Parsed during setup — plate-level batch prefixes |

#### Output table columns

The QC table contains 27 columns organized into 6 groups:

```
SAMPLE_ID  MEAN_AUTOSOMAL_COV  X_COV_RATIO  Y_COV_RATIO  INFERRED_SEX  MITO_COV_RATIO
MEDIAN_GENOME_COV  PCT_GENOME_COV_10X  PCT_GENOME_COV_20X
SD_COV  MAD_COV  IQR_COV  MEDIAN_BIN_COV
HQ_MEDIAN_COV  HQ_SD_COV  HQ_MAD_COV  HQ_IQR_COV  MTDNA_CN
POPULATION  SUPERPOPULATION  REPORTED_SEX  RELATEDNESS
RELEASE_BATCH  CENTER_NAME  STUDY_ID  INSTRUMENT_MODEL  LIBRARY_NAME
```

##### Column descriptions

**Sample identifier**
- `SAMPLE_ID` — Unique sample identifier (e.g. `NA12718`, `HG00096`)

**Coverage metrics (mosdepth summary)**
- `MEAN_AUTOSOMAL_COV` — Mean coverage across autosomal chromosomes (chr1–22). Weighted average of per-chromosome mean coverages. **Use case:** identify sequencing depth batch effects.
- `X_COV_RATIO` — Ratio of chrX coverage to autosomal coverage. Males have ~0.5 (one X copy), females have ~1.0 (two X copies). **Use case:** sex inference, sample swap detection.
- `Y_COV_RATIO` — Ratio of chrY coverage to autosomal coverage. Males have detectable Y coverage (>0.1), females have minimal Y coverage (~0). **Use case:** sex inference.
- `INFERRED_SEX` — Genetic sex inferred from coverage ratios. `M` if Y_COV_RATIO > 0.1 and X_COV_RATIO < 0.75, otherwise `F`. **Use case:** validate against reported sex, detect sample swaps.
- `MITO_COV_RATIO` — Ratio of chrM (mitochondrial) coverage to autosomal coverage. Typically 10–100× higher than nuclear genome. **Use case:** QC flag for mitochondrial enrichment or depletion.

**Coverage distribution metrics (mosdepth global distribution)**
- `MEDIAN_GENOME_COV` — Median depth across the entire genome (all chromosomes). More robust to outliers than mean coverage. **Use case:** assess overall sequencing depth quality.
- `PCT_GENOME_COV_10X` — Percentage of the genome with ≥10× coverage. Standard threshold for variant calling. **Use case:** assess breadth of coverage for variant calling pipelines.
- `PCT_GENOME_COV_20X` — Percentage of the genome with ≥20× coverage. Higher-confidence variant calling threshold. **Use case:** assess high-quality coverage breadth.

**Coverage variability metrics (03a_mosdepth_coverage_summary.sh)**
- `SD_COV` — Standard deviation of per-bin coverage across autosomal chromosomes. High SD indicates uneven coverage. **Use case:** identify library prep or sequencing quality batch effects.
- `MAD_COV` — Median absolute deviation of per-bin coverage. Robust alternative to SD, less sensitive to outliers. **Use case:** robust measure of coverage uniformity.
- `IQR_COV` — Interquartile range (Q3 - Q1) of per-bin coverage. Captures spread of the middle 50% of bins. **Use case:** assess coverage variability, less sensitive to extreme outliers.
- `MEDIAN_BIN_COV` — Median coverage across autosomal bins (1 kb bins from mosdepth). Similar to `MEDIAN_GENOME_COV` but computed from binned data. **Use case:** compare to mean coverage to assess skew.

**High-quality coverage metrics (03a, requires bedtools + BED_EXCLUDE)**
- `HQ_MEDIAN_COV` — Median coverage across autosomal bins that do **not** overlap the exclusion BED. Provides a cleaner baseline by excluding blacklisted, low-mappability, and structurally variant regions. **Use case:** baseline for mtDNA CN estimation and robust coverage assessment.
- `HQ_SD_COV` — Standard deviation of per-bin coverage for non-excluded autosomal bins. **Use case:** assess coverage uniformity in high-quality regions only.
- `HQ_MAD_COV` — MAD of per-bin coverage for non-excluded autosomal bins. **Use case:** robust variability measure excluding problematic regions.
- `HQ_IQR_COV` — IQR of per-bin coverage for non-excluded autosomal bins. **Use case:** spread of coverage in high-quality regions.
- `MTDNA_CN` — Estimated mitochondrial DNA copy number, computed as 2 × chrM_mean_coverage / HQ_MEDIAN_COV. Diploid correction (×2) accounts for the autosomal reference being diploid. **Use case:** detect mitochondrial enrichment/depletion, sample QC.

**Sample metadata (IGSR panel)**
- `POPULATION` — 3-letter population code (e.g. `GBR` = British, `YRI` = Yoruba, `CHB` = Han Chinese). 26 populations total. **Use case:** color PCA plots by ancestry, validate population stratification on PC1/PC2.
- `SUPERPOPULATION` — Continental ancestry group: `AFR` (African), `AMR` (Admixed American), `EAS` (East Asian), `EUR` (European), `SAS` (South Asian). **Use case:** demonstrate population structure in PCA.
- `REPORTED_SEX` — Self-reported sex from IGSR metadata (`M` or `F`). **Use case:** compare to inferred sex for QC.
- `RELATEDNESS` — `unrelated` if both parents unknown in the pedigree (PED), otherwise `related`. **Use case:** identify related individuals (parent-offspring, siblings) that may cluster in PCA.

**Sequencing batch metadata (manifest)**
- `RELEASE_BATCH` — Data release batch: `2504` (unrelated samples) or `698` (related samples). **Use case:** identify batch effects between the two sequencing releases.
- `CENTER_NAME` — Sequencing center (expected `NYGC` for all 1000G 30x samples). **Use case:** multi-center batch effects (not applicable for 1000G 30x).
- `STUDY_ID` — ENA study accession (e.g. `PRJEB31736`, `PRJEB36890`). **Use case:** trace data provenance.
- `INSTRUMENT_MODEL` — Sequencing platform (e.g. `Illumina NovaSeq 6000`). **Use case:** identify instrument-specific batch effects.
- `LIBRARY_NAME` — Library preparation identifier, typically includes plate/batch prefix (e.g. `LP6005442-DNA_A01`). **Use case:** identify library prep batch effects.

#### Joining QC with PCA results

```bash
# Join the QC table with PCA scores on SAMPLE_ID
join -t $'\t' -1 1 -2 1 \
  <(sort $WORK_DIR/qc_output/sample_qc.tsv) \
  <(sort $WORK_DIR/ngspca_output/svd.pcs.txt) \
  > $WORK_DIR/qc_output/pcs_with_qc.tsv
```

The merged table can then be used in R, Python, or any plotting tool to color PCA scatter plots by any QC metric. For example, in R:

```r
library(ggplot2)
d <- read.table("pcs_with_qc.tsv", header = TRUE, sep = "\t")

# Color by superpopulation (population stratification on PC1/PC2)
ggplot(d, aes(PC1, PC2, color = SUPERPOPULATION)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by superpopulation")

# Color by mean coverage (sequencing depth batch effect)
ggplot(d, aes(PC1, PC2, color = MEAN_AUTOSOMAL_COV)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_viridis_c() +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by mean coverage")

# Color by coverage variability (SD — library/sequencing quality batch effect)
ggplot(d, aes(PC1, PC2, color = SD_COV)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_viridis_c() +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by coverage SD")

# Color by coverage IQR (another measure of coverage uniformity)
ggplot(d, aes(PC1, PC2, color = IQR_COV)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_viridis_c() +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by coverage IQR")

# Color by release batch (2504 unrelated vs 698 related)
ggplot(d, aes(PC1, PC2, color = RELEASE_BATCH)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by release batch")

# Color by relatedness (unrelated vs related)
ggplot(d, aes(PC1, PC2, color = RELATEDNESS)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() + labs(title = "1000G 30x — PC1 vs PC2 by relatedness")
```

---

## Configuration reference

All parameters are set in [`config.sh`](config.sh). Override any variable by exporting it before running a script:

```bash
export WORK_DIR=/project/mylab/1000G_highcov
export NUM_PC=300
export NGSPCA_THREADS=64
bash 00_setup.sh
```

| Variable | Default | Description |
|---|---|---|
| `WORK_DIR` | `/scratch/$USER/1000G_highcov` | Root working directory |
| `SIF_IMAGE` | `$WORK_DIR/ngs-pca.sif` | Path to the Apptainer image |
| `REF_FASTA` | `$WORK_DIR/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa` | GRCh38 reference genome |
| `MOSDEPTH_BIN_SIZE` | `1000` | Bin size in bp for mosdepth |
| `MOSDEPTH_THREADS` | `2` | Threads per mosdepth task |
| `NUM_PC` | `200` | Number of PCs to compute |
| `ITERS` | `10` | Power iterations for randomized SVD |
| `OVERSAMPLE` | `200` | Oversampling parameter |
| `RANDOM_SEED` | `42` | Random seed for reproducibility |
| `NGSPCA_THREADS` | `32` | Threads for loading BED files |
| `ASPERA_BANDWIDTH` | `300m` | Aspera transfer speed limit |

---

## Download methods

The pipeline downloads CRAMs from the ENA FTP and reference data from the EBI 1000G FTP. Three download methods are available (in order of speed):

### 1. Aspera (optional, fastest)

[IBM Aspera Connect](https://www.ibm.com/products/aspera/downloads) uses the FASP protocol for high-speed transfers typically 10–100× faster than FTP/HTTP. Install it on your HPC system or load it as a module:

```bash
# Common HPC module names (varies by site):
module load aspera-connect
module load ibm-aspera-connect

# Or install to your home directory from IBM:
# https://www.ibm.com/products/aspera/downloads
```

When `ascp` is on `PATH` and `ASPERA_SSH_KEY` points to a valid key, the pipeline uses it automatically for all downloads. If `ascp` is not found, the pipeline falls back to `wget` without any manual intervention.

```bash
# NYGC CRAMs are on the ENA — use the ENA Aspera user:
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -Tr -Q -l 300m -P33001 -L- \
  era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR323/ERR3239480/NA12718.final.cram \
  ./

# Reference genome is on the 1000G FTP — use the 1000G Aspera user:
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -Tr -Q -l 300m -P33001 -L- \
  fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  ./
```

**Firewall requirements:** TCP and UDP port 33001 must be open. See [IGSR FAQ](https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/) for IP ranges.

### 2. Globus (alternative for restricted networks)

If Aspera ports are blocked, [Globus](https://www.globus.org/) is a reliable alternative. See [IGSR Globus instructions](https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/). This requires manual setup outside the pipeline.

### 3. wget/FTP (fallback)

The pipeline automatically falls back to `wget` if Aspera fails. This is slower but works on any network:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram
```

---

## Adapting for other schedulers

The `#SBATCH` directives in `01_download_and_mosdepth.sh` and `02_run_ngspca.sh` target SLURM. To adapt for PBS/Torque, SGE, or LSF, replace the `#SBATCH` header. The rest of each script runs unchanged.

<details>
<summary>PBS/Torque example</summary>

```bash
#PBS -N 1kG_mosdepth
#PBS -o logs/mosdepth.out
#PBS -e logs/mosdepth.err
#PBS -t 1-3202
#PBS -l nodes=1:ppn=2,mem=4gb,walltime=04:00:00

# Replace SLURM_ARRAY_TASK_ID with PBS_ARRAYID
SLURM_ARRAY_TASK_ID=${PBS_ARRAYID}
```

</details>

---

## Troubleshooting

| Problem | Solution |
|---|---|
| Aspera download fails | Check that TCP/UDP port 33001 is open. The script falls back to wget automatically. |
| `mosdepth: error: could not load index` | Ensure the CRAI file was downloaded alongside the CRAM. |
| `OutOfMemoryError` in step 02 | Increase `--mem` in the SLURM directive and/or set `-Xmx` via `JAVA_TOOL_OPTIONS`. |
| Fewer than 3,202 mosdepth files | Re-run `sbatch --array=<missing_ids> 01_download_and_mosdepth.sh` for failed tasks. |
| `Array task ... exceeds task count` | Use the Step 1 batching submission snippet above to recompute `TOTAL_TASKS` and resubmit with `%` throttling. |
| Manifest is empty or has too few samples | Re-download the indexes: `rm $WORK_DIR/manifest.tsv && bash 00_setup.sh` |
| First ~25–30 tasks run, then many download failures | This is often remote/network connection saturation. Keep `%` throttling on array submissions (for example `%10` to `%30`) and/or increase `SAMPLES_PER_TASK` to reduce concurrent Aspera/wget sessions. |
| Container image pull fails | Check internet access and try: `apptainer pull --force ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest` |

---

## Data sources

- **1000 Genomes 30x data portal:** https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
- **IGSR download FAQ:** https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/
- **NYGC 30x sequence indexes (EBI FTP):** ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
- **ENA projects:** [PRJEB31736](https://www.ebi.ac.uk/ena/browser/view/PRJEB31736) (2,504 unrelated), [PRJEB36890](https://www.ebi.ac.uk/ena/browser/view/PRJEB36890) (698 related), [PRJEB55077](https://www.ebi.ac.uk/ena/browser/view/PRJEB55077) (3,202 combined)
- **GRCh38 reference genome:** ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
- **Byrska-Bishop et al. (2022):** https://doi.org/10.1016/j.cell.2022.08.004
