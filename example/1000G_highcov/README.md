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
| — | `install_aspera.sh` | Builds `aspera.def` into `$WORK_DIR` — a current `ascp` plus the EMBL-EBI public key — and verifies by logging in to ENA. **Run automatically by stage 1**; only invoke it directly to re-provision | Interactive / login node |
| **1** | `01_download_and_mosdepth.sh` | For each sample: download CRAM (Aspera → aria2c → parallel curl → wget) → run mosdepth → remove CRAM | Array job (3,202 tasks) |
| **2** | `02_run_ngspca.sh` | Run NGS-PCA on all mosdepth results → ~200 PCs | Single large-memory job |
| **3a** | `03a_mosdepth_coverage_summary.sh` | Compute autosomal coverage stats (mean, median, SD, MAD, IQR) and HQ statistics (non-excluded bins) from mosdepth output | Parallelized (all cores) |
| **3** | `03_collect_qc.sh` | Aggregate per-sample QC into one table for batch-effect overlay | Interactive or short job |

All three stages use the **same container image** (`ghcr.io/jlanej/ngs-pca:latest`), which bundles:

- **NGS-PCA** — randomized SVD for coverage PCA
- **mosdepth** v0.3.9 — fast BAM/CRAM depth calculation
- Pre-built **exclusion BED files** for GRCh38

No additional software installation is required beyond [Apptainer](https://apptainer.org/) (Singularity). IBM Aspera Connect (`ascp`) is an optional system-level tool for faster downloads; if it is missing, too old, or disabled via `USE_ASPERA=0`, the pipeline falls back to parallel HTTPS downloads (`aria2c`, else `curl` byte ranges) and finally to single-stream `wget`. See [step 1](#step-1-download-manager--spawned-mosdepth-jobs) for the Aspera client version requirement.

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

# 3. Start the download manager; it spawns one mosdepth job per verified sample
bash 01_download_and_mosdepth.sh

# 4. Once all 3,202 tasks finish, submit the NGS-PCA job
sbatch 02_run_ngspca.sh

# 5. Compute coverage summary statistics (parallelized across all cores)
bash 03a_mosdepth_coverage_summary.sh

# 6. Collect QC metrics for batch-effect overlay
bash 03_collect_qc.sh
```

---

## Pre-computed example results

The `output/` directory in this repository contains the full results of running the pipeline on all 3,202 samples, committed here so you can explore the outputs without re-running the pipeline.

> **Note:** The bin-loadings file (`svd.loadings.txt`) is not committed because it is very large (~1 GB). All other output files are included.

### `output/ngspca_output/`

| File | Rows | Description |
|---|---|---|
| `svd.pcs.txt` | 3,203 (1 header + 3,202 samples) | Sample-by-PC matrix: 200 principal components for each of the 3,202 samples |
| `svd.singularvalues.txt` | 201 (1 header + 200 PCs) | Singular values for each of the 200 computed PCs (proxy for variance explained) |
| `svd.bins.txt` | 142,070 | Genomic bins (mosdepth regions) retained after filtering, in the row order used by the loadings matrix |
| `svd.samples.txt` | 3,202 | Sample identifiers in the row order of `svd.pcs.txt` |

> **ID format note:** Sample IDs in `svd.pcs.txt` / `svd.samples.txt` include a `.by1000.` suffix (e.g., `HG00096.by1000.`), while `sample_qc.tsv` below uses bare IDs (e.g., `HG00096`). To join QC ↔ PCA outputs, strip the suffix from the PCA sample IDs (for example, remove a trailing `.by1000.` substring before merging).
### `output/qc_output/`

| File | Rows | Description |
|---|---|---|
| `sample_qc.tsv` | 3,203 (1 header + 3,202 samples) | Per-sample QC table with 28 columns: mosdepth summary stats, coverage dispersion, HQ autosomal coverage, mtDNA CN, population/sex panel metadata, and sequencing batch annotations. See [QC table columns](#output-table-columns) below. |
| `mosdepth_coverage_summary.tsv` | 3,203 (1 header + 3,202 samples) | Per-sample autosomal coverage statistics computed by `03a_mosdepth_coverage_summary.sh`: mean, median, SD, MAD, and IQR for all autosomal bins and for high-quality (non-excluded) bins only |

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

### Step 1: Download manager + spawned mosdepth jobs

```bash
bash 01_download_and_mosdepth.sh
```

This submits **one long-running download manager job**. The manager holds `DOWNLOAD_SLOTS`
(default 6) concurrent transfers — the WAN is the scarce, failure-prone resource, and few fast
streams beat many starved ones — and the moment a sample's CRAM verifies, its mosdepth run is
submitted as an independent job (`01b_mosdepth_sample.sh`, name `1kG_mosdepth_<sample>`), so
mosdepth concurrency is whatever the scheduler grants and never waits on a download slot.

The sweep is idempotent three ways: samples with mosdepth output are skipped, samples whose
mosdepth job is already queued or running are skipped (one `squeue` snapshot at start), and a
verified CRAM already on disk is reused without re-downloading. Re-running the script continues
where the last sweep stopped. The manager pauses when `MAX_LOCAL_CRAMS` CRAMs sit on disk
awaiting mosdepth, so a stalled queue cannot fill scratch, and it exits non-zero if any sample
failed, listing them.

```bash
# Smoke-test with the first 10 samples that need work
DOWNLOAD_LIMIT=10 bash 01_download_and_mosdepth.sh

# Run the manager in place instead of as a job (login/DTN node; tmux advised)
DOWNLOADER_LOCAL=1 bash 01_download_and_mosdepth.sh
```

Each sample is processed as:

1. **Download** the CRAM and CRAI from the ENA, trying each transport in descending order of speed and stopping at the first that succeeds:

   | Order | Method | Notes |
   |---|---|---|
   | 1 | **Aspera** (`ascp`) | FASP; fastest when the client is new enough (see below). Skip with `USE_ASPERA=0`. |
   | 2 | **aria2c** | `DOWNLOAD_CONNECTIONS` (default 4) parallel HTTPS streams, resumable. |
   | 3 | **parallel `curl`** | Same idea using only `curl` byte-range requests, assembled and size-checked. |
   | 4 | **`wget`** | Single-stream last resort. |

   ```
   ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
     -Tr -Q -l 50m -P33001 -L- -k 2 \
     era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR323/ERR3239480/NA12718.final.cram \
     /download/
   ```

   The rate is per task, and `ASPERA_BANDWIDTH` × concurrent tasks is the aggregate the site asks
   of EBI — size it to a few Gbit/s (the defaults, 60 × 50m, are ~3 Gbit/s). Asking for far more
   is not faster: 200 × 300m collapsed into packet loss, sessions crawling at ~13 Mbit/s, and
   `Session data transfer timeout` errors — while also getting the HTTPS fallbacks throttled.
   ascp retries `ASPERA_RETRIES` times before falling back, resuming its partial with `-k 2`.

   > **Network requirements for Aspera:** TCP port 33001 (outgoing) and UDP port 33001 (both directions) must be open. See [IGSR download FAQ](https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/) for details.

   > **Aspera client version:** `fasp.sra.ebi.ac.uk` now runs OpenSSH 8.7 and offers only modern SSH transport algorithms — `curve25519-sha256`, `ecdh-sha2-nistp*` and `diffie-hellman-group-exchange-sha256` for key exchange, and *encrypt-then-MAC* MACs only (`umac-128-etm@`, `hmac-sha2-{256,512}-etm@`). The libssh2 bundled with **ascp 3.9.x supports none of these**, so its handshake dies during algorithm negotiation:
   >
   > ```
   > [libssh2] Failure Event: -5 - Unable to exchange encryption keys
   > ERR [asssh] SSH connection startup failed, err = -5
   > ascp: failed to authenticate, exiting.
   > ```
   >
   > The final line is misleading: the session never reaches authentication, so **changing the `-i` key file does not help**. This is a *server-side* change — a client that worked for years breaks without being touched.
   >
   > **Two things changed, and you need both fixes.** Upgrading the client alone is not enough: the anonymous key `asperaweb_id_dsa.openssh` is **no longer accepted by ENA** either. Verified directly — same client and server, old key rejected, new key accepted. EMBL-EBI replaced it with an RSA key for the public accounts (`fasp-public`, `fasp-ml`, `era-fasp`), documented in [KB0011597](https://embl.service-now.com/kb?id=kb_article_view&sysparm_article=KB0011597) and [KB0011565](https://embl.service-now.com/kb?id=kb_article_view&sysparm_article=KB0011565).
   >
   > **This is handled for you.** Running `bash 01_download_and_mosdepth.sh` provisions Aspera on the submit host before it submits the manager job, so there is no extra step to remember.
   >
   > Under the hood that runs `install_aspera.sh`, which builds `aspera.def` into `$WORK_DIR/aspera/aspera.sif`, pairing a checksum-pinned [Aspera Connect 4.2.13](https://www.ibm.com/products/aspera/downloads) with the current EBI key and verifying with a real `--mode=test-login`. `config.sh` then picks up `ASPERA_BIN` and `ASPERA_SSH_KEY`. Invoke it directly only to re-provision. If provisioning fails for any reason, submission continues and downloads use the parallel-HTTPS paths.
   >
   > **No key is stored in this repository.** EBI publishes it unauthenticated, but only inside a JavaScript-rendered knowledge-base page: plain HTTP cannot reach it (the layout API carries no article text; the documented KB APIs return 401). The image's first build stage therefore renders the page with headless Chromium and extracts the PEM. Chromium stays in that stage and is not part of the final image. The extraction refuses to guess if the article ever contains more than one distinct key, and the `%test` turns a rotation or page change into a build failure rather than a silent fallback.
   >
   > Building from a definition file needs `apptainer build --fakeroot`, which some sites disable. If that is your case, use `USE_ASPERA=0` and the parallel HTTPS paths.
   >
   > Until Aspera is provisioned, set `USE_ASPERA=0` to skip the handshake attempts and the FASP log noise — the parallel HTTPS paths are used instead.

2. **Verify** the downloaded CRAM's MD5 checksum against the value in the NYGC sequence index.
   A mismatch deletes the pair and re-downloads once — over HTTPS, since a completed Aspera
   transfer that fails MD5 means systematic corruption, which also sets a persistent distrust
   flag (`$DOWNLOAD_STATE_DIR/aspera_disabled`) that keeps the whole run and later sweeps off
   Aspera until the file is removed.

3. **Spawn mosdepth** as its own job the moment verification passes (twice, timed, when
   `COMPARE_FAST_MODE=1`):

   ```
   mosdepth -n -t 2 --by 1000 --fasta GRCh38.fa output_prefix input.cram
   ```

4. The mosdepth job **removes** the CRAM and CRAI when it finishes. Only the mosdepth output
   (`*.regions.bed.gz`, ~15 MB per sample) is retained.

**Monitor progress:**

```bash
# Check how many samples are complete
ls $WORK_DIR/mosdepth_output/*.regions.bed.gz | wc -l

# The manager and the mosdepth jobs it has spawned
squeue -u $USER -n 1kG_download
squeue -u $USER | grep 1kG_mosdepth

# The manager's narration, one line per dispatch
tail -f $WORK_DIR/logs/download_manager_<JOBID>.out

# One sample's download story, or its mosdepth job
cat $WORK_DIR/logs/download_<SAMPLE>.log
cat $WORK_DIR/logs/mosdepth_<SAMPLE>_<JOBID>.out
```

**Triage failures.** The manager's summary lists failed samples and exits non-zero when any
exist; each failure's full story is in its own `download_<SAMPLE>.log` (ascp logs its FASP noise
there too — an ERR line does not mean the sample failed if a fallback then succeeded). Re-running
`bash 01_download_and_mosdepth.sh` retries exactly the samples still missing output.

```bash
# Which transport completed each download
grep -h ": .* download complete" $WORK_DIR/logs/download_*.log | sort | uniq -c

# Samples that failed every transport in the last sweep
grep -rlx "fail" $WORK_DIR/download_state/results | wc -l
```

No aria2c completions at all usually means it is not installed on the node running the manager
(`command -v aria2c`; it is the best HTTPS transport, worth a `module load`). MD5 mismatches are
attributed the same way:

```bash
for f in $(grep -l "MD5 mismatch" $WORK_DIR/logs/download_*.log); do
  grep -m1 -o "CRAM: .* download complete" "$f"
done | sort | uniq -c
```

A scattering across transports is transfer noise the in-sweep re-download absorbs. Mismatches on
essentially **every** Aspera transfer while HTTPS runs clean is a different animal — systematic
payload corruption with ascp exiting success (measured once at 506 of 506). The distrust flag
handles it automatically from the first occurrence; set `USE_ASPERA=0` to skip even the first
wasted transfer while it is investigated on the host.

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

### Optional: mosdepth fast-mode comparison

mosdepth's `--fast-mode` skips CIGAR-aware depth and mate-overlap correction, and its README
recommends it for most use cases. Whether it changes the PCs is an empirical question this
pipeline can answer about itself: both effects are largely per-sample multiplicative shifts,
which the log₂ fold change against each sample's own median absorbs, so near-perfect PC
correlation is the expectation — and if it holds, fast mode is a free speedup on the most
expensive stage.

```bash
# 1. Rerun step 1 with the comparison on: each sample's mosdepth job runs
#    mosdepth twice (standard -> mosdepth_output/, fast -> mosdepth_output_fast/),
#    recording each run's wall time. Downloads are shared and untimed; the two
#    runs alternate order so cache warming cancels out in aggregate.
export COMPARE_FAST_MODE=1
bash 01_download_and_mosdepth.sh

# 2. NGS-PCA on each tree — same parameters, seed, and exclusion bed
sbatch 02_run_ngspca.sh
MOSDEPTH_DIR="${WORK_DIR}/mosdepth_output_fast" \
NGSPCA_OUTPUT="${WORK_DIR}/ngspca_output_fast" sbatch 02_run_ngspca.sh

# 3. Optional: QC tables for each tree, so depth-derived phenotypes -
#    MTDNA_CN, coverage ratios, inferred sex - are compared between modes too
sbatch 03a_mosdepth_coverage_summary.sh && bash 03_collect_qc.sh
MOSDEPTH_DIR="${WORK_DIR}/mosdepth_output_fast" QC_OUTPUT="${WORK_DIR}/qc_output_fast" \
  sbatch 03a_mosdepth_coverage_summary.sh
MOSDEPTH_DIR="${WORK_DIR}/mosdepth_output_fast" QC_OUTPUT="${WORK_DIR}/qc_output_fast" \
  bash 03_collect_qc.sh

# 4. Correlations, runtime boxplots, and a written summary
#    (host python3; the figures need matplotlib, the tables do not)
bash 04_fast_mode_eval.sh
```

The evaluation lands in `$QC_OUTPUT/fast_mode_eval/`: per-PC Pearson r with singular-value
ratios (`pc_correlation.tsv`), wall-time statistics and per-sample paired speedups
(`timing_summary.tsv`), a six-panel figure (`fast_mode_summary.png`), and
`fast_mode_report.md` with the headline numbers, including an order check that would expose
cache-warming bias. Correlations are evaluated as |r|, since a singular vector's sign is
arbitrary. When both QC tables exist, the eval also writes `qc_concordance.tsv` and
`fast_mode_qc.png` — per-column Pearson r plus the median relative difference, so a uniform
bias (the signature of skipped mate-overlap correction, e.g. in `MTDNA_CN`) is visible even
where correlation is perfect.

Enabling the comparison reprocesses any sample missing either tree, so every timed pair comes
from one node and one download; a cohort already processed without it will re-download. Timing
rows record the mosdepth version — keep one version per comparison (the image pins one, see the
Dockerfile).

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
| **Superpopulation** (AFR/AMR/EAS/EUR/SAS) | Derived from population code | Pop→superpop lookup in `03_collect_qc.sh` |
| **Reported sex** | IGSR sample panel | Downloaded once during setup |
| **Family role** (father/mother/child/unrel…) | IGSR sample panel (PED col 8) | PED "Relationship" field |
| **Relatedness** (unrelated/related) | Derived from PED parental IDs | Both parents "0" → unrelated |
| **Release batch** (2504 or 698) | Manifest (`manifest.tsv`) | Tagged by source sequence.index file |
| **Sequencing center** | Manifest (sequence.index col 6) | Parsed during setup — expected `NYGC` for all |
| **Study ID** | Manifest (sequence.index col 4) | Parsed during setup — study accession |
| **Instrument model** | Manifest (sequence.index col 14) | Parsed during setup (e.g. `Illumina NovaSeq 6000`) |
| **Library name** | Manifest (sequence.index col 15) | Parsed during setup — plate-level batch prefixes |

#### Output table columns

The QC table contains 28 columns organized into 6 groups:

```
SAMPLE_ID  MEAN_AUTOSOMAL_COV  X_COV_RATIO  Y_COV_RATIO  INFERRED_SEX  MITO_COV_RATIO
MEDIAN_GENOME_COV  PCT_GENOME_COV_10X  PCT_GENOME_COV_20X
SD_COV  MAD_COV  IQR_COV  MEDIAN_BIN_COV
HQ_MEDIAN_COV  HQ_SD_COV  HQ_MAD_COV  HQ_IQR_COV  MTDNA_CN
POPULATION  SUPERPOPULATION  REPORTED_SEX  FAMILY_ROLE  RELATEDNESS
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
- `MTDNA_CN` — Estimated mitochondrial DNA copy number, computed as 2 × chrM_mean_coverage / HQ_MEDIAN_COV. Diploid correction (×2) accounts for the autosomal reference being diploid. Computed in `03_collect_qc.sh` using MITO_COV_RATIO × MEAN_AUTOSOMAL_COV as chrM_mean. `NA` when HQ_MEDIAN_COV is unavailable (requires `bedtools` + `BED_EXCLUDE`). **Use case:** detect mitochondrial enrichment/depletion, sample QC.

**Sample metadata (IGSR panel)**
- `POPULATION` — 3-letter population code (e.g. `GBR` = British, `YRI` = Yoruba, `CHB` = Han Chinese). 26 populations total. **Use case:** color PCA plots by ancestry, validate population stratification on PC1/PC2.
- `SUPERPOPULATION` — Continental ancestry group derived from `POPULATION` via a lookup table: `AFR` (African), `AMR` (Admixed American), `EAS` (East Asian), `EUR` (European), `SAS` (South Asian). The IGSR PED file does not contain a superpopulation column; it is mapped from the 26 population codes. **Use case:** demonstrate population structure in PCA.
- `REPORTED_SEX` — Self-reported sex from IGSR metadata (`M` or `F`). **Use case:** compare to inferred sex for QC.
- `FAMILY_ROLE` — PED Relationship field (column 8). Describes the individual's declared role in the family pedigree: `unrel` (unrelated, no family members in the dataset), `father`, `mother`, `child`, `pat` (paternal grandfather/uncle), `mat` (maternal grandmother/aunt), or other family descriptors. **Use case:** identify family structure and understand sample composition.
- `RELATEDNESS` — Derived from PED parental IDs: `unrelated` if both paternal and maternal IDs are `0` (founders or unrelated individuals), otherwise `related` (individuals with known parents in the pedigree). **Note:** `FAMILY_ROLE` and `RELATEDNESS` are complementary but not fully concordant — a `father` typically has `RELATEDNESS = unrelated` because founders have `0/0` parental IDs, while his `child` has `RELATEDNESS = related` because parental IDs are filled in. A small number of `unrel` individuals have non-zero parental IDs, giving `RELATEDNESS = related`. **Use case:** identify related individuals (parent-offspring, siblings) that may cluster in PCA.

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
| `BED_EXCLUDE` | `../../resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz` (from `config.sh` dir, fallback `/app/resources/...`) | Exclusion BED for NGS-PCA and HQ autosomal coverage stats |
| `ASPERA_BANDWIDTH` | `100m` | Per-transfer FASP target rate; × `DOWNLOAD_SLOTS` = aggregate asked of EBI |
| `ASPERA_RETRIES` | `3` | ascp attempts per file, resuming the partial with `-k 2`, before HTTPS fallback |
| `DOWNLOAD_SLOTS` | `6` | Concurrent transfers held by the download manager |
| `DOWNLOAD_CONNECTIONS` | `16` | HTTPS range streams per download; × `DOWNLOAD_SLOTS` = connections at ENA |
| `MAX_LOCAL_CRAMS` | `60` | Downloads pause while this many CRAMs await mosdepth (~16 GB each) |
| `DOWNLOAD_LIMIT` | `0` | Stop the sweep after this many dispatches (0 = whole manifest) |
| `MOSDEPTH_MEM` / `MOSDEPTH_TIME` | `8G` / `06:00:00` | Resources for each spawned mosdepth job |

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

The `-l 300m` in the manual examples below suits a *single* transfer; the pipeline's per-task
rate is `ASPERA_BANDWIDTH`, sized so that rate × concurrent tasks stays within a few Gbit/s.

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

### 2. Globus (bulk staging — best when direct downloads crawl)

When the manager's own transports are slow or unreliable from your site, [Globus](https://www.globus.org/)
is usually the fastest path: EBI exposes ENA through the **"EMBL-EBI Public Data"** collection
(UUID `47772002-3e5b-4fd3-b97c-18cee38d6df2`), with `ftp.sra.ebi.ac.uk` mounted under `/vol1` —
so every manifest URL is a Globus path with the host stripped
([ENA docs](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html),
[IGSR FAQ](https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/)).
Globus transfers are checksummed, resumed, and retried by the service itself, and run
endpoint-to-endpoint on data-transfer infrastructure rather than through a compute node.

The pipeline is built to receive this. Staged files named `<SAMPLE>.cram` in `$CRAM_DIR` are
treated as "already present": the manager skips the download, **still verifies the manifest
MD5**, and goes straight to spawning mosdepth. `globus_batch_from_manifest.sh` writes **two**
batch files that do the renaming during the transfer — indexes and CRAMs apart, because Globus
chooses its own file order *within* a task (a mixed task was observed running the CRAMs first,
completing no pair for hours). Submitted as two tasks, the CRAI task lands in minutes, and from
then on every CRAM that finishes completes its pair:

```bash
# 1. Generate the batch files (every sample still needing mosdepth output;
#    files already on disk are handled by --sync-level checksum below)
bash globus_batch_from_manifest.sh

# 2. Find your site's collection UUID, then submit BOTH - CRAIs first, and
#    always with --sync-level checksum: complete files are skipped nearly
#    free, partial or corrupt ones re-transferred, which makes regenerate +
#    resubmit the whole restart story after a cancel or failure.
#    (globus-cli: try 'module spider globus' first, else pip install
#    globus-cli - or run it from a laptop; transfers are endpoint-to-endpoint
#    and the submitting machine never touches the data. Then: globus login)
globus endpoint search "<your site>"
globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
  --batch globus_batch_crai.txt --sync-level checksum --label "1kG crai" --notify failed,inactive
globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
  --batch globus_batch_cram.txt --sync-level checksum --label "1kG cram" --notify failed,inactive

# 3. While the transfer runs, dispatch mosdepth as pairs land. The watcher
#    polls for CRAM+CRAI pairs whose mtimes have settled, submits each one's
#    mosdepth job with the manifest MD5 (verified on the compute node, in
#    parallel), and never downloads or deletes anything - so it is safe
#    alongside the transfer, unlike the manager. A pair dispatched while
#    secretly still in flight fails its MD5 loudly and is re-dispatched once
#    it settles again; content that stops changing and still fails is parked.
bash 01c_dispatch_staged.sh

# 4. Watch the transfer (or the web app's Activity page) and the watcher
globus task list
tail -f $WORK_DIR/logs/stagewatch_<JOBID>.out

# 5. When the transfer completes, one final sweep: verifies and dispatches
#    whatever remains, and downloads anything Globus missed
bash 01_download_and_mosdepth.sh
```

Two operational notes. **Stop the manager while Globus writes into `$CRAM_DIR`**
(`scancel -n 1kG_download`) — two writers racing on the same sample is how partial files
happen; resweep after the transfer instead. And mind disk: the full cohort staged at once peaks
near 50 TB in `$CRAM_DIR`, draining as mosdepth jobs finish (the manager's `MAX_LOCAL_CRAMS`
backpressure deliberately ignores pre-staged files — they consume no new space and dispatching
them frees it). To stage in waves, `split -l 1000 globus_batch.txt` and submit the pieces as
separate transfers.

### 3. wget/FTP (fallback)

The pipeline automatically falls back to `wget` if Aspera fails. This is slower but works on any network:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram
```

---

## Adapting for other schedulers

The `#SBATCH` headers target SLURM, and the download manager submits its per-sample mosdepth
jobs with one `sbatch` call (in `download_sample`, `01_download_and_mosdepth.sh`) plus a
`squeue` snapshot for its in-flight check. On PBS/SGE/LSF: run the manager in place with
`DOWNLOADER_LOCAL=1` (a login or DTN node suffices — it only downloads), replace that `sbatch`
call with your scheduler's equivalent (`qsub`/`bsub` wrapping `01b_mosdepth_sample.sh SAMPLE
LINE`), and either translate the in-flight check or leave it — without `squeue` the manager
still skips every sample whose output exists, it just cannot see jobs that are queued but not
yet finished. `02_run_ngspca.sh` needs only its header translated.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| Aspera download fails | Check that TCP/UDP port 33001 is open. The script falls back to wget automatically. |
| `mosdepth: error: could not load index` | Ensure the CRAI file was downloaded alongside the CRAM. |
| `OutOfMemoryError` in step 02 | Increase `--mem` in the SLURM directive and/or set `-Xmx` via `JAVA_TOOL_OPTIONS`. |
| Fewer than 3,202 mosdepth files | Re-run `bash 01_download_and_mosdepth.sh` — the sweep retries exactly the samples still missing output. |
| Manifest is empty or has too few samples | Re-download the indexes: `rm $WORK_DIR/manifest.tsv && bash 00_setup.sh` |
| Widespread download failures or MD5 mismatches | See the triage in Step 1: transports and mismatches are attributed per sample in `logs/download_<sample>.log`, and a corrupt Aspera payload flips the run to HTTPS on its own. Lower `DOWNLOAD_SLOTS` or `ASPERA_BANDWIDTH` if timeouts dominate. |
| Container image pull fails | Check internet access and try: `apptainer pull --force ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest` |

---

## Data sources

- **1000 Genomes 30x data portal:** https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
- **IGSR download FAQ:** https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/
- **NYGC 30x sequence indexes (EBI FTP):** ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
- **ENA projects:** [PRJEB31736](https://www.ebi.ac.uk/ena/browser/view/PRJEB31736) (2,504 unrelated), [PRJEB36890](https://www.ebi.ac.uk/ena/browser/view/PRJEB36890) (698 related), [PRJEB55077](https://www.ebi.ac.uk/ena/browser/view/PRJEB55077) (3,202 combined)
- **GRCh38 reference genome:** ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
- **Byrska-Bishop et al. (2022):** https://doi.org/10.1016/j.cell.2022.08.004
