# NGS-PCA

Principal component analysis of next-generation sequencing coverage data via randomized singular value decomposition.

## Overview

NGS-PCA computes PCs from sequencing coverage across fixed-width genomic bins. A bin size of 1 kb has been used historically and is recommended as a starting point, but other sizes (e.g., 500 bp or 5 kb) are equally supported. The pipeline operates in three stages:

1. **Region selection** — Retain autosomal bins that do not overlap a user-provided exclusion BED file (e.g., structural variant blacklists, low-mappability regions, segmental duplications).
2. **Normalization** — Within each sample, compute log₂ fold change relative to the sample's median bin coverage. Then center each bin to a median of zero across all samples.
3. **Randomized SVD** — Approximate the truncated SVD using the randomized algorithm of [Halko, Martinsson, and Tropp (2011)](https://doi.org/10.1137/090771806), with the power iteration scheme of [Rokhlin, Szlam, and Tygert (2009)](https://doi.org/10.1137/080736417). The implementation is analogous to the [rSVD](https://github.com/erichson/rSVD) R package.

## Prerequisites

Install [mosdepth](https://github.com/brentp/mosdepth) (see [installation instructions](https://github.com/brentp/mosdepth#installation)).

## Step 1: Run mosdepth

Compute coverage in fixed-width bins. 1 kb (`--by 1000`) is the recommended starting point but other bin sizes are supported:

```bash
mosdepth -n -t 1 --by 1000 \
  --fasta /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  output_prefix input.bam
```

To process a directory of BAM/CRAM files in parallel:

```bash
ref=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa
dirOfBams=/path/to/bams/
mosdepthResultsDir=/path/to/mosdepthOutput/

find "$dirOfBams" -type f -name "*.bam" \
  | parallel -j 24 "mosdepth -n -t 1 --by 1000 --fasta $ref $mosdepthResultsDir{/.}.by1000 {}"
```

## Step 2: Run NGS-PCA

The parameters below reflect those used in a real-world run of **142,000 samples** on a `ram2t` node (1800 GB RAM, 120 threads, ~60 hours walltime). Adjust memory (`-Xmx`), `-threads`, `-numPC`, and other parameters to match your cohort size and available resources. The defaults are more conservative (see table below).

### On HPC with Apptainer (recommended)

Pull the pre-built image and run directly — no Java or Maven installation required. Exclusion BED files are bundled inside the image at `/app/resources/`:

```bash
apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest

numPC=500
bedExclude=/app/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
iters=10
oversample=0

apptainer run \
  --bind /path/to/data:/data \
  ngs-pca.sif \
  -input /data/input.files.txt \
  -outputDir /data/ngsPCA/ \
  -numPC $numPC \
  -sampleEvery 0 \
  -threads 120 \
  -bedExclude $bedExclude \
  -iters $iters \
  -oversample $oversample
```

### Download pre-built JAR

Pre-built fat JARs are published automatically by GitHub Actions on every tagged release. Download the latest `ngspca-vX.Y.jar` from the [Releases page](https://github.com/jlanej/NGS-PCA/releases). Exclusion BED files are provided in the `resources/` directory of this repository — clone or download the repo alongside the JAR:

```bash
git clone https://github.com/jlanej/NGS-PCA.git

numPC=500
bedExclude=NGS-PCA/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
iters=10
oversample=0

java -Xmx1800G -jar ngspca-vX.Y.jar \
  -input /path/to/mosdepthOutput/ \
  -outputDir /path/to/ngsPCA/ \
  -numPC $numPC \
  -sampleEvery 0 \
  -threads 120 \
  -bedExclude $bedExclude \
  -iters $iters \
  -oversample $oversample
```

### Building from source

Requires Java 11+ and Maven 3.6+:

```bash
git clone https://github.com/jlanej/NGS-PCA.git
cd NGS-PCA
mvn -B package --file ngspca/pom.xml

numPC=500
bedExclude=resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
iters=10
oversample=0

java -Xmx1800G -jar ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  -input /path/to/mosdepthOutput/ \
  -outputDir /path/to/ngsPCA/ \
  -numPC $numPC \
  -sampleEvery 0 \
  -threads 120 \
  -bedExclude $bedExclude \
  -iters $iters \
  -oversample $oversample
```

### Key parameters

| Parameter | Description | Default |
|---|---|---|
| `-input` | Directory of mosdepth result files (`*.regions.bed.gz`), a file listing one path per line, or (with `-matrix`) a pre-computed matrix. | — |
| `-outputDir` | Directory where PCA results and auxiliary files will be written. | — |
| `-numPC` | Number of PCs to compute. A value around 5 % of the sample size is a reasonable starting point. | 20 |
| `-iters` | Power (subspace) iterations. More iterations improve accuracy at the cost of compute time. 10 is typically sufficient for large cohorts (10 K+); for smaller sample sizes, evaluate a range (e.g., 10–100). | 10 |
| `-oversample` | Oversampling parameter; improves the approximation. At least 10 is recommended. | 200 |
| `-distribution` | Distribution for the initial random matrix: `UNIFORM` or `GAUSSIAN`. | UNIFORM |
| `-randomSeed` | Random seed for reproducibility. | 42 |
| `-threads` | Threads for loading mosdepth files. | 4 |
| `-bedExclude` | BED file of regions to exclude before PCA. | — |
| `-sampleEvery` | Keep every *n*-th bin (0 or 1 = use all bins). | 1 |
| `-overwrite` | Flag: overwrite existing temporary (cached) files and recompute each step. | false |

### Exclusion BED files

Pre-built exclusion BED files are provided in `resources/`:

- **WGS (GRCh38):** [`ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz`](https://github.com/PankratzLab/NGS-PCA/blob/master/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz)
- **WES (UKB, GRCh38):** [`ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.xgen.sorted.merge.contig.bed.gz`](https://github.com/PankratzLab/NGS-PCA/blob/master/resources/GRCh38/UKB_WES/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.xgen.sorted.merge.contig.bed.gz) — combines the WGS exclusion regions with 20 kb-buffered xgen exome targets ([source](http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=3801)).

For custom WES analyses, concatenate the WGS exclusion BED with a 20 kb-buffered version of the capture target BED.

### Output files

| File | Contents |
|---|---|
| `svd.pcs.txt` | Sample × PC matrix |
| `svd.loadings.txt` | Bin × loading matrix |
| `svd.singularvalues.txt` | Singular values per PC |
| `svd.bins.txt` | Genomic bins retained after filtering |
| `svd.samples.txt` | Sample identifiers |

## Random matrix distribution

This implementation defaults to a **uniform** distribution (`-distribution UNIFORM`); `-distribution GAUSSIAN` is available for strict adherence to the Halko et al. paper. In practice the choice has negligible effect on results when power iterations are used. See [docs/random-matrix-distribution.md](docs/random-matrix-distribution.md) for full details and validation results.

## Integration Testing

CI validates deterministic output on every pull request and push to `main`/`master`. See [docs/integration-testing.md](docs/integration-testing.md) for how to reproduce results locally or on HPC with Apptainer.

## Future Work

See [docs/future-work.md](docs/future-work.md) for a summary of known performance, scientific, and code-quality improvements that are planned or under consideration.

