# NGS-PCA

Principal component analysis of next-generation sequencing coverage data via randomized singular value decomposition.

## Overview

NGS-PCA computes PCs from sequencing coverage across 1 kb genomic bins. The pipeline operates in three stages:

1. **Region selection** — Retain autosomal bins that do not overlap a user-provided exclusion BED file (e.g., structural variant blacklists, low-mappability regions, segmental duplications).
2. **Normalization** — Within each sample, compute log₂ fold change relative to the sample's median bin coverage. Then center each bin to a median of zero across all samples.
3. **Randomized SVD** — Approximate the truncated SVD using the randomized algorithm of [Halko, Martinsson, and Tropp (2011)](https://doi.org/10.1137/090771806), with the power iteration scheme of [Rokhlin, Szlam, and Tygert (2009)](https://doi.org/10.1137/100804139). The implementation is analogous to the [rSVD](https://github.com/erichson/rSVD) R package.

### Random matrix distribution

The Halko et al. algorithm specifies a standard Gaussian random test matrix. This implementation defaults to a **uniform** distribution (`--distribution UNIFORM`), which has been empirically validated to produce equivalent PCs when combined with power (subspace) iterations. A `--distribution GAUSSIAN` option is available for strict adherence to the published algorithm. In practice, the choice has minimal effect on the resulting PCs because the power iteration scheme rapidly converges the column space regardless of the initial random matrix distribution.

## Prerequisites

Install [mosdepth](https://github.com/brentp/mosdepth) (see [installation instructions](https://github.com/brentp/mosdepth#installation)).

## Step 1: Run mosdepth

Compute coverage in 1 kb bins (`--by 1000`):

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

The JAR can be downloaded from a [release](https://github.com/PankratzLab/NGS-PCA/releases) or run via Docker / Apptainer.

```bash
java -Xmx60G -jar ngspca.jar \
  -input /path/to/mosdepthOutput/ \
  -outputDir /path/to/ngsPCA/ \
  -numPC 100 \
  -sampleEvery 0 \
  -threads 24 \
  -iters 40 \
  -randomSeed 42 \
  -oversample 100 \
  -bedExclude ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

### Key parameters

| Parameter | Description | Default |
|---|---|---|
| `-numPC` | Number of PCs to compute. A value around 5 % of the sample size is a reasonable starting point. | 20 |
| `-iters` | Power (subspace) iterations. More iterations improve accuracy at the cost of compute time. 10 is typically sufficient for large cohorts (10 K+); for smaller sample sizes, evaluate a range (e.g., 10–100). | 10 |
| `-oversample` | Oversampling parameter; improves the approximation. At least 10 is recommended. | 200 |
| `-distribution` | Distribution for the initial random matrix: `UNIFORM` or `GAUSSIAN`. | UNIFORM |
| `-randomSeed` | Random seed for reproducibility. | 42 |
| `-threads` | Threads for loading mosdepth files. | 4 |
| `-bedExclude` | BED file of regions to exclude before PCA. | — |
| `-sampleEvery` | Keep every *n*-th bin (0 or 1 = use all bins). | 1 |

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

## Running with Apptainer (HPC)

```bash
apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest

apptainer run \
  --bind /path/to/data:/data \
  ngs-pca.sif \
  -input /data/input.files.txt \
  -outputDir /data/output/ \
  -numPC 100 \
  -sampleEvery 0 \
  -threads 24 \
  -iters 40 \
  -randomSeed 42 \
  -oversample 100 \
  -bedExclude /data/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

## Integration Testing

CI validates deterministic output on every pull request and push to `main`/`master`. See [docs/integration-testing.md](docs/integration-testing.md) for how to reproduce results locally or on HPC with Apptainer.

