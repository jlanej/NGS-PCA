# NGS-PCA
Methods for running PCA on NGS data


## Install [mosdepth](https://github.com/brentp/mosdepth)

Install mosdepth using the instructions from https://github.com/brentp/mosdepth#installation.
There are lots of ways to do this, including downloading a linux executable, as a Docker image, or using brew.

## Run mosdepth on bams or crams
The "--by 1000" (compute coverage on 1000bp bins) is really the only important argument, and each run is going to look something like this:
`mosdepth -n -t 1 --by 1000 --fasta /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa output_filename input_filename.bam`

but here is an example script to iterate over all BAM files in a directory (could be run on CRAM files as well)

```bash
ref=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa

dirOfBams=/path/to/bams/
mosdepthResultsDir=/path/to/mosdepthOutput/
mosdepthThreads=1
parallelThreads=24

find "$dirOfBams" -type f -name "*.bam" \
|parallel -j $parallelThreads "mosdepth -n -t $mosdepthThreads --by 1000 --fasta $ref $mosdepthResultsDir{/.}.by1000 {}"
```

## Run ngs pca

The number of PCs to compute should be in the range of 5% of your sample size and likely far more than you'll actually use.
We're still working on optimizing the number of iterations and oversample parameter - but this should be reasonable. For smaller sample sizes it may be worth testing a range of `-iters` arguments (10,20,30,40,50,100,etc). More iterations increases the accuracy of the PCs, but also increases compute time. For larger sample sizes (10K+), 10 iterations appears to be sufficient.


This will generate svd.pcs.txt in the output directory

```bash
ngsPCAOutputDir=/path/to/ngsPCA/
ngsPCAThreads=24
# number of PCs to compute, will likely only use ~10 for 700 samples, so computing 100 should be plenty to play with
numPCs=100

ngsPCAExcludeRegions=ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed
jar=$HOME/ngspca.jar

java -Xmx60G -jar "$jar" \
-input $mosDepthResultsDir \
-outputDir $ngsPCAOutputDir \
-numPC $numPCs \
-sampleEvery 0 \
-threads $ngsPCAThreads \
-iters 	40 \
-randomSeed 42 \
-oversample 100 \
-bedExclude $ngsPCAExcludeRegions

```
### Exclude bed

`ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed` can be found [here](https://github.com/PankratzLab/NGS-PCA/blob/master/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz). This bed file is suitable for analysis of GRCh38/hg38 WGS samples. 

For GRCh38/hg38 WES analysis, the WGS exclude bed file can be concatenated with the bed file that defines the exome targets, where the targets have first been buffered by 20kb. A pre-made WES exclude bed suitable for UKB samples can be found [here](https://github.com/PankratzLab/NGS-PCA/blob/master/resources/GRCh38/UKB_WES/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.xgen.sorted.merge.contig.bed.gz). The original targets used to generate this file are sourced from http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=3801 and can be retrieved with `wget  -nd  biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/xgen_plus_spikein.b38.bed`


### Brief pipeline description

The jar can be downloaded from a release https://github.com/PankratzLab/NGS-PCA/releases or be run from Docker / Apptainer


The ngspca jar will essentially:

1. Select autosomal bins that do not overlap any region in the excluded bed
2. Normalize input data
	- Normalize within sample by computing fold change 
		- log2(coverage of bin / median coverage of all selected bins)
	- Center each bin to median fold-change of 0 across all samples 
3. Perform Randomized PCA
	- Described in https://epubs.siam.org/doi/abs/10.1137/090771806 and https://epubs.siam.org/doi/abs/10.1137/100804139
  	- Similar to the https://github.com/erichson/rSVD R package

## Running with Apptainer (HPC)

The published Docker image can be used directly on HPC systems via [Apptainer](https://apptainer.org/) (formerly Singularity).

**Pull the image:**

```bash
apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest
```

**Run NGS-PCA:**

```bash
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

> **Note:** Bind-mount (`--bind`) any directories containing your input files, output directory, and exclude BED file so they are accessible inside the container.

## Integration Testing

CI runs automatically on every pull request and push to `main`/`master` via [`.github/workflows/integration-test.yml`](.github/workflows/integration-test.yml). It builds the JAR from source and validates deterministic output against precomputed MD5 checksums using 18 1000 Genomes chr1 samples.

### Reproducing the integration test locally

**Prerequisites:** Java 11+, Maven 3.6+ (or use the Docker/Apptainer image instead of building from source).

**1. Clone the repository and build:**

```bash
git clone https://github.com/jlanej/NGS-PCA.git
cd NGS-PCA
mvn -B package --file ngspca/pom.xml
```

**2. Run NGS-PCA on the bundled example data:**

```bash
mkdir -p example/test_output

java -jar ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  -input example/input.files.txt \
  -outputDir example/test_output/ \
  -numPC 2 \
  -sampleEvery 0 \
  -iters 10 \
  -randomSeed 42 \
  -oversample 2 \
  -bedExclude resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

**3. Verify output checksums:**

```bash
cd example
FAILED=0; CHECKED=0
while IFS='  ' read -r expected_md5 filepath; do
  [ -z "$expected_md5" ] && continue
  case "$filepath" in
    ./exampleOutput_1000G_chr1/svd.*)
      CHECKED=$((CHECKED + 1))
      filename=$(basename "$filepath")
      actual_md5=$(md5sum "test_output/$filename" | awk '{print $1}')
      if [ "$actual_md5" = "$expected_md5" ]; then
        echo "PASS: $filename"
      else
        echo "FAIL: $filename (expected: $expected_md5, got: $actual_md5)"
        FAILED=1
      fi
      ;;
  esac
done < 1000G.chr1.md5
echo "$CHECKED files checked; FAILED=$FAILED"
```

Expected output checksums (5 `svd.*.txt` files):

| File | MD5 |
|---|---|
| `svd.bins.txt` | `538a93fe43bc451d915e6672a2e8f357` |
| `svd.loadings.txt` | `4d83db0d1dde4993a11c58d0fcdeb23f` |
| `svd.pcs.txt` | `8a193cc0dc050c9bf527461fb32dc43e` |
| `svd.samples.txt` | `5130c6c74b45485b0869be765245fc68` |
| `svd.singularvalues.txt` | `e09c0dac18f6896e5c00260ec55a933a` |

### Reproducing the integration test with Apptainer

```bash
mkdir -p example/test_output

apptainer run \
  --bind "$(pwd)":/repo \
  ngs-pca.sif \
  -input /repo/example/input.files.txt \
  -outputDir /repo/example/test_output/ \
  -numPC 2 \
  -sampleEvery 0 \
  -iters 10 \
  -randomSeed 42 \
  -oversample 2 \
  -bedExclude /repo/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

Then verify checksums using step 3 above.

