# Integration Testing

CI runs automatically on every pull request and push to `main`/`master` via [`.github/workflows/integration-test.yml`](../.github/workflows/integration-test.yml). It builds the JAR from source and validates deterministic output against precomputed MD5 checksums using 18 1000 Genomes chr1 samples. The workflow also runs both UNIFORM and GAUSSIAN distributions and verifies that the resulting comparison summary is reproducible.

## Reproducing the integration test locally

**Prerequisites:** Java 11+, Maven 3.6+, Python 3

**1. Clone the repository and build:**

```bash
git clone https://github.com/jlanej/NGS-PCA.git
cd NGS-PCA
mvn -B package --file ngspca/pom.xml
```

**2. Run NGS-PCA with both distributions:**

```bash
mkdir -p example/test_output example/test_output_gaussian

java -jar ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  -input example/input.files.txt \
  -outputDir example/test_output/ \
  -numPC 2 -sampleEvery 0 -iters 10 -randomSeed 42 -oversample 2 \
  -distribution UNIFORM \
  -bedExclude resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz

java -jar ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  -input example/input.files.txt \
  -outputDir example/test_output_gaussian/ \
  -numPC 2 -sampleEvery 0 -iters 10 -randomSeed 42 -oversample 2 \
  -distribution GAUSSIAN \
  -bedExclude resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

**3. Generate the distribution comparison summary:**

```bash
python3 example/compare_distributions.py \
  example/test_output/ \
  example/test_output_gaussian/ \
  example/test_output/distribution_comparison.txt
```

**4. Verify output checksums:**

```bash
cd example
FAILED=0; CHECKED=0
while IFS='  ' read -r expected_md5 filepath; do
  [ -z "$expected_md5" ] && continue
  case "$filepath" in
    ./exampleOutput_1000G_chr1/svd.*|./exampleOutput_1000G_chr1/distribution_comparison.txt)
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

Expected output checksums (5 `svd.*.txt` files + distribution comparison):

| File | MD5 |
|---|---|
| `svd.bins.txt` | `538a93fe43bc451d915e6672a2e8f357` |
| `svd.loadings.txt` | `4d83db0d1dde4993a11c58d0fcdeb23f` |
| `svd.pcs.txt` | `8a193cc0dc050c9bf527461fb32dc43e` |
| `svd.samples.txt` | `5130c6c74b45485b0869be765245fc68` |
| `svd.singularvalues.txt` | `e09c0dac18f6896e5c00260ec55a933a` |
| `distribution_comparison.txt` | `deb7110ca864a28ebcdd883616cb81bb` |

## Distribution comparison results (18 1000G chr1 samples, 2 PCs)

The comparison confirms that UNIFORM and GAUSSIAN distributions produce effectively identical PCs on these data when power (subspace) iterations are used:

```
PC   UNIFORM_SV   GAUSSIAN_SV   UNIFORM_MEAN  ...  ABS_PEARSON_R  ABS_SPEARMAN_R
PC1  125.538359   125.538359    0.043120      ...  1.000000       1.000000
PC2  48.121044    48.122006    -0.016546      ...  0.999967       1.000000
```

The full summary (including mean, median, and SD for each distribution) is in [`example/exampleOutput_1000G_chr1/distribution_comparison.txt`](../example/exampleOutput_1000G_chr1/distribution_comparison.txt).

## Reproducing the integration test with Apptainer (HPC)

**Pull the image:**

```bash
apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest
```

**Run on the bundled example data:**

```bash
mkdir -p example/test_output example/test_output_gaussian

apptainer run \
  --bind "$(pwd)":/repo \
  ngs-pca.sif \
  -input /repo/example/input.files.txt \
  -outputDir /repo/example/test_output/ \
  -numPC 2 -sampleEvery 0 -iters 10 -randomSeed 42 -oversample 2 \
  -distribution UNIFORM \
  -bedExclude /repo/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz

apptainer run \
  --bind "$(pwd)":/repo \
  ngs-pca.sif \
  -input /repo/example/input.files.txt \
  -outputDir /repo/example/test_output_gaussian/ \
  -numPC 2 -sampleEvery 0 -iters 10 -randomSeed 42 -oversample 2 \
  -distribution GAUSSIAN \
  -bedExclude /repo/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz
```

Then generate the comparison and verify checksums using steps 3–4 above.
