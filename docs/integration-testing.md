# Integration Testing

CI runs automatically on every pull request and push to `main`/`master` via [`.github/workflows/integration-test.yml`](../.github/workflows/integration-test.yml). It builds the JAR from source and validates deterministic output against precomputed MD5 checksums using 18 1000 Genomes chr1 samples.

## Reproducing the integration test locally

**Prerequisites:** Java 11+, Maven 3.6+

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
| `svd.loadings.txt` | `0af00491e808c5621cb7c441fba9de44` |
| `svd.pcs.txt` | `8b8fdac27bcdd14cb2be1e646d00a8a1` |
| `svd.samples.txt` | `5130c6c74b45485b0869be765245fc68` |
| `svd.singularvalues.txt` | `6a793316e79e460ed38400fe257cba2e` |

## Reproducing the integration test with Apptainer (HPC)

**Pull the image:**

```bash
apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest
```

**Run on the bundled example data:**

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
