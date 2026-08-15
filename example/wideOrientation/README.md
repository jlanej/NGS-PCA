# Wide-orientation fixture

A cohort with more samples than bins, which the 1000G example is not.

`RandomizedSVD` transposes its input when there are fewer rows than columns, and then takes `A_t`
from the matrix as given rather than copying it. The 1000G example has 15,593 bins against 18
samples, so neither branch is reached by it — leaving the orientation that a selected-bin analysis
of a large cohort runs in with no coverage. This fixture covers it.

```bash
./generate.sh /tmp/wide 2000 1500
java -jar ../../ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  -input /tmp/wide -outputDir /tmp/wide_output \
  -numPC 20 -iters 5 -oversample 60 -randomSeed 42 -threads 4 -distribution UNIFORM
```

The cohort is generated rather than committed — 2,000 files is a lot to carry, and none of them
needs reading by hand. Depths come from a MINSTD generator formatted by integer arithmetic, so
every awk emits the same bytes.

`expected.md5` holds the checksums of the six output files. Note `-oversample 60` is chosen so
that `numPC + oversample` exceeds `BlockRealMatrix.BLOCK_SIZE`, which is what makes
`ParallelMultiply` divide the products rather than fall back to the serial path.

The job also re-runs with `-threads 1` and requires byte-identical output. That takes the serial
fallback in both `ParallelMultiply` and `ThinQR`, so it compares the parallel paths against the
serial ones — and unlike the checksums it holds on any platform, because it compares two runs on
the same machine rather than against a recorded value.

**These checksums are x86_64, and were taken from a CI run.** They have to be: values that pass
through `Math.log` depend on the platform, and a reference generated anywhere else does not match
the runner. What makes them trustworthy is not where they came from but that this branch is
bit-identical to upstream `PankratzLab@2ffbcfc` in this orientation on a single machine — a
property that does not depend on which machine, so these are the runner's rounding of an output
already checked against upstream. As with `example/1000G.chr1.md5`, the numeric outputs depend on the
platform's `Math.log`, which differs between x86_64 and aarch64 — so they reproduce on CI and on a
cluster, and not on an Apple Silicon workstation. `svd.bins.txt`, `svd.samples.txt` and
`autosomal.median.txt` match everywhere, since none of them involves a transcendental.
