# Benchmark

Wall time and peak RSS for this build against another, on a generated cohort.

**Indicative, not a guide.** Both numbers move with the machine, the JVM, the heap, and above all
the shape of the cohort. Nothing here should be quoted as the speed of NGS-PCA.

```bash
git clone https://github.com/PankratzLab/NGS-PCA.git /tmp/upstream
git -C /tmp/upstream checkout 2ffbcfc
(cd /tmp/upstream/ngspca && mvn -B package -DskipTests)

NUMPC=100 OVERSAMPLE=400 THREADS=8 ./compare.sh \
  /tmp/upstream/ngspca/target/ngspca-0.02-SNAPSHOT.jar \
  ../../ngspca/target/ngspca-0.02-SNAPSHOT.jar 1200 3000
```

`NUMPC + OVERSAMPLE` is the width of the decomposition and matters more than the cohort size: the
QR scales with rows times width squared, and the products with bins times samples times width. A
run left at the small defaults says little about one at `-numPC 500`.

The script refuses to report timings unless both builds produced identical output. A speed
comparison between two different answers is not a comparison.

## A recorded run

Apple M-series laptop, 10 cores, JDK 21, `-Xmx6g`, 1,200 samples × 3,000 bins, width 500,
`-threads 8`:

| build | seconds | peak RSS MB |
|---|---|---|
| upstream `2ffbcfc` | 44 | 1075 |
| this build | 11 | 1002 |

Output identical across all five files.

## What these numbers do not show

**The memory saving is invisible here.** Removing both transposes takes peak from two copies of the
cohort matrix to one, but at this shape that matrix is 29 MB against a gigabyte of JVM and htsjdk
overhead, so it disappears into the noise. The saving is proportional to the matrix — it is worth
157 GB at 140k × 140k and nothing at all here. Only a real cohort will show it.

**Nor does the speed generalise.** 4× at this shape is mostly the QR, which dominates when the
matrix is small relative to the width. At near-square cohort shapes the products dominate instead,
and the speedup there is bounded by how many block-columns the right hand matrix has.

The measurement worth having is a single instrumented production run — `/usr/bin/time -v` on a real
cohort — checked against those two claims.
