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

## A recorded sweep

`sweep.sh` crosses the two axes that decide which regime a run is in. Apple M-series laptop,
10 cores, JDK 21, `-Xmx4g`, `-threads 8`, against upstream `2ffbcfc`. Output identical to upstream
on every row.

| shape | samples x bins | width | upstream s | this build s | upstream MB | this build MB | speedup |
|---|---|---|---|---|---|---|---|
| bins>samples | 1200 x 3000 | 500 | 43 | 11 | 802 | 776 | 3.9x |
| samples>bins | 3000 x 1200 | 500 | 43 | 11 | 839 | 774 | 3.9x |
| near-square | 2000 x 2000 | 500 | 44 | 13 | 805 | 775 | 3.4x |
| samples>bins | 3000 x 1200 | 100 | 7 | 5 | 754 | 748 | 1.4x |
| bins>samples | 1200 x 3000 | 100 | 7 | 5 | 745 | 748 | 1.4x |

**Width drives the speedup; orientation does not.** 3.9x at width 500 in both orientations, 1.4x at
width 100 in both. That is what the arithmetic says should happen: the QR costs (rows + columns)
times width squared, and rows + columns does not change when the matrix is transposed, so neither
does the cost. Widening the decomposition makes the QR a larger share of the run, and the QR is
what got faster. A run at `-numPC 500` is in the top group; one at `-numPC 20` is in the bottom.

**The memory column matches the copy count, at this scale.** Removing both transposes takes peak
from several copies of the matrix to one, and how many were removed depends on the orientation:

| shape | matrix | copies removed | predicted | observed |
|---|---|---|---|---|
| bins>samples | 28.8 MB | 1 (upstream held A and A_t) | 28.8 MB | 26 MB |
| samples>bins | 28.8 MB | 2 (upstream also held the shape transpose) | 57.6 MB | 65 MB |
| near-square | 32.0 MB | 1 | 32.0 MB | 30 MB |

Three shapes, each within about 10% of prediction, and twice the saving in the orientation where
upstream held three copies rather than two. That is a check on the model behind the 314 GB to
157 GB claim, not a measurement of it - the saving is proportional to the matrix, and here the
matrix is 29 MB against 750 MB of JVM and htsjdk baseline. The two narrow rows show nothing, which
is expected: those runs last five seconds and peak RSS never gets far from the baseline.

## What these numbers do not show

**They do not generalise by size.** The largest matrix here is 32 MB; a production one is 157 GB. Nothing in
this table constrains what happens when the products stop fitting in cache, and at cohort shapes
the products dominate the QR, with a speedup bounded by how many block-columns the right hand
matrix has rather than by how many cores were asked for.

The measurement worth having is a single instrumented production run — `/usr/bin/time -v` on a real
cohort — checked against those two claims.
