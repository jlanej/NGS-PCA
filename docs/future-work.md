# Future Work / TODOs

This document summarizes known areas for improvement in NGS-PCA, covering performance, scientific correctness, and code quality. Contributions are welcome.

**Where the time goes.** Measured on a synthetic cohort of 200,000 bins by 200 samples, 20 PCs,
200 oversamples, 10 power iterations, 8 threads — before and after the decomposition was
parallelised. The before run logged no separate finish phase; that cell is the remainder of its
total:

| phase | before | after |
|---|---|---|
| load mosdepth files | 7 s | 7 s |
| serialize raw matrix (gzip) | 14 s | 14 s |
| normalize | 14 s | 14 s |
| power iterations | 2540 s | 50 s |
| finish (QᵀA, SVD, write) | 324 s | 7 s |
| **total** | **2899 s** | **95 s** |

Every output byte-identical between the two. The decomposition is no longer the whole cost, so
serialization and normalization — items 3 and 5 below — are now roughly 15% each rather than the
rounding errors they were, and are the next things worth measuring. Note that this cohort is
tall and thin; at the near-square shapes of a large cohort the matrix products dominate further
and the QR matters less.

---

## 1. Eliminate Jama dependency in Randomized SVD — DONE

Jama is gone from the shipped jar; it remains a test-scope dependency, as the reference
`ThinQR` is held to bit-for-bit.

**Not by the route this item recommended.** Commons Math's `QRDecomposition.getQ()` returns an
m by m matrix, and the Q wanted here is the thin one — at 140,000 bins the square Q is not
representable, let alone cheap. `ThinQR` performs Jama's Householder arithmetic, operation for
operation, over column-major storage with the independent per-column updates run in parallel.
Measured at 140000 x 500: 8.2 s against Jama's 992.6 s, every one of 70,000,000 entries identical.

The array copies this item blamed were not the cost — they measured at a tenth of a second.
Jama's row-major layout was, because each inner loop strides across as many separate arrays as
the matrix has rows, and that penalty grows with row count rather than staying proportional.

---

## 2. Clarify median vs. mean centering in normalization (scientific) — DONE

Documented rather than changed: the README's overview now states that median centering makes this
a robust variant of PCA whose components are not exactly the covariance eigenvectors, and suggests
wording for methods text. Switching to mean centering would change every output and break parity
with upstream for no clear scientific gain; if classical PCA is ever wanted, it belongs behind a
flag.

---

## 3. Reduce per-row array allocation in normalization (performance)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/NormalizationOperations.java`

Inside `centerRowsToMedian`, a new `double[]` array is allocated for every row:

```java
for (int row = 0; row < dm.getRowDimension(); row++) {
  double[] tmp = new double[dm.getColumnDimension()];
  ...
}
```

For datasets with millions of genomic bins, this creates millions of short-lived heap objects.

**Recommended fix:** Hoist the array allocation outside the loop and reuse it:

```java
double[] tmp = new double[dm.getColumnDimension()];
for (int row = 0; row < dm.getRowDimension(); row++) { ... }
```

Additionally, consider using Commons Math matrix visitor patterns (`RealMatrixChangingVisitor`) for more cache-friendly block-level iteration.

---

## 4. Avoid deep-copy transposition for large matrices — DONE

Neither transpose is materialised any more. `fit` needs a matrix with at least as many rows as
columns and also its transpose; it now holds only the matrix as given, and takes every product
either one appears in through it, using `Mᵀ X = (Xᵀ M)ᵀ`. Only the narrow matrix and the result are
transposed, and both are a fraction of the size.

Peak drops to a single copy of the matrix in both orientations. With more bins than samples that
removes one copy, the cached `A_t` — roughly 314 GB to 157 GB for a 140k by 140k cohort. With
more samples than bins it removes two, since the shape transpose was materialised first and `A_t`
cached from it. `example/benchmark/README.md` measures both.

It is also faster, which was not the point but is the larger effect: 2.0 s against 3.7 s for one
product at 8000 by 8000. Dividing the large matrix by column-blocks yields hundreds of tasks where
dividing the narrow one yielded ten, so the ceiling in item 4b does not bind on this product.

The identity holds to the bit — multiplication commutes exactly in IEEE 754 and both forms sum over
the same shared index in the same block order — and is pinned by
`ParallelMultiplyTest.testTransposedProductMatchesAMaterialisedTranspose`. Verified end to end in
both orientations, including against upstream `PankratzLab@2ffbcfc`.

---

## 4b. Parallelism ceiling in the matrix products (performance)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/ParallelMultiply.java`

The products are divided by output block-column, so concurrency follows the width of the
right-hand matrix. Each subspace iteration multiplies once in each direction: the product that
takes the narrow matrix on the right is bounded at `ceil(k / 52)` tasks, where `k` is `-numPC`
plus `-oversample` — five at the defaults, ten for `-numPC 500 -oversample 0` — while the other
slices the cohort matrix and follows a cohort dimension, enough for any node. The two directions
cost the same FLOPs, so roughly half the product work runs at the ceiling, and beyond `k / 52`
threads the remaining cores are idle during that half. The one-off products — `A Ω` before the
loop, `Qᵀ A` after it — split the same way, one on each side, with orientation deciding which.

Dividing by rows of the left-hand matrix instead would lift the ceiling, but `BlockRealMatrix`
keeps `blocks`, `blockRows` and `blockColumns` private, so it needs either reflection or a
row-slice that does not copy — see item 4. Raising `-oversample` also raises the ceiling, at the
cost of a wider decomposition.

---

## 5. I/O Memory Bloat: HTSJDK BEDFeature Instantiation (performance) — DONE

The coverage extraction no longer goes through htsjdk: `loadSpecificRegions` reads the gzip
stream itself and returns one `double[]` per file rather than a `List` of three million
`BEDFeature`s. `BedUtilsTest` pins the replacement to what `BEDCodec` produced — the 1-based
keys, the fourth column, the skipped lines, bgzf content, and the refusal of content that is not
gzip. Region selection from the first file still uses htsjdk, so the keys the filter matches
against remain the codec's own.

Measured at 150 samples by 200,000 bins, 8 threads, `-Xmx4g`: loading drops from 7 s to 2 s and
peak RSS from 4.0 GB to 1.6 GB, with every output byte-identical, serialized matrices included.
The in-flight cost per queued file is now the bin count times eight bytes, so `-threads 24` no
longer implies tens of gigabytes of feature objects at cohort scale.

---

## 6. Silent NullPointerException in FileOps.gzLines() (correctness)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/FileOps.java`

The `gzLines()` method catches an `IOException` if a file is missing or corrupted, but then falls through to use the potentially-null `gzipIs` stream:

```java
} catch (IOException e) {
  log.log(Level.SEVERE, "an exception was thrown while reading " + path.toString(), e);
  closeSafely(gzipIs, log);
  // ...
}
BufferedReader reader = new BufferedReader(new InputStreamReader(gzipIs)); // NPE if gzipIs is null
```

If `IOException` is triggered, `gzipIs` is `null`. The code catches the exception but doesn't abort, proceeding directly to `new InputStreamReader(gzipIs)` and throwing a raw `NullPointerException` that destroys the thread and provides a confusing stack trace for the end-user.

**Recommended fix:** If an `IOException` is caught, either throw an `UncheckedIOException` or `return Stream.empty()` immediately.

---

## 7. Thread Pool Starvation Risk (performance/correctness)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/MosdepthUtils.java`

The producer task is submitted to the same fixed thread pool it feeds:

```java
ExecutorService executor = Executors.newFixedThreadPool(Math.max(threads, 2));
Runnable producerTask = () -> { ... };
executor.submit(producerTask);
```

The producer task occupies a thread in the pool indefinitely while feeding the blocking queue. With `-threads 2`, one thread is locked as the producer, leaving only one thread for actual file-reading work. In more complex scenarios, mixing producers and consumers in the same small, fixed thread pool can result in deadlocks if the queue dynamics back up.

**Recommended fix:** Run the producer on a separate, dedicated thread (`new Thread(producerTask).start()`) so the entire `ExecutorService` is 100% dedicated to parsing files.

---

## 8. Fragile Matrix Row Assumption (correctness)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/MosdepthUtils.java`

The `setColumnData()` method assumes the filtered regions from every sample's file are in the exact same order:

```java
for (int row = 0; row < features.size(); row++) {
  dm.addToEntry(row, col, Double.parseDouble(features.get(row).getName()));
}
```

The `ucscRegions` variable is a `HashSet<String>`. While reading files iteratively preserves their top-to-bottom order, the code blindly trusts that every mosdepth file has the exact same bins listed in the exact same order. If one sample has a missing region line, or regions are sorted differently, coverage will be silently assigned to the wrong genomic bins in the final matrix.

**Recommended fix:** Map bin coordinates to a matrix row index during initialization (e.g., `Map<String, Integer> regionToRowMap`). When reading the coverage value for a bin, look up its exact row index to guarantee mathematical alignment across the entire cohort matrix.

