# Future Work / TODOs

This document summarizes known areas for improvement in NGS-PCA, covering performance, scientific correctness, and code quality. Contributions are welcome.

**Where the time goes.** Measured on a synthetic cohort of 200,000 bins by 200 samples, 20 PCs,
200 oversamples, 10 power iterations, 8 threads — before and after the decomposition was
parallelised:

| phase | before | after |
|---|---|---|
| load mosdepth files | 7 s | 7 s |
| serialize raw matrix (gzip) | 14 s | 14 s |
| normalize | 14 s | 14 s |
| power iterations | 2540 s | 50 s |
| finish (QᵀA, SVD, write) | — | 7 s |
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

## 2. Clarify median vs. mean centering in normalization (scientific)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/NormalizationOperations.java`

The row-centering step uses the **median** rather than the **mean**:

```java
private static void centerRowsToMedian(RealMatrix dm) { ... }
```

Standard PCA is defined in terms of deviations from the mean, so the resulting components do not correspond to the eigenvectors of the sample covariance matrix. The median-centering approach is a form of robust preprocessing that reduces the influence of coverage outliers, but this diverges from classical PCA.

**Options:**
- Switch to mean centering for standard, statistically interpretable PCA.
- Keep median centering but clearly document that this is a robust variant and that the singular vectors do not strictly align with the covariance eigenvectors.

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

## 4. Avoid deep-copy transposition for large matrices (memory)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/RandomizedSVD.java`

The algorithm caches an explicit transposed copy of the input matrix:

```java
BlockRealMatrix A_t = A.transpose();
```

Apache Commons Math's `transpose()` physically allocates and copies all data into a new matrix. For large inputs (e.g., 10⁶ bins × 2000 samples), this roughly doubles RAM usage.

**It is worse than doubling, and this is now the ceiling that binds first.** `fit` may hold three
full copies at once: the caller's matrix, which `computeSVD` still references throughout; the
shape transpose taken when there are fewer rows than columns; and this cached `A_t`. At 140k by
140k that is roughly 471 GB for a 157 GB matrix. Since the decomposition was parallelised, a
cohort will run out of memory before it runs out of time.

**Recommended fix:** a lazy transpose wrapper (`RealMatrix` subclass) translating row/column
indices without copying. Note it must keep `BlockRealMatrix.multiply`'s blocked path reachable, or
the products fall back to the generic one and every committed checksum moves — the parity tests in
`ParallelMultiplyTest` are the gate for any attempt.

---

## 4b. Parallelism ceiling in the matrix products (performance)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/ParallelMultiply.java`

The products are divided by output block-column, which bounds them at `ceil(k / 52)` concurrent
tasks, where `k` is `-numPC` plus `-oversample` — ten for `-numPC 500 -oversample 0`. Beyond that
a node's remaining cores go unused during the phase that dominates a large run.

Dividing by rows of the left-hand matrix instead would lift the ceiling, but `BlockRealMatrix`
keeps `blocks`, `blockRows` and `blockColumns` private, so it needs either reflection or a
row-slice that does not copy — see item 4. Raising `-oversample` also raises the ceiling, at the
cost of a wider decomposition.

`QᵀA` after the loop is still serial. It is one product against the full matrix where the loop
runs twenty, and parallelising it means slicing the large matrix rather than the small one.

---

## 5. I/O Memory Bloat: HTSJDK BEDFeature Instantiation (performance)

**Files:** `ngspca/src/main/java/org/pankratzlab/ngspca/MosdepthUtils.java`, `ngspca/src/main/java/org/pankratzlab/ngspca/BedUtils.java`

Coverage data from mosdepth output files is currently parsed using HTSJDK's `BEDFileReader`:

```java
BEDFileReader reader = new BEDFileReader(file, false);
CloseableIterator<BEDFeature> iter = reader.iterator();
List<BEDFeature> result = iter.stream()
    .filter(bf -> ucscRegions.contains(getBedUCSC(bf)))
    .collect(Collectors.toList());
```

A whole-genome BAM processed into 1 kb bins produces roughly 3 million lines. HTSJDK creates a heavy `BEDFeature` object for every single line. Processing a cohort with `threads = 4` instantiates ~12 million `BEDFeature` objects at any given moment just to extract the 4th column (coverage), generating immense garbage-collection (GC) pressure.

**Recommended fix:** Because mosdepth outputs predictable, tab-delimited files, drop HTSJDK for the coverage extraction phase. Read the file iteratively with standard `BufferedReader` or `Files.lines()`, split by `\t`, and directly parse the 4th element as a double. This completely bypasses the creation of intermediate objects.

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

