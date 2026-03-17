# Future Work / TODOs

This document summarizes known areas for improvement in NGS-PCA, covering performance, scientific correctness, and code quality. Contributions are welcome.

---

## 1. Eliminate Jama dependency in Randomized SVD (performance)

**File:** `ngspca/src/main/java/org/pankratzlab/ngspca/RandomizedSVD.java`

The subspace iteration loop currently converts between Apache Commons Math `RealMatrix` objects and Jama `Matrix` objects in order to perform QR decomposition:

```java
QRDecomposition qr = new QRDecomposition(new Matrix(Y.getData()));
Y = MatrixUtils.createRealMatrix(qr.getQ().getArray());
```

Each conversion allocates a full `double[][]` copy of the matrix. With large datasets (e.g., 10⁶ genomic bins × 300-column subspace) and many iterations (10–40), this causes repeated GC pressure and memory bloat.

**Recommended fix:** Use Apache Commons Math's own `QRDecomposition` directly, avoiding all intermediate array copies:

```java
org.apache.commons.math3.linear.QRDecomposition qr =
    new org.apache.commons.math3.linear.QRDecomposition(Y);
Y = qr.getQ();
```

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

**Recommended fix:** If memory is a concern, compute `A.transpose().multiply(...)` inline and rely on any library-level optimization, or implement a lazy transpose wrapper (`RealMatrix` subclass) that translates row/column indices without copying data.

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

