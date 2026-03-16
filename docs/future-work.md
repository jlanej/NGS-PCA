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


