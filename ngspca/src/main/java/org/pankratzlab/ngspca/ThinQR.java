package org.pankratzlab.ngspca;

import java.util.concurrent.ForkJoinPool;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Orthonormal basis for the columns of a tall matrix, by Householder QR.
 * <p>
 * This is the power iteration's inner step, and was the whole cost of {@link RandomizedSVD#fit}
 * before it was parallelised. It is the thin Q that is wanted - m by n, not m by m - which rules out
 * {@link org.apache.commons.math3.linear.QRDecomposition}, whose getQ() is m by m and so is not
 * representable at the row counts this is used for.
 * <p>
 * The arithmetic is Jama's, operation for operation, and {@code ThinQRTest} holds it to that
 * bit-for-bit. What differs is only how the work is arranged: columns are stored contiguously
 * rather than as rows of a row-major array, so the inner loops walk memory in order, and the
 * columns updated by each Householder step are independent of each other, so they are updated in
 * parallel. Both leave every value's arithmetic untouched.
 */
class ThinQR {

  private ThinQR() {

  }

  /**
   * @param matrix a matrix with at least as many rows as columns
   * @param pool the column updates run here
   * @return the thin Q of its QR decomposition, with the same dimensions
   */
  static RealMatrix orthonormalBasis(RealMatrix matrix, ForkJoinPool pool) {
    int m = matrix.getRowDimension();
    int n = matrix.getColumnDimension();
    if (m < n) {
      throw new IllegalArgumentException("Cannot take a thin Q of a " + m + " by " + n
                                         + " matrix - it has fewer rows than columns");
    }
    double[][] columns = new double[n][];
    for (int column = 0; column < n; column++) {
      columns[column] = matrix.getColumn(column);
    }
    double[][] q = decompose(columns, m, n, pool);

    // back to row-major, and built through MatrixUtils so the implementation chosen for the result
    // is the one the caller would have got before: it picks by size, and RealMatrix implementations
    // multiply by different algorithms, so the type reaches the values through the next multiply
    double[][] rowMajor = new double[m][n];
    for (int column = 0; column < n; column++) {
      double[] values = q[column];
      for (int row = 0; row < m; row++) {
        rowMajor[row][column] = values[row];
      }
    }
    return MatrixUtils.createRealMatrix(rowMajor);
  }

  /**
   * @param columns the matrix, one column per entry, overwritten with the Householder vectors
   * @param m rows
   * @param n columns
   * @param pool the column updates run here
   * @return the thin Q, one column per entry
   */
  private static double[][] decompose(double[][] columns, int m, int n, ForkJoinPool pool) {
    for (int k = 0; k < n; k++) {
      final double[] pivot = columns[k];
      final int from = k;

      double nrm = 0;
      for (int i = k; i < m; i++) {
        nrm = hypot(nrm, pivot[i]);
      }
      if (nrm == 0.0) {
        continue;
      }
      if (pivot[k] < 0) {
        nrm = -nrm;
      }
      for (int i = k; i < m; i++) {
        pivot[i] /= nrm;
      }
      pivot[k] += 1.0;

      // each column is reflected using the pivot alone, so they do not interact
      inPool(pool, k + 1, n, j -> reflect(pivot, columns[j], from, m));
    }

    double[][] q = new double[n][];
    for (int column = 0; column < n; column++) {
      q[column] = new double[m];
    }
    for (int k = n - 1; k >= 0; k--) {
      final double[] pivot = columns[k];
      final int from = k;
      // q[k] is still the zeros it was allocated with: step k is the first to write it, since
      // earlier steps reflect only columns j >= their own k, all of which are above this one
      q[k][k] = 1.0;
      if (pivot[k] != 0) {
        inPool(pool, k, n, j -> reflect(pivot, q[j], from, m));
      }
    }
    return q;
  }

  /**
   * Run {@code action} over [from, to) in the given pool. Submitting the stream is what confines it
   * there; a bare parallel stream would use the common pool and every core on the machine.
   */
  private static void inPool(ForkJoinPool pool, int from, int to, IntConsumer action) {
    pool.submit(() -> IntStream.range(from, to).parallel().forEach(action)).join();
  }

  /**
   * Apply the Householder reflection held in {@code pivot} to one column
   */
  private static void reflect(double[] pivot, double[] column, int from, int m) {
    double s = 0.0;
    for (int i = from; i < m; i++) {
      s += pivot[i] * column[i];
    }
    s = -s / pivot[from];
    for (int i = from; i < m; i++) {
      column[i] += s * pivot[i];
    }
  }

  /**
   * sqrt(a^2 + b^2) without under/overflow, as Jama computes it - the norm below has to accumulate
   * exactly as it did there
   */
  private static double hypot(double a, double b) {
    double r;
    if (Math.abs(a) > Math.abs(b)) {
      r = b / a;
      r = Math.abs(a) * Math.sqrt(1 + r * r);
    } else if (b != 0) {
      r = a / b;
      r = Math.abs(b) * Math.sqrt(1 + r * r);
    } else {
      r = 0.0;
    }
    return r;
  }
}
