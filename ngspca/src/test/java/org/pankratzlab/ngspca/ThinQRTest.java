package org.pankratzlab.ngspca;

import java.util.Random;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import Jama.Matrix;
import Jama.QRDecomposition;

/**
 * {@link ThinQR} rearranges how Jama's QR is computed but not what it computes. Jama is kept as a
 * test dependency to hold it to that: these compare the two bit-for-bit, since anything less would
 * let the PCs drift.
 */
public class ThinQRTest extends TestCase {

  private static BlockRealMatrix random(int m, int n, long seed) {
    Random random = new Random(seed);
    BlockRealMatrix matrix = new BlockRealMatrix(m, n);
    for (int row = 0; row < m; row++) {
      for (int column = 0; column < n; column++) {
        matrix.setEntry(row, column, random.nextGaussian());
      }
    }
    return matrix;
  }

  private static void assertMatchesJama(int m, int n, long seed) {
    BlockRealMatrix matrix = random(m, n, seed);
    double[][] expected = new QRDecomposition(new Matrix(matrix.getData())).getQ().getArray();

    RealMatrix actual = ThinQR.orthonormalBasis(matrix);

    assertEquals(m, actual.getRowDimension());
    assertEquals(n, actual.getColumnDimension());
    for (int row = 0; row < m; row++) {
      for (int column = 0; column < n; column++) {
        assertEquals("entry " + row + "," + column + " of a " + m + " by " + n + " matrix",
                     Double.doubleToRawLongBits(expected[row][column]),
                     Double.doubleToRawLongBits(actual.getEntry(row, column)));
      }
    }
  }

  public void testMatchesJamaOnTallMatrices() {
    assertMatchesJama(500, 20, 42);
    assertMatchesJama(97, 13, 7);
  }

  /**
   * The power iteration reaches a square matrix whenever the requested components plus oversampling
   * meet the sample count
   */
  public void testMatchesJamaOnASquareMatrix() {
    assertMatchesJama(40, 40, 11);
  }

  public void testMatchesJamaOnASingleColumn() {
    assertMatchesJama(50, 1, 3);
  }

  /**
   * A duplicated column makes a Householder norm zero, which is the branch the decomposition skips
   */
  public void testMatchesJamaWithARankDeficientColumn() {
    BlockRealMatrix matrix = random(60, 5, 5);
    matrix.setColumn(3, matrix.getColumn(1));
    matrix.setColumn(4, new double[60]);
    double[][] expected = new QRDecomposition(new Matrix(matrix.getData())).getQ().getArray();

    RealMatrix actual = ThinQR.orthonormalBasis(matrix);

    for (int row = 0; row < 60; row++) {
      for (int column = 0; column < 5; column++) {
        assertEquals(Double.doubleToRawLongBits(expected[row][column]),
                     Double.doubleToRawLongBits(actual.getEntry(row, column)));
      }
    }
  }

  /**
   * A thin Q of a wide matrix is not defined, and silently returning something of the wrong shape
   * would be discovered much later
   */
  public void testRefusesAWideMatrix() {
    try {
      ThinQR.orthonormalBasis(random(4, 9, 1));
      fail("expected an IllegalArgumentException for a matrix with fewer rows than columns");
    } catch (IllegalArgumentException expected) {
      // expected
    }
  }
}
