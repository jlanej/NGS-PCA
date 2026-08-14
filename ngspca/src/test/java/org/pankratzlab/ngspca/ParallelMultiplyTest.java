package org.pankratzlab.ngspca;

import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Dividing the work must not move a single entry, so every case is compared against
 * {@link BlockRealMatrix#multiply} bit-for-bit rather than within a tolerance.
 */
public class ParallelMultiplyTest extends TestCase {

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

  private static void assertMatchesSerial(int m, int n, int k, int threads) {
    BlockRealMatrix a = random(m, n, 17);
    BlockRealMatrix b = random(n, k, 23);
    RealMatrix expected = a.multiply(b);

    ForkJoinPool pool = new ForkJoinPool(threads);
    try {
      RealMatrix actual = ParallelMultiply.multiply(a, b, pool);

      assertEquals(m, actual.getRowDimension());
      assertEquals(k, actual.getColumnDimension());
      for (int row = 0; row < m; row++) {
        for (int column = 0; column < k; column++) {
          assertEquals("entry " + row + "," + column + " with " + threads + " threads",
                       Double.doubleToRawLongBits(expected.getEntry(row, column)),
                       Double.doubleToRawLongBits(actual.getEntry(row, column)));
        }
      }
    } finally {
      pool.shutdown();
    }
  }

  /**
   * Shapes either side of a block boundary, since the slices are cut on it
   */
  public void testMatchesSerialAcrossShapes() {
    assertMatchesSerial(200, 150, 260, 4);   // exactly five block-columns
    assertMatchesSerial(200, 150, 261, 4);   // one column into a sixth
    assertMatchesSerial(97, 53, 209, 4);     // nothing aligned
    assertMatchesSerial(400, 60, 53, 4);     // just past one block-column
  }

  public void testIsUnaffectedByThreadCount() {
    for (int threads : new int[] {2, 3, 8, 17}) {
      assertMatchesSerial(150, 120, 217, threads);
    }
  }

  /**
   * Too few columns to divide, or a pool that cannot divide them - the serial multiply is the
   * answer, and must still be the same answer
   */
  public void testFallsBackWithoutChangingTheResult() {
    assertMatchesSerial(120, 90, 52, 8);
    assertMatchesSerial(120, 90, 300, 1);
  }

  /**
   * The right hand matrix is whatever the previous step returned, which for a small one is an
   * Array2DRowRealMatrix rather than a BlockRealMatrix
   */
  public void testAcceptsANonBlockRightHandMatrix() {
    BlockRealMatrix a = random(300, 40, 5);
    RealMatrix b = MatrixUtils.createRealMatrix(random(40, 60, 6).getData());
    RealMatrix expected = a.multiply(b);

    ForkJoinPool pool = new ForkJoinPool(4);
    try {
      RealMatrix actual = ParallelMultiply.multiply(a, b, pool);
      for (int row = 0; row < 300; row++) {
        for (int column = 0; column < 60; column++) {
          assertEquals(Double.doubleToRawLongBits(expected.getEntry(row, column)),
                       Double.doubleToRawLongBits(actual.getEntry(row, column)));
        }
      }
    } finally {
      pool.shutdown();
    }
  }
}
