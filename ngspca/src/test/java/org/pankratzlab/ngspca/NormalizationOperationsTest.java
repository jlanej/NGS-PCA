package org.pankratzlab.ngspca;

import java.util.logging.Logger;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.BlockRealMatrix;

/**
 * Tests the per-sample medians that normalization computes and now reports
 */
public class NormalizationOperationsTest extends TestCase {

  private static final double DELTA = 1e-12;

  private static final Logger LOG = Logger.getLogger(NormalizationOperationsTest.class.getName());

  private static double log2(double num) {
    return Math.log(num) / Math.log(2);
  }

  public void testColumnMediansOddAndEvenRowCounts() {
    // 3 rows: the median is the middle value
    BlockRealMatrix odd = new BlockRealMatrix(new double[][] {{1, 10}, {2, 40}, {3, 30}});
    double[] oddMedians = NormalizationOperations.columnMedians(odd);
    assertEquals(2, oddMedians.length);
    assertEquals(2.0, oddMedians[0], DELTA);
    assertEquals(30.0, oddMedians[1], DELTA);

    // 4 rows: the median is the mean of the two middle values
    BlockRealMatrix even = new BlockRealMatrix(new double[][] {{1, 10}, {2, 40}, {3, 30}, {4, 20}});
    double[] evenMedians = NormalizationOperations.columnMedians(even);
    assertEquals(2.5, evenMedians[0], DELTA);
    assertEquals(25.0, evenMedians[1], DELTA);
  }

  /**
   * The reported medians must be the ones the fold-change was actually taken against, and the
   * matrix must be normalized exactly as before
   */
  public void testFoldChangeReportsTheMediansItUsed() {
    BlockRealMatrix dm = new BlockRealMatrix(new double[][] {{1, 10}, {2, 40}, {3, 30}});
    double[] medians = NormalizationOperations.foldChangeAndCenterRows(dm, LOG);

    assertEquals(2.0, medians[0], DELTA);
    assertEquals(30.0, medians[1], DELTA);

    // row centering shifts a whole row, so the within-row difference is the fold-change difference
    assertEquals(log2(1.0 / 2.0) - log2(10.0 / 30.0), dm.getEntry(0, 0) - dm.getEntry(0, 1), DELTA);
    assertEquals(log2(2.0 / 2.0) - log2(40.0 / 30.0), dm.getEntry(1, 0) - dm.getEntry(1, 1), DELTA);

    // and every row ends up centered on a median of 0
    for (int row = 0; row < dm.getRowDimension(); row++) {
      assertEquals(0.0, dm.getEntry(row, 0) + dm.getEntry(row, 1), DELTA);
    }
  }

  /**
   * The floor exists so an empty sample does not divide by zero. It must not reach the reported
   * median, which is the sample's actual coverage - flooring is the reader's decision.
   */
  public void testEmptySampleReportsAnUnflooredMedian() {
    BlockRealMatrix dm = new BlockRealMatrix(new double[][] {{0, 4}, {0, 8}, {0, 16}});
    double[] medians = NormalizationOperations.foldChangeAndCenterRows(dm, LOG);

    assertEquals(0.0, medians[0], DELTA);
    assertEquals(8.0, medians[1], DELTA);

    // the floor still applies where it matters: no NaN or infinity from dividing by zero
    for (int row = 0; row < dm.getRowDimension(); row++) {
      for (int col = 0; col < dm.getColumnDimension(); col++) {
        assertTrue("entry " + row + "," + col + " is not finite",
                   !Double.isNaN(dm.getEntry(row, col)) && !Double.isInfinite(dm.getEntry(row, col)));
      }
    }
  }
}
