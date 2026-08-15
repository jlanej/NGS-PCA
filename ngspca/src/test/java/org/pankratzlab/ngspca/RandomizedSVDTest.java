package org.pankratzlab.ngspca;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.pankratzlab.ngspca.RandomizedSVD.DISTRIBUTION;

/**
 * The decomposition transposes its input when there are fewer rows than columns, and that branch
 * is not reached by a cohort with more bins than samples - which the example data is. Selected-bin
 * analyses of large cohorts do reach it, so both orientations are exercised here.
 */
public class RandomizedSVDTest extends TestCase {

  private static final Logger LOG = Logger.getLogger(RandomizedSVDTest.class.getName());
  static {
    LOG.setLevel(Level.WARNING);
  }

  private static List<String> names(String prefix, int count) {
    List<String> names = new ArrayList<>();
    for (int i = 0; i < count; i++) {
      names.add(prefix + i);
    }
    return names;
  }

  /**
   * With oversampling that reaches full rank the randomized decomposition should land on the exact
   * singular values, whichever way round the matrix is
   */
  private static void assertRecoversSingularValues(int rows, int columns, int components) {
    Random random = new Random(rows * 31L + columns);
    BlockRealMatrix a = new BlockRealMatrix(rows, columns);
    for (int row = 0; row < rows; row++) {
      for (int column = 0; column < columns; column++) {
        a.setEntry(row, column, random.nextGaussian());
      }
    }
    double[] expected = new SingularValueDecomposition(a).getSingularValues();

    RandomizedSVD svd = new RandomizedSVD(names("sample", columns), names("bin", rows), LOG);
    svd.fit(a.copy(), components, 10, Math.min(rows, columns), 42, DISTRIBUTION.UNIFORM, 4);

    assertEquals(components, svd.getW().getRowDimension());
    for (int component = 0; component < components; component++) {
      double got = svd.getW().getEntry(component, 0);
      assertEquals("singular value " + component + " of a " + rows + " by " + columns + " matrix",
                   expected[component], got, 1e-8 * expected[0]);
    }
  }

  public void testTallMatrix() {
    assertRecoversSingularValues(80, 20, 5);
  }

  /**
   * Fewer rows than columns, so the input is transposed and A_t becomes the matrix as given
   */
  public void testWideMatrix() {
    assertRecoversSingularValues(20, 80, 5);
  }

  public void testSquareMatrix() {
    assertRecoversSingularValues(40, 40, 5);
  }
}
