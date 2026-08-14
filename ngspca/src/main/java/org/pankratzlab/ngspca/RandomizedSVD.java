package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.ForkJoinPool;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.random.MersenneTwister;

public class RandomizedSVD {

  /**
   * Distribution used to seed the initial random matrix for the randomized SVD.
   * <p>
   * The Halko et al. algorithm (https://arxiv.org/pdf/0909.4061.pdf) specifies a standard Gaussian
   * random matrix, but in practice a uniform distribution also produces good results, particularly
   * when combined with power (subspace) iterations. The default is UNIFORM for backward
   * compatibility; GAUSSIAN is available for strict adherence to the published algorithm.
   */
  enum DISTRIBUTION {
    UNIFORM, GAUSSIAN;
  }

  //  https://arxiv.org/pdf/0909.4061.pdf
  //  Compute a (truncated) randomized SVD of a BlockRealMatrix
  // Implementation is similar to http://arxiv.org/abs/1608.02148
  private int numComponents;
  static final int DEFAULT_NITERS = 10;
  static final int DEFAULT_OVERSAMPLES = 200;
  private boolean transpose = false;
  private RealMatrix[] rsvd = new RealMatrix[3];
  private final Logger log;

  /**
   * Column names of the original input data
   */
  private final List<String> originalColNames;
  /**
   * Row names of the original input data
   */
  private final List<String> originalRowNames;

  public RandomizedSVD(List<String> originalColNames, List<String> originalRowNames, Logger log) {
    this.originalColNames = originalColNames;
    this.originalRowNames = originalRowNames;
    this.log = log;
  }

  List<String> getColumnNames() {
    return originalColNames;
  }

  List<String> getRowNames() {
    return originalRowNames;
  }

  /**
   * @param A matrix to perform randomized PCA on. Its column order is part of the computation
   *          rather than labelling: the projection below multiplies the columns as they stand by a
   *          random matrix fixed by {@code randomSeed}, so the same data in a different column
   *          order gives results that differ within the approximation error of the method. Callers
   *          are responsible for column order being something that reproduces
   * @param numberOfComponentsToStore number of PCs to compute
   * @param niters specifies the number of power (subspace) iterations to reduce the approximation
   *          error. The power scheme is recommended, if the singular values decay slowly. In
   *          practice, 2 or 3 iterations achieve good results, however, computing power iterations
   *          increases the computational costs
   * @param numOversamples is an oversampling parameter to improve the approximation. A value of at
   *          least 10 is recommended,
   * @param randomSeed random seed for sampling matrix
   * @param d distribution to use for generating the initial random matrix
   * @param threads how many threads the decomposition may use
   */
  public void fit(BlockRealMatrix A, int numberOfComponentsToStore, int niters, int numOversamples,
                  int randomSeed, DISTRIBUTION d, int threads) {
    ForkJoinPool pool = new ForkJoinPool(Math.max(1, threads));
    try {
      fit(A, numberOfComponentsToStore, niters, numOversamples, randomSeed, d, pool);
    } finally {
      pool.shutdown();
    }
  }

  private void fit(BlockRealMatrix A, int numberOfComponentsToStore, int niters, int numOversamples,
                   int randomSeed, DISTRIBUTION d, ForkJoinPool pool) {
    this.numComponents = Math.min(numberOfComponentsToStore,
                                  Math.min(A.getColumnDimension(), A.getRowDimension()));
    if (numComponents < numberOfComponentsToStore) {
      log.info(numberOfComponentsToStore + " PCs requested, but only be able to compute "
               + numComponents);
    }
    log.info("Initializing matrices");

    int n = A.getColumnDimension();
    transpose = A.getRowDimension() < n;
    rsvd[0] = MatrixUtils.createRealMatrix(A.getRowDimension(), numComponents);
    rsvd[1] = MatrixUtils.createRealMatrix(numComponents, 1);
    rsvd[2] = MatrixUtils.createRealMatrix(A.getColumnDimension(), numComponents);

    // The decomposition wants a matrix with at least as many rows as columns, and its transpose.
    // Neither is materialised: every product either one appears in is taken through the matrix as
    // given, using (X^T M)^T = M^T X. At cohort scale each transpose would be a second copy of
    // something measured in hundreds of gigabytes, and taking the product this way is also faster,
    // since dividing the large matrix by column-blocks yields far more tasks than dividing the
    // small one does.
    if (transpose) {
      log.info("Treating the input as transposed, since row N < column N");
      n = A.getRowDimension();
    }

    log.info("Selecting randomized Q using distribution " + d.toString());

    RealMatrix Y = times(A, randn(n, Math.min(n, numComponents + numOversamples), randomSeed, d),
                         pool);

    log.info("Beginning LU decomp iterations");
    for (int i = 0; i < niters; i++) {
      log.info("Subspace iteration: " + Integer.toString(i));
      log.info("Y QR decomp");
      Y = ThinQR.orthonormalBasis(Y, pool);
      log.info("Computing A Y cross prod");
      RealMatrix Z = transposeTimes(A, Y, pool);
      log.info("Z QR decomp");
      Z = ThinQR.orthonormalBasis(Z, pool);
      log.info("A %*% Z");
      Y = times(A, Z, pool);
    }

    RealMatrix Q = ThinQR.orthonormalBasis(Y, pool);
    log.info("Q^T %*% A");
    RealMatrix B = transposeTimes(A, Q, pool).transpose();
    log.info("SVD of reduced matrix");
    SingularValueDecomposition svd = new SingularValueDecomposition(B);

    log.info("Finding W");
    RealMatrix W = Q.multiply(svd.getU());

    log.info("Setting SVD V/W/U results");
    if (transpose) {
      for (int i = 0; i < numComponents; i++) {

        rsvd[0].setColumn(i, svd.getV().getColumn(i));
        rsvd[1].setEntry(i, 0, svd.getSingularValues()[i]);
        rsvd[2].setColumn(i, W.getColumn(i));
      }
    } else {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].setColumn(i, W.getColumn(i));
        rsvd[1].setEntry(i, 0, svd.getSingularValues()[i]);
        rsvd[2].setColumn(i, svd.getV().getColumn(i));
      }

      log.info("Finished SVD");
    }
  }

  /**
   * A X, where A is the input when it has at least as many rows as columns and its transpose
   * otherwise
   */
  private RealMatrix times(BlockRealMatrix input, RealMatrix x, ForkJoinPool pool) {
    return transpose ? transposeProduct(input, x, pool) : ParallelMultiply.multiply(input, x, pool);
  }

  /**
   * A^T X, for the same A
   */
  private RealMatrix transposeTimes(BlockRealMatrix input, RealMatrix x, ForkJoinPool pool) {
    return transpose ? ParallelMultiply.multiply(input, x, pool) : transposeProduct(input, x, pool);
  }

  /**
   * M^T X as (X^T M)^T, so M^T never exists. X is the narrow matrix, so transposing it and the
   * result costs a fraction of what transposing M would.
   */
  private static RealMatrix transposeProduct(BlockRealMatrix m, RealMatrix x, ForkJoinPool pool) {
    // through a BlockRealMatrix regardless of what x is, so the product takes the blocked path it
    // would have taken with M^T on the left
    BlockRealMatrix xt = new BlockRealMatrix(x.transpose().getData());
    return ParallelMultiply.multiply(xt, m, pool).transpose();
  }

  /**
   * @param rows
   * @param columns
   * @param randomSeed seed for deterministic random generation
   * @param d the {@link DISTRIBUTION} to sample from
   * @return a {@link RealMatrix} populated with random values (deterministic random using
   *         {@link MersenneTwister})
   */
  private static RealMatrix randn(int rows, int columns, int randomSeed, DISTRIBUTION d) {
    RealMatrix m = MatrixUtils.createRealMatrix(rows, columns);
    MersenneTwister twister = new MersenneTwister(randomSeed);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m.setEntry(i, j, switch (d) {
          case UNIFORM -> twister.nextDouble();
          case GAUSSIAN -> twister.nextGaussian();
        });
      }
    }
    return m;
  }

  public RealMatrix getV() {
    return (transpose ? rsvd[0] : rsvd[2]);
  }

  public RealMatrix getW() {
    return rsvd[1];
  }

  /**
   * @param file dump the PCs to this text file
   * @param log
   */
  void dumpPCsToText(String file, Logger log) {
    //
    RealMatrix v = rsvd[2];
    v = v.transpose();
    List<String> pcNames = getNumberedColumnHeader("PC", v.getRowDimension());

    dumpMatrix(file, v, "SAMPLE", pcNames, originalColNames, true, log);
  }

  /**
   * @param type name each column carries, e.g. "PC"
   * @param num how many columns to name
   * @return the names, numbered from one
   */
  private static List<String> getNumberedColumnHeader(String type, int num) {
    List<String> names = new ArrayList<>();
    for (int i = 0; i < num; i++) {
      names.add(type + (i + 1));
    }
    return names;
  }

  private static void dumpMatrix(String file, RealMatrix m, String rowTitle,
                                 List<String> outputColumnNames, List<String> outputRowNames,
                                 boolean transposed, Logger log) {
    try (PrintWriter writer = new PrintWriter(new File(file))) {

      StringJoiner joiner = new StringJoiner("\t");
      joiner.add(rowTitle);
      for (String colName : outputColumnNames) {
        joiner.add(colName);
      }
      writer.println(joiner.toString());
      log.info(outputRowNames.size() + " output rows by " + outputColumnNames.size()
               + " output columns");

      for (int outputRow = 0; outputRow < outputRowNames.size(); outputRow++) {
        StringJoiner sample = new StringJoiner("\t");
        sample.add(outputRowNames.get(outputRow));
        for (int outputColumn = 0; outputColumn < outputColumnNames.size(); outputColumn++) {
          if (transposed) {
            sample.add(Double.toString(m.getEntry(outputColumn, outputRow)));
          } else {
            sample.add(Double.toString(m.getEntry(outputRow, outputColumn)));
          }
        }
        writer.println(sample.toString());

      }
    } catch (FileNotFoundException e) {

      log.log(Level.SEVERE, "unable to write to file " + file, e);
      throw new UncheckedIOException("unable to write " + file, e);

    }
  }

  //  private static void printDims(RealMatrix m, Logger log) {
  //    log.info("Row:" + m.getRowDimension() + " Col: " + m.getColumnDimension());
  //  }

  /**
   * @param file loadings will be computed and dumped to this file
   * @param log
   */
  void computeAndDumpLoadings(String file, Logger log) {
    RealMatrix loadingData = rsvd[0];
    List<String> loadingNames = getNumberedColumnHeader("Loading",
                                                            loadingData.getColumnDimension());
    dumpMatrix(file, loadingData, "MARKER", loadingNames, originalRowNames, false, log);

  }

  /**
   * @param file singular values will be dumped to this file
   * @param log
   */
  void dumpSingularValuesToText(String file, Logger log) {
    StringJoiner joiner = new StringJoiner("\n");
    joiner.add("PC\tSINGULAR_VALUES");

    for (int component = 0; component < numComponents; component++) {
      joiner.add(component + 1 + "\t" + Double.toString(rsvd[1].getEntry(component, 0)));
    }

    try {
      FileUtils.writeStringToFile(new File(file), joiner.toString(), Charset.defaultCharset(),
                                  false);
    } catch (IOException e) {
      log.log(Level.SEVERE, "unable to write to file " + file, e);
      throw new UncheckedIOException("unable to write " + file, e);
    }

  }
}
