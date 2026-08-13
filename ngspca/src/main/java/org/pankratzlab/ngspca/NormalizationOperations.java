package org.pankratzlab.ngspca;

import java.util.logging.Logger;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.ranking.NaNStrategy;

/**
 * class to perform input matrix normalization prior to PCA
 */
public class NormalizationOperations {

  /**
   * 0.005 is half the lowest value given by mosdepth
   */
  static final double MIN_DEPTH = 0.005;

  private NormalizationOperations() {

  }

  /**
   * compute fold-change (by column) , and then center the matrix so each row has median of 0;
   *
   * @param m an {@link RealMatrix} that has been FC-ed by column and centered by row
   * @return the median of each column, as computed from the input values and prior to the
   *         {@link #MIN_DEPTH} floor - this is each sample's median depth across the rows of the
   *         matrix, and the value the fold-change was taken against
   */
  static double[] foldChangeAndCenterRows(RealMatrix dm, Logger log) {
    // compute fold change
    double[] medians = computeFoldChangeByColumn(dm, log);
    // center rows to median of 0
    centerRowsToMedian(dm);
    return medians;
  }

  /**
   * Set the values of the matrix to the log 2 fold change (computed by column)
   *
   * @param dm the {@link RealMatrix} that will be converted
   * @return the un-floored median of each column
   */
  private static double[] computeFoldChangeByColumn(RealMatrix dm, Logger log) {
    double[] medians = columnMedians(dm);

    // the value actually divided by - a sample with no coverage would otherwise divide by zero
    double[] denominators = new double[medians.length];
    for (int column = 0; column < medians.length; column++) {
      denominators[column] = Math.max(medians[column], MIN_DEPTH);
    }

    // convert columns to log2 fold-change from median
    for (int row = 0; row < dm.getRowDimension(); row++) {
      for (int column = 0; column < dm.getColumnDimension(); column++) {
        double entry = dm.getEntry(row, column);
        double standard = log2(Math.max(entry, MIN_DEPTH) / denominators[column]);
        if (Double.isNaN(standard)) {
          throw new IllegalArgumentException("Invalid sample normalized value ("
                                             + Double.toString(Double.NaN) + ") detected");
        }
        dm.setEntry(row, column, standard);
      }
    }
    return medians;
  }

  /**
   * @param dm the {@link RealMatrix} to summarize, holding one sample per column
   * @return the median of each column, un-floored
   */
  static double[] columnMedians(RealMatrix dm) {
    double[] medians = new double[dm.getColumnDimension()];
    for (int column = 0; column < dm.getColumnDimension(); column++) {
      medians[column] = median(dm.getColumn(column));
    }
    return medians;
  }

  private static double median(double[] tmp) {
    return new Median().withNaNStrategy(NaNStrategy.REMOVED).evaluate(tmp);
  }

  private static double log2(double num) {
    return Math.log(num) / Math.log(2);
  }

  /**
   * @param dm Center the rows of this {@link RealMatrix} to a median of 0
   */
  private static void centerRowsToMedian(RealMatrix dm) {
    for (int row = 0; row < dm.getRowDimension(); row++) {
      double[] tmp = new double[dm.getColumnDimension()];
      for (int col = 0; col < dm.getColumnDimension(); col++) {
        tmp[col] = dm.getEntry(row, col);
      }
      double median = median(tmp);
      for (int col = 0; col < dm.getColumnDimension(); col++) {
        dm.setEntry(row, col, tmp[col] - median);
      }
    }
  }

}
