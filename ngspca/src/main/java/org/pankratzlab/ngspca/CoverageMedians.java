package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Reports the per-sample median depth that {@link NormalizationOperations} computes on its way to
 * PCA.
 * <p>
 * The median is a by-product of normalization rather than anything extra: the fold-change step
 * already takes the median of every column, so writing it out costs one small file and no
 * additional pass over the data. Downstream tools that need each sample's median autosomal
 * coverage can read it here instead of re-deriving it from the mosdepth files.
 * <p>
 * The output is a tab-delimited table with a single header line and no comment lines, so it can be
 * read directly by tools that do not skip comments.
 */
class CoverageMedians {

  /**
   * Written when NGS-PCA selected the regions itself, so the rows are known to be autosomal and
   * exclusion-filtered
   */
  static final String AUTOSOMAL_FILE = "autosomal.median.txt";

  /**
   * Written for a user-supplied {@link CmdLine#MATRIX_INPUT_ARG}, whose rows NGS-PCA cannot vouch
   * for
   */
  static final String MATRIX_FILE = "matrix.median.txt";

  /**
   * Median across autosomal bins that survived the exclusion bed - what depth-based association
   * methods normalize against
   */
  static final String AUTOSOMAL_COLUMN = "AUTO_HQ_median";

  /**
   * Median across whatever rows the input matrix held - deliberately not
   * {@link #AUTOSOMAL_COLUMN}, since NGS-PCA did not choose those rows
   */
  static final String MATRIX_COLUMN = "MEDIAN";

  private CoverageMedians() {

  }

  /**
   * @param file write the table here
   * @param valueColumn name of the median column, which states what the rows of the matrix were
   * @param samples sample identifiers, in matrix column order - the same identifiers written to
   *          svd.samples.txt and svd.pcs.txt, so the tables join
   * @param medians per-sample medians, in matrix column order
   * @param numBins number of rows the medians were computed over, reported for every sample so the
   *          table records what it was derived from
   * @param log
   */
  static void write(String file, String valueColumn, List<String> samples, double[] medians,
                    int numBins, Logger log) {
    if (samples.size() != medians.length) {
      throw new IllegalArgumentException("Invalid number of medians (" + medians.length
                                         + ") for " + samples.size() + " samples");
    }
    int flagged = 0;
    try (PrintWriter writer = new PrintWriter(new File(file))) {
      writer.println("SAMPLE\t" + valueColumn + "\tN_BINS");
      for (int i = 0; i < samples.size(); i++) {
        // catches NaN as well as zero and negative medians
        if (!(medians[i] > NormalizationOperations.MIN_DEPTH)) {
          flagged++;
        }
        writer.println(samples.get(i) + "\t" + Double.toString(medians[i]) + "\t"
                       + Integer.toString(numBins));
      }
    } catch (FileNotFoundException e) {
      log.log(Level.SEVERE, "unable to write to file " + file, e);
      throw new UncheckedIOException("unable to write " + file, e);
    }
    log.info("Wrote " + samples.size() + " per-sample medians, computed over " + numBins
             + " bins, to " + file);
    if (flagged > 0) {
      log.warning(flagged + " sample(s) have a median depth at or below "
                  + Double.toString(NormalizationOperations.MIN_DEPTH)
                  + " in " + file
                  + " - these are empty or failed inputs, and any normalization against them is meaningless");
    }
  }

}
