package org.pankratzlab.ngspca;

import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * Tests the per-sample median coverage table
 */
public class CoverageMediansTest extends TestCase {

  private static final Logger LOG = Logger.getLogger(CoverageMediansTest.class.getName());

  private static List<String> writeAndRead(List<String> samples, double[] medians,
                                           int numBins) throws Exception {
    File out = File.createTempFile("ngspca.medians", ".txt");
    out.deleteOnExit();
    CoverageMedians.write(out.getAbsolutePath(), CoverageMedians.AUTOSOMAL_COLUMN, samples, medians,
                          numBins, LOG);
    return Files.readAllLines(out.toPath(), Charset.defaultCharset());
  }

  /**
   * A single header line, no comments, and one row per sample in matrix column order - what a
   * downstream reader selecting SAMPLE and AUTO_HQ_median expects to find
   */
  public void testWritesAReadableTable() throws Exception {
    List<String> lines = writeAndRead(Arrays.asList("sampleA", "sampleB"),
                                      new double[] {30.5, 0.25}, 1234);

    assertEquals(3, lines.size());
    assertEquals("SAMPLE\t" + CoverageMedians.AUTOSOMAL_COLUMN + "\tN_BINS", lines.get(0));
    assertEquals("sampleA\t30.5\t1234", lines.get(1));
    assertEquals("sampleB\t0.25\t1234", lines.get(2));
  }

  /**
   * A sample with no coverage is reported rather than dropped or floored, so the reader can see it
   */
  public void testEmptySampleIsReported() throws Exception {
    List<String> lines = writeAndRead(Arrays.asList("empty"), new double[] {0.0}, 10);

    assertEquals(2, lines.size());
    assertEquals("empty\t0.0\t10", lines.get(1));
  }

  /**
   * The table is one of the run's outputs, so failing to write it is a failed run rather than a
   * missing file for a reader to trip over
   */
  public void testFailsLoudlyWhenItCannotBeWritten() throws Exception {
    File dir = java.nio.file.Files.createTempDirectory("ngspca.unwritable").toFile();
    dir.deleteOnExit();
    try {
      CoverageMedians.write(dir.getAbsolutePath(), CoverageMedians.AUTOSOMAL_COLUMN,
                            Arrays.asList("sampleA"), new double[] {1.0}, 10, LOG);
      fail("expected an UncheckedIOException when the file cannot be written");
    } catch (java.io.UncheckedIOException expected) {
      // expected
    }
  }

  /**
   * Sample names and medians are matched by position, so a length mismatch is a bug rather than
   * something to write out and let a reader discover
   */
  public void testRefusesMismatchedLengths() {
    try {
      writeAndRead(Arrays.asList("sampleA", "sampleB"), new double[] {30.5}, 10);
      fail("expected an IllegalArgumentException for 2 samples and 1 median");
    } catch (IllegalArgumentException expected) {
      // expected
    } catch (Exception e) {
      fail("expected an IllegalArgumentException, got " + e);
    }
  }
}
