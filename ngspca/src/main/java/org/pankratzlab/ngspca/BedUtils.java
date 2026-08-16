package org.pankratzlab.ngspca;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

public class BedUtils {

  private BedUtils() {

  }

  /**
   * @param file load autosomal {@link BEDFeature}s from this file and convert to UCSC format
   * @param excluder autosomal regions that return true for
   *          {@link BEDOverlapDetector#overlapsAny(Locatable)} will not be included
   * @return {@link List} of ucsc regions
   */
  static List<String> loadAutosomalUCSC(String file, BEDOverlapDetector excluder) {
    return BedUtils.loadAll(file).stream().filter(BedUtils::autosomal)
                   .filter(excluder::overlapsNone).map(BedUtils::getBedUCSC)
                   .collect(Collectors.toList());
  }

  /**
   * @param file load all {@link BEDFeature}s in this file
   * @return {@link List} of {@link BEDFeature}s
   */
  static List<BEDFeature> loadAll(String file) {
    BEDFileReader reader = new BEDFileReader(file, false);
    List<BEDFeature> result = reader.iterator().stream().collect(Collectors.toList());
    reader.close();
    return result;

  }

  private static final int BUFFER_BYTES = 1 << 16;

  /**
   * Reads the file directly rather than through htsjdk: a cohort is up to hundreds of thousands of
   * these files, and {@link BEDCodec} builds a {@link BEDFeature} per line only for the fourth
   * column to be read off it. Equivalence with what the codec produced is pinned by
   * {@code BedUtilsTest}: keys are contig:start+1-end, since tribble is 1-based inclusive where
   * bed is 0-based half-open; the coverage is the fourth column exactly; and blank, {@code #},
   * {@code track} and {@code browser} lines are skipped. The region list this filters against is
   * still built by htsjdk from the first file, so a drift in key format fails the consumer's count
   * check rather than misassigning rows.
   *
   * @param file mosdepth bed file to read
   * @param ucscRegions only rows for these regions are kept
   * @return the coverage column this file holds, in file order
   */
  static BedRegionResult loadSpecificRegions(String file, Set<String> ucscRegions) {
    double[] coverage = new double[ucscRegions.size()];
    int matched = 0;
    try (BufferedReader reader = open(file)) {
      String line;
      while ((line = reader.readLine()) != null) {
        if (line.isEmpty() || line.charAt(0) == '#' || line.startsWith("track")
            || line.startsWith("browser")) {
          continue;
        }
        int t1 = line.indexOf('\t');
        int t2 = line.indexOf('\t', t1 + 1);
        int t3 = line.indexOf('\t', t2 + 1);
        if (t1 < 0 || t2 < 0 || t3 < 0) {
          throw new IllegalArgumentException("Invalid bed line in " + file + ": " + line);
        }
        String key;
        try {
          // parsed and re-printed rather than substringed, exactly as the codec normalised them
          int start = Integer.parseInt(line, t1 + 1, t2, 10) + 1;
          int end = Integer.parseInt(line, t2 + 1, t3, 10);
          key = line.substring(0, t1) + ':' + start + '-' + end;
        } catch (NumberFormatException nfe) {
          throw new IllegalArgumentException("Invalid bed line in " + file + ": " + line, nfe);
        }
        if (ucscRegions.contains(key)) {
          if (matched < coverage.length) {
            int t4 = line.indexOf('\t', t3 + 1);
            String value = t4 < 0 ? line.substring(t3 + 1) : line.substring(t3 + 1, t4);
            try {
              coverage[matched] = Double.parseDouble(value);
            } catch (NumberFormatException nfe) {
              throw new IllegalArgumentException("Invalid (non-numeric) coverage value in file "
                                                 + file + " in row " + matched, nfe);
            }
          }
          matched++;
        }
      }
    } catch (IOException e) {
      throw new UncheckedIOException("unable to read " + file, e);
    }
    return new BedRegionResult(file, coverage, matched);
  }

  /**
   * Buffered reader over the gzip file. Content that is not gzip fails here, as it did through
   * htsjdk, which trusts the extension the same way.
   */
  private static BufferedReader open(String file) throws IOException {
    InputStream in = new GZIPInputStream(new BufferedInputStream(new FileInputStream(file),
                                                                 BUFFER_BYTES),
                                         BUFFER_BYTES);
    return new BufferedReader(new InputStreamReader(in, StandardCharsets.ISO_8859_1),
                              BUFFER_BYTES);
  }

  /**
   * The coverage column one bed file yielded, in file order, alongside the file it came from so a
   * consumer taking results off a queue can check it got the one it expected. {@code matched}
   * counts every line whose region was asked for, including any beyond the expected count, so the
   * consumer's count check sees duplicated regions rather than a silently truncated column.
   */
  static record BedRegionResult(String file, double[] coverage, int matched) {

  }

  /**
   * @param bedFeature
   * @return UCSC representation of this {@link BEDFeature}
   */
  static String getBedUCSC(BEDFeature bedFeature) {
    return bedFeature.getContig() + ":" + bedFeature.getStart() + "-" + bedFeature.getEnd();

  }

  private static boolean autosomal(BEDFeature bedFeature) {
    String contig = bedFeature.getContig().replaceAll("chr", "");
    return StringUtils.isNumeric(contig) && Integer.parseInt(contig) < 23;
  }

  /**
   * @param file bed file
   * @param requireIndex , if the .tbi index is required for query
   * @return {@link BEDFileReader}
   */
  static BEDFileReader getReader(String file, boolean requireIndex) {
    return new BEDFileReader(file, requireIndex);
  }

  static class BEDFileReader implements Iterable<BEDFeature> {

    private final FeatureReader<BEDFeature> reader;

    /**
     * TODO : currently only handles .tbi indices
     */
    private BEDFileReader(final String file, final boolean requireIndex) {
      reader = AbstractFeatureReader.getFeatureReader(file, new BEDCodec(), requireIndex);
    }

    // /** Queries for records within the region specified. */
    // public CloseableIterator<BEDFeature> query(final String chrom, final int
    // start, final int end) {
    // try {
    // return reader.query(chrom, start, end);
    // } catch (final IOException ioe) {
    // throw new TribbleException("Could not create an iterator from a feature
    // reader.", ioe);
    // }
    // }

    public void close() {
      try {
        reader.close();
      } catch (final IOException ioe) {
        throw new TribbleException("Could not close a bed context feature reader.", ioe);
      }
    }

    /** Returns an iterator over all records in this bed file */
    public CloseableIterator<BEDFeature> iterator() {
      try {
        return reader.iterator();
      } catch (final IOException ioe) {
        throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
      }
    }

  }

  /**
   * Helper to test for overlapping regions sourced from a bed file
   */
  static class BEDOverlapDetector {

    private OverlapDetector<BEDFeature> detector;
    private int numOverlapped;

    /**
     * @param bedFile
     */
    BEDOverlapDetector(String bedFile, Logger log) {
      super();
      numOverlapped = 0;
      if (bedFile != null) {
        if (!FileOps.fileExists(bedFile)) {
          String err = "BED file to exclude did not exist: " + bedFile;
          log.severe(err);
          throw new IllegalArgumentException(err);
        }
        List<BEDFeature> regions = loadAll(bedFile);
        log.info("Loaded " + regions.size() + " regions to detect overlaps from " + bedFile);
        detector = OverlapDetector.create(regions);
      } else {
        detector = OverlapDetector.create(new ArrayList<>());
      }
    }

    boolean overlapsAny(Locatable query) {
      boolean olap = detector.overlapsAny(query);
      if (olap) {
        numOverlapped++;
      }
      return olap;
    }

    boolean overlapsNone(Locatable query) {
      return !overlapsAny(query);
    }

    /**
     * @return the numExcluded
     */
    int getNumExcluded() {
      return numOverlapped;
    }

  }
}
