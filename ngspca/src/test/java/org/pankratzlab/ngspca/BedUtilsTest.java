package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
import junit.framework.TestCase;
import org.pankratzlab.ngspca.BedUtils.BedRegionResult;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.bed.BEDFeature;

/**
 * {@link BedUtils#loadSpecificRegions} reads mosdepth beds directly rather than through htsjdk's
 * {@link htsjdk.tribble.bed.BEDCodec}, and these hold it to what the codec produced - the same
 * regions matched against the same 1-based keys, and the same coverage values to the bit.
 */
public class BedUtilsTest extends TestCase {

  /**
   * Bins across three contigs with decimal coverages, the shape mosdepth emits
   */
  private static List<String> bedLines() {
    List<String> lines = new ArrayList<>();
    int start = 0;
    for (int i = 0; i < 120; i++) {
      String contig = i % 3 == 2 ? "chrX" : "chr" + (1 + i % 2);
      lines.add(contig + "\t" + start + "\t" + (start + 1000) + "\t" + (i + i * 0.01));
      start += 1000;
    }
    return lines;
  }

  private static File write(List<String> lines, boolean gzip) throws Exception {
    File file = File.createTempFile("ngspca.bedutils", ".regions.bed.gz");
    file.deleteOnExit();
    try (OutputStream out = gzip ? new GZIPOutputStream(new FileOutputStream(file))
                                 : new FileOutputStream(file)) {
      out.write(String.join("\n", lines).concat("\n").getBytes(StandardCharsets.US_ASCII));
    }
    return file;
  }

  /**
   * Every other region of the file, keyed the way the region list is built in production - by
   * htsjdk, so the keys the parser must match are the codec's own
   */
  private static Set<String> everyOtherRegion(File file) {
    Set<String> regions = new LinkedHashSet<>();
    List<BEDFeature> features = BedUtils.loadAll(file.getAbsolutePath());
    for (int i = 0; i < features.size(); i += 2) {
      regions.add(BedUtils.getBedUCSC(features.get(i)));
    }
    return regions;
  }

  /**
   * The loader this replaced: htsjdk features filtered by key, coverage parsed from the name
   */
  private static double[] htsjdkReference(File file, Set<String> regions) {
    List<Double> values = new ArrayList<>();
    BedUtils.BEDFileReader reader = BedUtils.getReader(file.getAbsolutePath(), false);
    for (BEDFeature feature : reader) {
      if (regions.contains(BedUtils.getBedUCSC(feature))) {
        values.add(Double.parseDouble(feature.getName()));
      }
    }
    reader.close();
    double[] result = new double[values.size()];
    for (int i = 0; i < result.length; i++) {
      result[i] = values.get(i);
    }
    return result;
  }

  private static void assertMatchesReference(File file, Set<String> regions) {
    double[] expected = htsjdkReference(file, regions);
    BedRegionResult actual = BedUtils.loadSpecificRegions(file.getAbsolutePath(), regions);

    assertEquals(file.getAbsolutePath(), actual.file());
    assertEquals(expected.length, actual.matched());
    for (int i = 0; i < expected.length; i++) {
      assertEquals("row " + i, Double.doubleToRawLongBits(expected[i]),
                   Double.doubleToRawLongBits(actual.coverage()[i]));
    }
  }

  public void testMatchesTheHtsjdkPathItReplaced() throws Exception {
    File file = write(bedLines(), true);
    assertMatchesReference(file, everyOtherRegion(file));
  }

  /**
   * Real mosdepth output is bgzf - concatenated gzip members - which a plain GZIPInputStream must
   * also read to the end of
   */
  public void testReadsBgzf() throws Exception {
    File file = File.createTempFile("ngspca.bedutils.bgzf", ".regions.bed.gz");
    file.deleteOnExit();
    try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(file)) {
      out.write(String.join("\n", bedLines()).concat("\n").getBytes(StandardCharsets.US_ASCII));
    }
    assertMatchesReference(file, everyOtherRegion(file));
  }

  /**
   * htsjdk trusts the .gz extension and fails on uncompressed content; the replacement has to
   * fail such a file too, not quietly read it where the old loader stopped the run
   */
  public void testRefusesPlainTextNamedGz() throws Exception {
    File file = write(bedLines(), false);
    File gzipped = write(bedLines(), true);
    Set<String> regions = everyOtherRegion(gzipped);
    try {
      BedUtils.loadSpecificRegions(file.getAbsolutePath(), regions);
      fail("expected a plain-text file named .gz to be refused");
    } catch (RuntimeException expected) {
      // expected - the old htsjdk path threw for the same input
    }
  }

  /**
   * Lines the codec skipped stay skipped
   */
  public void testSkipsHeaderAndBlankLines() throws Exception {
    List<String> lines = new ArrayList<>();
    lines.add("# a comment");
    lines.add("track name=ignored");
    lines.add("browser position chr1");
    lines.add("");
    lines.addAll(bedLines());
    File file = write(lines, true);
    assertMatchesReference(file, everyOtherRegion(file));
  }

  /**
   * The consumer's count check needs duplicates counted and misses missing, not clamped away
   */
  public void testCountsDuplicatesAndMissesForTheConsumerCheck() throws Exception {
    List<String> lines = bedLines();
    File file = write(lines, true);
    Set<String> regions = everyOtherRegion(file);

    List<String> duplicated = new ArrayList<>(lines);
    duplicated.add(lines.get(0));
    BedRegionResult duplicate = BedUtils.loadSpecificRegions(write(duplicated, true).getAbsolutePath(),
                                                             regions);
    assertEquals(regions.size() + 1, duplicate.matched());

    BedRegionResult missing = BedUtils.loadSpecificRegions(write(lines.subList(2, lines.size()),
                                                                 true).getAbsolutePath(),
                                                           regions);
    assertEquals(regions.size() - 1, missing.matched());
  }
}
