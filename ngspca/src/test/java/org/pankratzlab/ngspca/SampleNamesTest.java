package org.pankratzlab.ngspca;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * Tests the sample identifiers every output names samples by
 */
public class SampleNamesTest extends TestCase {

  private static final Logger LOG = Logger.getLogger(SampleNamesTest.class.getName());

  public void testNoSuffixLeavesNamesAlone() {
    List<String> names = Arrays.asList("sampleA.cram.by1000.", "sampleB.cram.by1000.");

    assertEquals(names, SampleNames.resolve(names, null, LOG));
    assertEquals(names, SampleNames.resolve(names, "", LOG));
  }

  /**
   * The stripped extension carries no leading dot, so the names it leaves behind end in one - the
   * suffix has to be able to take that with it
   */
  public void testRemovesALiteralTrailingSuffix() {
    List<String> resolved = SampleNames.resolve(Arrays.asList("sampleA.cram.by1000.",
                                                              "sampleB.cram.by1000."),
                                                ".by1000.", LOG);

    assertEquals(Arrays.asList("sampleA.cram", "sampleB.cram"), resolved);
  }

  /**
   * A literal suffix, not a pattern, and only where it actually trails
   */
  public void testOnlyRemovesTheSuffixWhereItTrails() {
    List<String> resolved = SampleNames.resolve(Arrays.asList("by1000.sampleA.cram.by1000.",
                                                              "sampleB.by1000.cram"),
                                                ".by1000.", LOG);

    assertEquals("by1000.sampleA.cram", resolved.get(0));
    assertEquals("sampleB.by1000.cram", resolved.get(1));
  }

  /**
   * A suffix that matches nothing looks downstream exactly like not passing one, so it is a typo
   * rather than a no-op
   */
  public void testRefusesASuffixThatMatchesNothing() {
    try {
      SampleNames.resolve(Arrays.asList("sampleA.cram.", "sampleB.cram."), ".by500.", LOG);
      fail("expected an IllegalArgumentException for a suffix no sample name ends with");
    } catch (IllegalArgumentException expected) {
      assertTrue("message should quote the suffix: " + expected.getMessage(),
                 expected.getMessage().contains(".by500."));
    }
  }

  /**
   * Two samples sharing a name cannot be told apart in any output, so the collision has to stop the
   * run rather than reach a reader
   */
  public void testRefusesNamesMadeDuplicateBySuffixRemoval() {
    try {
      SampleNames.resolve(Arrays.asList("sample.cram.by1000.", "sample.cram"), ".by1000.", LOG);
      fail("expected an IllegalArgumentException for a collision created by suffix removal");
    } catch (IllegalArgumentException expected) {
      assertTrue("message should name the duplicate: " + expected.getMessage(),
                 expected.getMessage().contains("sample.cram"));
    }
  }

  public void testRefusesDuplicatesWithoutASuffix() {
    try {
      SampleNames.resolve(Arrays.asList("sample.cram.", "sample.cram."), null, LOG);
      fail("expected an IllegalArgumentException for duplicate input names");
    } catch (IllegalArgumentException expected) {
      assertTrue("message should name the duplicate: " + expected.getMessage(),
                 expected.getMessage().contains("sample.cram."));
    }
  }
}
