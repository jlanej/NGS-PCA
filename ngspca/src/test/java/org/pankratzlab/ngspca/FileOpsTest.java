package org.pankratzlab.ngspca;

import java.io.File;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * Tests that a write which did not happen is not mistaken for one that did
 */
public class FileOpsTest extends TestCase {

  private static final Logger LOG = Logger.getLogger(FileOpsTest.class.getName());

  /**
   * @return a path that cannot be written to, because a directory is already there
   */
  private static String unwritablePath() throws Exception {
    File dir = Files.createTempDirectory("ngspca.unwritable").toFile();
    dir.deleteOnExit();
    return dir.getAbsolutePath();
  }

  public void testWriteToTextFailsLoudly() throws Exception {
    try {
      FileOps.writeToText(Arrays.asList("a", "b"), unwritablePath(), LOG);
      fail("expected an UncheckedIOException when the file cannot be written");
    } catch (UncheckedIOException expected) {
      // expected
    }
  }

  /**
   * These are read back by a later run, so a failed write that is allowed to pass leaves the next
   * run to discover it
   */
  public void testWriteSerialFailsLoudly() throws Exception {
    try {
      FileOps.writeSerial(Arrays.asList("a", "b"), unwritablePath(), LOG);
      fail("expected an UncheckedIOException when the file cannot be written");
    } catch (UncheckedIOException expected) {
      // expected
    }
  }

  public void testWriteToTextRoundTrips() throws Exception {
    File out = File.createTempFile("ngspca.text", ".txt");
    out.deleteOnExit();

    FileOps.writeToText(Arrays.asList("a", "b"), out.getAbsolutePath(), LOG);

    assertEquals(Arrays.asList("a", "b"), Files.readAllLines(out.toPath()));
  }
}
