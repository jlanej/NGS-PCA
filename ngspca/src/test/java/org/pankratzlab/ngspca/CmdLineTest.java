package org.pankratzlab.ngspca;

import junit.framework.TestCase;

/**
 * Tests command line handling that has to happen before the arguments are parsed
 */
public class CmdLineTest extends TestCase {

  public void testRecognizesEveryFormOfTheHelpFlag() {
    assertTrue(CmdLine.isHelpRequested(new String[] {"-h"}));
    assertTrue(CmdLine.isHelpRequested(new String[] {"-help"}));
    assertTrue(CmdLine.isHelpRequested(new String[] {"--help"}));
  }

  /**
   * The required options are missing from a command line that only asks for help, so the request
   * has to be found wherever it appears rather than read off a parsed result
   */
  public void testFindsHelpAmongOtherArguments() {
    assertTrue(CmdLine.isHelpRequested(new String[] {"-input", "in.txt", "--help"}));
    assertTrue(CmdLine.isHelpRequested(new String[] {"--help", "-input", "in.txt"}));
  }

  public void testDoesNotMistakeOtherArgumentsForHelp() {
    assertFalse(CmdLine.isHelpRequested(new String[] {}));
    assertFalse(CmdLine.isHelpRequested(new String[] {"-input", "in.txt", "-outputDir", "out/"}));
    // a path that merely contains the word, and the flag glued to a value
    assertFalse(CmdLine.isHelpRequested(new String[] {"-input", "/data/help/in.txt"}));
    assertFalse(CmdLine.isHelpRequested(new String[] {"-helpful"}));
  }
}
