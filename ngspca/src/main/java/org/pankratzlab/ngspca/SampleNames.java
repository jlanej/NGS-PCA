package org.pankratzlab.ngspca;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Resolves the sample identifiers that every output file names samples by.
 * <p>
 * Sample names are the mosdepth file name with the {@link MosdepthUtils#MOSDEPTH_BED_EXT} extension
 * removed, which leaves whatever the file naming convention put in front of it - including, since
 * the extension carries no leading dot, a trailing one. Tools reading these outputs alongside a
 * depth matrix they built themselves have to match sample identifiers exactly, so
 * {@link CmdLine#SAMPLE_SUFFIX_ARG} removes a fixed suffix here rather than leaving every consumer
 * to rewrite the tables afterwards.
 */
class SampleNames {

  /**
   * how many duplicates to name before giving up on the list
   */
  private static final int MAX_REPORTED = 5;

  private SampleNames() {

  }

  /**
   * @param names sample identifiers, in matrix column order
   * @param suffix literal suffix to remove from each name, or null to use the names as they are
   * @param log
   * @return the identifiers to report, in the order given
   */
  static List<String> resolve(List<String> names, String suffix, Logger log) {
    List<String> resolved = names;
    if (suffix != null && !suffix.isEmpty()) {
      resolved = new ArrayList<>(names.size());
      int stripped = 0;
      for (String name : names) {
        if (name.endsWith(suffix)) {
          resolved.add(name.substring(0, name.length() - suffix.length()));
          stripped++;
        } else {
          resolved.add(name);
        }
      }
      if (stripped == 0) {
        // the flag was asked for and did nothing, which downstream looks the same as not passing it
        throw new IllegalArgumentException("No sample name ends with \"" + suffix + "\", so --"
                                           + CmdLine.SAMPLE_SUFFIX_ARG
                                           + " had no effect - sample names look like \""
                                           + names.get(0) + "\"");
      }
      log.info("Removed the suffix \"" + suffix + "\" from " + stripped + " of " + names.size()
               + " sample names");
      if (stripped < names.size()) {
        log.warning((names.size() - stripped) + " sample name(s) do not end with \"" + suffix
                    + "\" and were left as they are");
      }
    }
    requireUnique(resolved, suffix);
    return resolved;
  }

  /**
   * Every output addresses samples by name, so two samples sharing one is not recoverable later
   */
  private static void requireUnique(List<String> names, String suffix) {
    Set<String> seen = new HashSet<>();
    Set<String> duplicates = new LinkedHashSet<>();
    for (String name : names) {
      if (!seen.add(name)) {
        duplicates.add(name);
      }
    }
    if (duplicates.isEmpty()) {
      return;
    }
    StringBuilder message = new StringBuilder("Duplicate sample name(s) detected: ");
    int reported = 0;
    for (String duplicate : duplicates) {
      if (reported == MAX_REPORTED) {
        message.append(" and ").append(duplicates.size() - MAX_REPORTED).append(" more");
        break;
      }
      message.append(reported > 0 ? ", " : "").append(duplicate);
      reported++;
    }
    if (suffix != null && !suffix.isEmpty()) {
      message.append(" - removing the suffix \"").append(suffix)
             .append("\" made these names collide");
    }
    throw new IllegalArgumentException(message.toString());
  }
}
