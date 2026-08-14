package org.pankratzlab.ngspca;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.pankratzlab.ngspca.RandomizedSVD.DISTRIBUTION;

/**
 * Utility class to process cmdline arguments
 */
class CmdLine {

  static final String INPUT_ARG = "input";
  static final String MATRIX_INPUT_ARG = "matrix";
  static final String NORM_MATRIX_INPUT_ARG = "normalizeMatrix";
  static final String OUTPUT_DIR_ARG = "outputDir";
  static final String OVERWRITE_ARG = "overwrite";
  static final String NUM_COMPONENTS_ARG = "numPC";
  static final String NUM_THREADS_ARG = "threads";
  static final String NUM_SAMPLE_ARG = "sampleEvery";
  static final String SAMPLE_SUFFIX_ARG = "sampleSuffix";
  static final String EXCLUDE_BED_FILE = "bedExclude";
  static final String N_ITERS = "iters";
  static final String OVERSAMPLE = "oversample";
  static final String RANDOM_SEED = "randomSeed";
  static final String DISTRIBUTION_ARG = "distribution";

  static final int DEFAULT_RANDOM_SEED = 42;
  static final int DEFAULT_PCS = 20;
  static final int DEFAULT_SAMPLE = 1;
  static final String DEFAULT_EXCLUDE_BED_FILE = null;
  static final String DEFAULT_SAMPLE_SUFFIX = null;
  static final DISTRIBUTION DEFAULT_DISTRIBUTION = DISTRIBUTION.UNIFORM;

  static final int DEFAULT_THREADS = 24;

  static final String HELP = "help";

  private CmdLine() {
    // utility
  }

  /**
   * Get available cmd line options
   * 
   * @return the {@link Options} available
   */
  static Options generateOptions() {
    final Option help = Option.builder("h").required().longOpt(HELP).hasArg()
                              .desc("Print application usage and exit").hasArg(false)
                              .required(false).build();
    final Option overwrite = Option.builder("r").required().longOpt(OVERWRITE_ARG).hasArg()
                                   .desc("Overwrite existing temporary files and recompute each step")
                                   .hasArg(false).required(false).build();
    final Option inputOption = Option.builder("i").hasArg(true).longOpt(INPUT_ARG)
                                     .desc("""
                                           An existing directory containing mosdepth result files \
                                           (*%s extension) OR a file listing paths to mosdepth \
                                           result files, one result file per line OR a matrix of \
                                           values to directly use (see the %s argument)"""
                                           .formatted(MosdepthUtils.MOSDEPTH_BED_EXT,
                                                      MATRIX_INPUT_ARG))
                                     .required(true).build();

    final Option outputOption = Option.builder("o").hasArg(true).required().longOpt(OUTPUT_DIR_ARG)
                                      .hasArg()
                                      .desc("PCA results and auxiliary files will be placed here")
                                      .required().build();
    final Option numComponents = Option.builder("n").hasArg(true).required()
                                       .longOpt(NUM_COMPONENTS_ARG).hasArg()
                                       .desc("""
                                             The number of PCs to retain - the minimum of this and \
                                             the number of markers/samples will be retained. \
                                             Default is %d""".formatted(DEFAULT_PCS))
                                       .required(false).build();

    final Option numThreads = Option.builder("t").hasArg(true).required().longOpt(NUM_THREADS_ARG)
                                    .hasArg()
                                    .desc("""
                                          Number of threads to utilize, both when loading data and \
                                          in the parallel steps of the decomposition. Set this to \
                                          the cores your job was allocated - the default suits a \
                                          node, not a workstation, and asking for more than are \
                                          available is reported. Default is \
                                          %d""".formatted(DEFAULT_THREADS))
                                    .required(false).build();

    final Option sampleEvery = Option.builder("s").hasArg(true).required().longOpt(NUM_SAMPLE_ARG)
                                     .hasArg()
                                     .desc("""
                                           Sample mosdepth bins. For example, --%s 10 would select \
                                           every tenth bin. Default is %d (use every bin)"""
                                           .formatted(NUM_SAMPLE_ARG, DEFAULT_SAMPLE))
                                     .required(false).build();

    final Option sampleSuffix = Option.builder(SAMPLE_SUFFIX_ARG).hasArg(true)
                                      .longOpt(SAMPLE_SUFFIX_ARG).hasArg()
                                      .desc("""
                                            Optional: remove this literal suffix from every sample \
                                            name. Sample names are otherwise the mosdepth file name \
                                            with the %s extension removed, which leaves a trailing \
                                            '.' - so --%s .by1000. names \
                                            sample.cram.by1000.regions.bed.gz as sample.cram. \
                                            Applies to every output file that names samples, so \
                                            they stay consistent with each other. Mosdepth input \
                                            only - the column names of a matrix given with --%s are \
                                            used as they are"""
                                            .formatted(MosdepthUtils.MOSDEPTH_BED_EXT,
                                                       SAMPLE_SUFFIX_ARG, MATRIX_INPUT_ARG))
                                      .required(false).build();

    final Option bedExcludes = Option.builder("b").hasArg(true).required().longOpt(EXCLUDE_BED_FILE)
                                     .hasArg()
                                     .desc("""
                                           Optional: Provide a file to exclude specific regions \
                                           from PCA input, prior to sampling with %s"""
                                           .formatted(NUM_SAMPLE_ARG))
                                     .required(false).build();
    final Option niter = Option.builder(N_ITERS).hasArg(true).required().longOpt(N_ITERS).hasArg()
                               .desc("""
                                     Specifies the number of power (subspace) iterations to reduce \
                                     the approximation error. The power scheme is recommended, if \
                                     the singular values decay slowly. In practice, 2 or 3 \
                                     iterations achieve good results, however, computing power \
                                     iterations increases the computational costs. Default is \
                                     %d""".formatted(RandomizedSVD.DEFAULT_NITERS))
                               .required(false).build();
    final Option oversamples = Option.builder(OVERSAMPLE).hasArg(true).required()
                                     .longOpt(OVERSAMPLE).hasArg()
                                     .desc("""
                                           An oversampling parameter to improve the approximation \
                                           of the randomized PCA. A value of at least 10 is \
                                           recommended. Default is \
                                           %d""".formatted(RandomizedSVD.DEFAULT_OVERSAMPLES))
                                     .required(false).build();

    final Option randomSeed = Option.builder(RANDOM_SEED).hasArg(true).required()
                                    .longOpt(RANDOM_SEED).hasArg()
                                    .desc("""
                                          Random seed for generating sampling matrix for randomized \
                                          PCA (probably not worth changing the default). Default is \
                                          %d""".formatted(DEFAULT_RANDOM_SEED))
                                    .required(false).build();

    final Option distribution = Option.builder(DISTRIBUTION_ARG).hasArg(true)
                                      .longOpt(DISTRIBUTION_ARG).hasArg()
                                      .desc("""
                                            The distribution used to seed the initial random \
                                            matrix. Options are %s or %s. Default is %s"""
                                            .formatted(DISTRIBUTION.UNIFORM, DISTRIBUTION.GAUSSIAN,
                                                       DEFAULT_DISTRIBUTION))
                                      .required(false).build();
    final Option matrix = Option.builder(MATRIX_INPUT_ARG).hasArg(false).longOpt(MATRIX_INPUT_ARG)
                                .desc("""
                                      The input provided by %s is a matrix (i.e. SVD will be \
                                      performed directly on the matrix, without normalization, to \
                                      generate PCS""".formatted(INPUT_ARG))
                                .required(false).build();
    final Option normMatrix = Option.builder(NORM_MATRIX_INPUT_ARG).hasArg(false)
                                    .longOpt(NORM_MATRIX_INPUT_ARG)
                                    .desc("""
                                          The input provided by %s is a matrix and should be \
                                          normalized (log2 by sample/column, centered by \
                                          row)""".formatted(INPUT_ARG))
                                    .required(false).build();
    final Options options = new Options();
    options.addOption(help);

    options.addOption(inputOption);
    options.addOption(matrix);
    options.addOption(normMatrix);
    options.addOption(outputOption);
    options.addOption(numComponents);
    options.addOption(numThreads);
    options.addOption(sampleEvery);
    options.addOption(sampleSuffix);
    options.addOption(bedExcludes);
    options.addOption(niter);
    options.addOption(oversamples);
    options.addOption(randomSeed);
    options.addOption(distribution);
    options.addOption(overwrite);

    return options;
  }

  /**
   * {@link #INPUT_ARG} and {@link #OUTPUT_DIR_ARG} are required, so a command line that asks only
   * for the usage message does not parse - the request has to be spotted before parsing rather than
   * read off the parsed result
   *
   * @param commandLineArguments arguments provided to the cmd line
   * @return whether the usage message was asked for
   */
  static boolean isHelpRequested(final String[] commandLineArguments) {
    for (String argument : commandLineArguments) {
      if ("-h".equals(argument) || ("-" + HELP).equals(argument)
          || ("--" + HELP).equals(argument)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Parse cmd line options
   *
   * @param log {@link Logger}
   * @param options {@link Options} to select from
   * @param commandLineArguments arguments provided to the cmd line
   * @return parsed {@link CommandLine}
   */
  static CommandLine generateCommandLine(Logger log, final Options options,
                                         final String[] commandLineArguments) {
    final CommandLineParser cmdLineParser = new DefaultParser();
    CommandLine commandLine = null;
    try {
      commandLine = cmdLineParser.parse(options, commandLineArguments);
    } catch (ParseException parseException) {
      log.severe(parseException.getMessage());
    }
    return commandLine;
  }

  /**
   * Print the help usage to the log, based on {@link Options} available
   * 
   * @param log {@link Logger}
   * @param options {@link Options}
   */
  static void printHelp(Logger log, final Options options) {
    StringWriter out = new StringWriter();
    PrintWriter pw = new PrintWriter(out);

    HelpFormatter formatter = new HelpFormatter();
    formatter.printHelp(pw, 80, "ngspca", "USAGE:", options, 0, 0, "", true);
    pw.flush();
    pw.close();
    log.info(out.toString());
  }
}
