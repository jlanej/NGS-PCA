package org.pankratzlab.ngspca;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.pankratzlab.ngspca.BedUtils.BEDOverlapDetector;
import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;
import org.pankratzlab.ngspca.RandomizedSVD.DISTRIBUTION;

/**
 * A simplified version of BamImport that uses MosDepth output or a custom input matrix to generate
 * PCS
 */
public class NGSPCA {

  private static void runInputMatrix(String inputMatrixFile, String outputDir, int numPcs,
                                     int niters, int numOversamples, int randomSeed,
                                     boolean overwrite, boolean normMatrix, DISTRIBUTION d,
                                     Logger log) throws InterruptedException, ExecutionException,
                                                 IOException {
    new File(outputDir).mkdirs();

    String delim = "\t";
    boolean gz = inputMatrixFile.endsWith(".gz");
    if (gz) {
      log.info("Assuming " + inputMatrixFile + " is gzipped");
    }

    log.info("Determining number of samples in " + inputMatrixFile);

    // the column names of a supplied matrix are used as given: whatever produced it named the
    // samples, and is in a better position than NGS-PCA to have named them correctly
    List<String> samples = FileOps.getFileHeader(inputMatrixFile, gz, delim, log);
    samples.remove(0);
    log.info("Found a total of " + samples.size() + " samples in " + inputMatrixFile);

    log.info("Determining number of regions in " + inputMatrixFile);
    List<String> regions = FileOps.getColumn(inputMatrixFile, gz, delim, 0, log);

    regions.remove(0);
    log.info("Found a total of " + regions.size() + " regions in " + inputMatrixFile);

    BlockRealMatrix dm;

    String tmpNormDm = Paths.get(outputDir, "tmp.mat.ser.gz").toString();
    if (!FileOps.fileExists(tmpNormDm) || overwrite) {
      log.info("Populating matrix from " + inputMatrixFile);
      log.info("Initializing matrix to " + samples.size() + " columns and " + regions.size()
               + " rows");
      dm = new BlockRealMatrix(regions.size(), samples.size());
      int[] rowIndex = {0};
      //add data to matrix, skipping header and first column of file
      try (Stream<String> stream = gz ? FileOps.gzLines(Paths.get(inputMatrixFile), log)
                                      : Files.lines(Paths.get(inputMatrixFile))) {
        stream.skip(1).map(l -> l.split(delim))
              .forEach(a -> dm.setRow(rowIndex[0]++,
                                      Utils.convertToDoubleArray(Arrays.copyOfRange(a, 1, a.length),
                                                                 0, log)));
      }
      if (normMatrix) {
        log.info("Normalizing input matrix");
        // the medians are not reported here: NGS-PCA did not choose the rows of a supplied matrix,
        // so it cannot say they are the autosomal, exclusion-filtered bins that make the median
        // mean what a reader of that file would take it to mean
        NormalizationOperations.foldChangeAndCenterRows(dm, log);
      }
      FileOps.writeSerial(dm, tmpNormDm, log);
    } else {
      dm = readCachedMatrix(tmpNormDm, log);
    }
    computeSVD(outputDir, numPcs, niters, numOversamples, randomSeed, log, samples, regions, dm, d);

  }

  /**
   * A serialized matrix that cannot be read back is what an interrupted run leaves behind: the file
   * exists, so the next run reuses it rather than rebuilding it, and the null it comes back as
   * would otherwise surface later as something unrelated to the real problem
   *
   * @param file the serialized matrix to read
   * @param log
   * @return the matrix it holds
   */
  private static BlockRealMatrix readCachedMatrix(String file, Logger log) {
    log.info("Loading existing serialized file " + file);
    BlockRealMatrix dm = (BlockRealMatrix) FileOps.readSerial(file, log);
    if (dm == null) {
      throw new IllegalStateException("Unable to read " + file
                                      + " - an interrupted run can leave it incomplete; re-run with --"
                                      + CmdLine.OVERWRITE_ARG + " to rebuild it");
    }
    return dm;
  }

  /**
   * @param input directory or file listing full paths containing MosDepth results, with extension
   *          {@link MosdepthUtils#MOSDEPTH_BED_EXT}
   * @param outputDir where results will be written
   * @param bedExclude if not null, regions overlapping this bed file will not be included
   * @param regionStrategy how to select markers for PCA
   * @param numPcs number of PCs to retain in the output file
   * @param sampleAt sample the mosdepth bins, once per this number
   * @param randomSeed random seed for sampling matrix
   * @param overwrite overwrite any existing output
   * @param threads number of threads for loading bed files
   * @param log
   * @throws InterruptedException
   * @throws ExecutionException
   * @throws IOException
   */
  private static void runMosdepth(String input, String outputDir, String bedExclude,
                                  REGION_STRATEGY regionStrategy, int numPcs, int niters,
                                  int numOversamples, int sampleAt, int randomSeed,
                                  boolean overwrite, int threads, String sampleSuffix,
                                  DISTRIBUTION d,
                                  Logger log) throws InterruptedException, ExecutionException,
                                              IOException {
    new File(outputDir).mkdirs();

    String[] extensions = new String[] {MosdepthUtils.MOSDEPTH_BED_EXT};

    // get all files with mosdepth bed extension
    List<String> mosDepthResultFiles;
    if (FileOps.isDir(input) && FileOps.dirExists(input)) {
      log.info("Detected " + input + " is a directory, searching for "
               + MosdepthUtils.MOSDEPTH_BED_EXT + " extensions");
      mosDepthResultFiles = FileOps.listFilesWithExtension(input, extensions);
    } else if (FileOps.fileExists(input)) {
      log.info("Detected " + input + " is a file, reading mosdepth result file paths");
      mosDepthResultFiles = FileOps.readFile(input);
    } else {
      String err = "Invalid or non-existent input argument " + input;
      log.severe(err);
      throw new IllegalArgumentException(err);
    }

    if (mosDepthResultFiles.isEmpty()) {

      String err = "No input files were found";
      log.severe(err);
      throw new IllegalArgumentException(err);
    } else {
      log.info("Detected " + mosDepthResultFiles.size() + " mosdepth input files in " + input);
    }
    // parse sample names from files
    List<String> samples = SampleNames.resolve(mosDepthResultFiles.stream()
                                                                  .map(f -> FileOps.stripDirectoryAndExtension(f,
                                                                                                               MosdepthUtils.MOSDEPTH_BED_EXT))
                                                                  .collect(Collectors.toList()),
                                               sampleSuffix, log);

    // load ucsc regions to use

    BEDOverlapDetector overlapDetector = new BEDOverlapDetector(bedExclude, log);
    List<String> regions = MosdepthUtils.getRegionsToUse(mosDepthResultFiles.get(0), regionStrategy,
                                                         overlapDetector, log);
    log.info(overlapDetector.getNumExcluded() + " regions removed during up-front filtering");
    if (sampleAt > 1) {
      log.info("Sampling the " + regions.size() + " mosdepth regions once every " + sampleAt
               + " bins");
      regions = IntStream.range(0, regions.size()).filter(n -> n % sampleAt == 0)
                         .mapToObj(regions::get).collect(Collectors.toList());
      log.info("Sampled " + regions.size() + " bins");

    }
    // Store the raw input matrix
    String tmpRawDm = Paths.get(outputDir, "tmp.raw.ser.gz").toString();
    // Store the temporary input matrix
    String tmpNormDm = Paths.get(outputDir, "tmp.mat.ser.gz").toString();
    // Report the medians normalization computed, for tools that normalize against them too
    String medianDm = Paths.get(outputDir, CoverageMedians.AUTOSOMAL_FILE).toString();

    // populate input matrix and normalize
    BlockRealMatrix dm;
    if (!FileOps.fileExists(tmpNormDm) || overwrite) {
      MosdepthUtils.NormalizedResult normalized = MosdepthUtils.processFiles(mosDepthResultFiles,
                                                                             new HashSet<String>(regions),
                                                                             tmpRawDm, threads,
                                                                             log);
      dm = normalized.matrix;
      CoverageMedians.write(medianDm, samples, normalized.columnMedians, regions.size(), log);
      FileOps.writeSerial(dm, tmpNormDm, log);
    } else {
      dm = readCachedMatrix(tmpNormDm, log);
      if (!FileOps.fileExists(medianDm)) {
        // the medians come from reading the mosdepth files, which reusing the matrix skips
        log.warning("Reused " + tmpNormDm + ", so " + medianDm + " was not written - re-run with --"
                    + CmdLine.OVERWRITE_ARG + " to produce it");
      }
    }
    //    String inputMatrix = Paths.get(outputDir, "svd.norm.input.txt").toString();
    //    log.info("Writing to " + inputMatrix);
    //
    //    RandomizedSVD.dumpMatrix(inputMatrix, dm, "BIN", samples.toArray(new String[samples.size()]),
    //                             regions.toArray(new String[regions.size()]), false, log);

    computeSVD(outputDir, numPcs, niters, numOversamples, randomSeed, log, samples, regions, dm, d);
  }

  static void computeSVD(String outputDir, int numPcs, int niters, int numOversamples,
                         int randomSeed, Logger log, List<String> samples, List<String> regions,
                         BlockRealMatrix dm, DISTRIBUTION d) {
    RandomizedSVD svd = new RandomizedSVD(samples, regions, log);

    log.info("Oversampling set to: " + numOversamples);
    log.info("Subspace iterations set to: " + niters);
    log.info("Random seed set to: " + randomSeed);
    svd.fit(dm, numPcs, niters, numOversamples, randomSeed, d);
    // perform SVD

    String pcs = Paths.get(outputDir, "svd.pcs.txt").toString();
    String loadings = Paths.get(outputDir, "svd.loadings.txt").toString();
    String singularValues = Paths.get(outputDir, "svd.singularvalues.txt").toString();
    String binsUsed = Paths.get(outputDir, "svd.bins.txt").toString();
    String samplesUsed = Paths.get(outputDir, "svd.samples.txt").toString();

    log.info("Writing to " + pcs);
    svd.dumpPCsToText(pcs, log);
    log.info("Writing to " + loadings);
    svd.computeAndDumpLoadings(loadings, log);
    log.info("Writing to " + singularValues);
    svd.dumpSingularValuesToText(singularValues, log);
    log.info("Writing to " + binsUsed);
    FileOps.writeToText(svd.getRowNames(), binsUsed, log);
    log.info("Writing to " + samplesUsed);
    FileOps.writeToText(svd.getColumnNames(), samplesUsed, log);
  }

  public static void main(String[] args) {
    Logger log = Logger.getLogger(NGSPCA.class.getName());
    Options options = CmdLine.generateOptions();
    if (CmdLine.isHelpRequested(args)) {
      // asking for the usage message is not a failure
      CmdLine.printHelp(log, options);
      return;
    }
    CommandLine cmd = CmdLine.generateCommandLine(log, options, args);
    if (cmd == null) {
      CmdLine.printHelp(log, options);
      System.exit(1);
    }

    String input = cmd.getOptionValue(CmdLine.INPUT_ARG);
    String outputDir = cmd.getOptionValue(CmdLine.OUTPUT_DIR_ARG);

    try {
      int numPcs = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_COMPONENTS_ARG,
                                                       Integer.toString(CmdLine.DEFAULT_PCS)));

      int threads = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_THREADS_ARG,
                                                        Integer.toString(CmdLine.DEFAULT_THREADS)));

      int sampleAt = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_SAMPLE_ARG,
                                                         Integer.toString(CmdLine.DEFAULT_SAMPLE)));

      int niters = Integer.parseInt(cmd.getOptionValue(CmdLine.N_ITERS,
                                                       Integer.toString(RandomizedSVD.DEFAULT_NITERS)));
      int numOversamples = Integer.parseInt(cmd.getOptionValue(CmdLine.OVERSAMPLE,
                                                               Integer.toString(RandomizedSVD.DEFAULT_OVERSAMPLES)));

      int randomSeed = Integer.parseInt(cmd.getOptionValue(CmdLine.RANDOM_SEED,
                                                           Integer.toString(CmdLine.DEFAULT_RANDOM_SEED)));

      DISTRIBUTION d = DISTRIBUTION.valueOf(cmd.getOptionValue(CmdLine.DISTRIBUTION_ARG,
                                                               CmdLine.DEFAULT_DISTRIBUTION.toString()));
      String bedExclude = cmd.getOptionValue(CmdLine.EXCLUDE_BED_FILE,
                                             CmdLine.DEFAULT_EXCLUDE_BED_FILE);
      String sampleSuffix = cmd.getOptionValue(CmdLine.SAMPLE_SUFFIX_ARG,
                                               CmdLine.DEFAULT_SAMPLE_SUFFIX);
      if (cmd.hasOption(CmdLine.MATRIX_INPUT_ARG)) {
        if (sampleSuffix != null) {
          // it does nothing here, and a flag that quietly does nothing is what the same flag
          // refuses to be when it matches no mosdepth file name
          log.warning("--" + CmdLine.SAMPLE_SUFFIX_ARG + " does not apply to --"
                      + CmdLine.MATRIX_INPUT_ARG
                      + " input, whose column names are used as they are");
        }
        runInputMatrix(input, outputDir, numPcs, niters, numOversamples, randomSeed,
                       cmd.hasOption(CmdLine.OVERWRITE_ARG),
                       cmd.hasOption(CmdLine.NORM_MATRIX_INPUT_ARG), d, log);
      } else {
        runMosdepth(input, outputDir, bedExclude, REGION_STRATEGY.AUTOSOMAL, numPcs, niters,
                    numOversamples, sampleAt, randomSeed, cmd.hasOption(CmdLine.OVERWRITE_ARG),
                    threads, sampleSuffix, d, log);
      }
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      log.severe("An exception occurred while running\nFeel free to open an issue at https://github.com/PankratzLab/NGS-PCA after reviewing the help message below");
      CmdLine.printHelp(log, options);
      // a run that produced no results must not report success - a workflow engine or job array
      // has nothing else to go on
      System.exit(1);
    }
  }
}
