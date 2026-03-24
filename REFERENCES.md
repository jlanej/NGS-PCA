# References

This document lists all resources, software, methods, and datasets used by or referenced in this
repository, along with a description of how each is used.

---

## 1. Core Algorithms

### 1.1 Randomized Singular Value Decomposition

**Halko, N., Martinsson, P.-G., & Tropp, J. A. (2011).** Finding structure with randomness:
Probabilistic algorithms for constructing approximate matrix decompositions. *SIAM Review*, 53(2),
217–288.
<https://doi.org/10.1137/090771806>
Preprint: <https://arxiv.org/abs/0909.4061>

> **How used:** The core algorithm implemented in `ngspca/src/main/java/org/pankratzlab/ngspca/RandomizedSVD.java`.
> NGS-PCA computes an approximate truncated SVD of the coverage matrix using the randomized
> subspace method described in this paper. The random test matrix and the power (subspace) iteration
> scheme follow this reference. The `-distribution GAUSSIAN` option corresponds exactly to the
> Gaussian random matrix specified in this paper; `-distribution UNIFORM` (the default) is an
> alternative that produces equivalent results in practice when power iterations are used.

---

### 1.2 Power Iteration Scheme for Subspace Approximation

**Rokhlin, V., Szlam, A., & Tygert, M. (2010).** A randomized algorithm for principal component
analysis. *SIAM Journal on Matrix Analysis and Applications*, 31(3), 1100–1124.
<https://doi.org/10.1137/080736417>

> **How used:** The subspace (power) iteration loop inside `RandomizedSVD.fit()` follows the
> scheme described in this paper. Each power iteration refines the column space estimate, improving
> approximation accuracy for matrices whose singular values decay slowly. Controlled by the
> `-iters` parameter (default: 10).

---

### 1.3 Randomized Matrix Decompositions — R Reference Implementation

**Erichson, N. B., Voronin, S., Brunton, S. L., & Kutz, J. N. (2019).** Randomized matrix
decompositions using R. *Journal of Statistical Software*, 89(11), 1–48.
<https://doi.org/10.18637/jss.v089.i11>
Preprint: <https://arxiv.org/abs/1608.02148>
GitHub: <https://github.com/erichson/rSVD>

> **How used:** The NGS-PCA Java implementation of randomized SVD is described in source
> comments as analogous to the rSVD R package from this reference. Cited in
> `RandomizedSVD.java` as a companion implementation.

---

## 2. Reference Datasets and Data Portals

### 2.1 1000 Genomes Project 30× High-Coverage WGS (NYGC)

**Byrska-Bishop, M., Evani, U. S., Zhao, X., Basile, A. O., Abel, H. J., Regier, A. A., …
Zody, M. C. (2022).** High-coverage whole-genome sequencing of the expanded 1000 Genomes
Project cohort including 602 trios. *Cell*, 185(18), 3426–3440.e19.
<https://doi.org/10.1016/j.cell.2022.08.004>

Data portal: <https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>
ENA projects: [PRJEB31736](https://www.ebi.ac.uk/ena/browser/view/PRJEB31736) (2,504 unrelated),
[PRJEB36890](https://www.ebi.ac.uk/ena/browser/view/PRJEB36890) (698 related),
[PRJEB55077](https://www.ebi.ac.uk/ena/browser/view/PRJEB55077) (3,202 combined)

> **How used:** The `example/1000G_highcov/` directory contains a fully reproducible end-to-end
> SLURM pipeline for downloading all 3,202 CRAM files via Aspera, computing coverage with
> mosdepth, and running NGS-PCA to produce ~200 PCs. The 18-sample chr1 subset is used for
> integration testing (see `example/` and `.github/workflows/integration-test.yml`).
> Please cite this paper when publishing results derived from this pipeline.

Sequence index files (used by `example/1000G_highcov/00_setup.sh`):
- 2,504 unrelated: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index`
- 698 related: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index`

---

### 2.2 GRCh38 Reference Genome

**Genome Reference Consortium Human Build 38 (GRCh38/hg38).**
Source used by this pipeline: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa`

> **How used:** Required by mosdepth for CRAM decoding (Step 1 of the pipeline). Downloaded
> automatically by `example/1000G_highcov/00_setup.sh`.

---

### 2.3 IGSR Sample Panel (Population Metadata)

International Genome Sample Resource (IGSR). Sample panel with population assignment,
superpopulation (AFR/AMR/EAS/EUR/SAS), and reported sex.
<https://www.internationalgenome.org/>
Download: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped`
Download FAQ: <https://www.internationalgenome.org/faq/what-tools-can-i-use-to-download-igsr-data/>

> **How used:** Downloaded by `example/1000G_highcov/00_setup.sh`. Used in
> `example/1000G_highcov/03_collect_qc.sh` to add population and superpopulation labels to the
> QC table, enabling overlay of population structure on PCA plots.

---

### 2.4 UK Biobank Exome Sequencing Targets (xGen)

UK Biobank. Exome capture target BED file (xgen panel) used for WES analyses.
Resource field 3801: <http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=3801>

**Sudlow, C., Gallacher, J., Allen, N., Beral, V., Burton, P., Danesh, J., … Collins, R. (2015).**
UK Biobank: An open access resource for identifying the causes of a wide range of complex
diseases of middle and old age. *PLOS Medicine*, 12(3), e1001779.
<https://doi.org/10.1371/journal.pmed.1001779>

> **How used:** Used to generate the WES-specific exclusion BED for UK Biobank analyses
> (`resources/GRCh38/UKB_WES/`). The xgen capture targets are buffered by 20 kb and combined
> with the WGS exclusion regions to mask non-target areas in WES coverage data.

---

## 3. Genomic Annotation and Exclusion Resources

### 3.1 Database of Genomic Variants (DGV)

**MacDonald, J. R., Ziman, R., Yuen, R. K. C., Feuk, L., & Scherer, S. W. (2014).** The
Database of Genomic Variants: a curated collection of structural variation in the human genome.
*Nucleic Acids Research*, 42(D1), D986–D992.
<https://doi.org/10.1093/nar/gkt958>
Database: <http://dgv.tcag.ca/>
GRCh38 release used: 2016-08-31 (`GRCh38_hg38_variants_2016-08-31.txt`)
hg19 release used: 2016-05-15 (`GRCh37_hg19_variants_2016-05-15.txt`)

> **How used:** Variant coordinates from DGV are included as one of the four exclusion layers in
> the pre-built exclusion BED files (`resources/GRCh38/` and `resources/hg19/`). Bins overlapping
> reported structural variant loci are excluded from PCA to reduce noise from true biological
> variation in the coverage signal.

---

### 3.2 10x Genomics Structural Variant Blacklist

10x Genomics. Structural variant blacklist BED files.
- GRCh38: `https://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed`
- hg19: `https://cf.10xgenomics.com/supp/genome/hg19/sv_blacklist.bed`

> **How used:** The SV blacklist is the first exclusion layer in the pre-built exclusion BED
> files (see `resources/GRCh38/generateExcludeBed.sh` and `resources/hg19/generateExcludeBed.hg19.sh`).
> It masks regions known to be problematic for structural variant analysis, including
> low-mappability and highly repetitive loci.

---

### 3.3 UCSC Genome Browser — Genomic Super-Duplications

**Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M., &
Haussler, D. (2002).** The human genome browser at UCSC. *Genome Research*, 12(6), 996–1006.
<https://doi.org/10.1101/gr.229102>
Browser: <https://genome.ucsc.edu/> (Table Browser → `genomicSuperDups` track)

> **How used:** Genomic super-duplication coordinates (`genomicSuperDups` table, hg38 and hg19)
> downloaded from the UCSC Table Browser and included as the fourth exclusion layer in the
> pre-built exclusion BED files.

---

### 3.4 GEM Mappability Track

**Marco-Sola, S., Sammeth, M., Guigó, R., & Ribeca, P. (2012).** The GEM mapper: fast, accurate
and versatile alignment by filtration. *Nature Methods*, 9(12), 1185–1188.
<https://doi.org/10.1038/nmeth.2221>
Software: <https://github.com/gemtools/gemtools>

> **How used:** GEM tools (`gem-indexer`, `gem-mappability`, `gem-2-wig`) were used to compute
> 50-mer (and 100/125/150-mer) mappability tracks for the reference genome. Regions with
> mappability < 1.0 form the second exclusion layer in the pre-built exclusion BED files
> (see `resources/GRCh38/generateExcludeBed.sh`).

---

## 4. Bioinformatics Software

### 4.1 mosdepth

**Pedersen, B. S., & Quinlan, A. R. (2018).** Mosdepth: quick coverage calculation for genomes
and exomes. *Bioinformatics*, 34(5), 867–868.
<https://doi.org/10.1093/bioinformatics/btx699>
GitHub: <https://github.com/brentp/mosdepth>
Version bundled in container: **v0.3.9**

> **How used:** Step 1 of the NGS-PCA pipeline. mosdepth computes per-bin sequencing coverage
> depth from BAM/CRAM files using fixed-width 1 kb bins (`--by 1000`). Its `*.regions.bed.gz`
> output files are the direct input to NGS-PCA. mosdepth v0.3.9 is bundled in the container
> image (`Dockerfile`) and is available in the pre-built container at `ghcr.io/jlanej/ngs-pca:latest`.

---

### 4.2 bedtools

**Quinlan, A. R., & Hall, I. M. (2010).** BEDTools: a flexible suite of utilities for comparing
genomic features. *Bioinformatics*, 26(6), 841–842.
<https://doi.org/10.1093/bioinformatics/btq033>
Website: <https://bedtools.readthedocs.io/>

> **How used:** `bedtools merge` and `bedtools sort` are used in the exclusion BED generation
> scripts (`resources/GRCh38/generateExcludeBed.sh`, `resources/hg19/generateExcludeBed.hg19.sh`,
> and `resources/GRCh38/UKB_WES/genUKB.exclude.sh`) to merge overlapping genomic intervals and
> produce the final sorted, merged exclusion BED files.

---

### 4.3 SAMtools

**Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009).**
The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079.
<https://doi.org/10.1093/bioinformatics/btp352>
Website: <https://www.htslib.org/>

> **How used:** Referenced in pipeline documentation for FASTA indexing (e.g., `samtools faidx`)
> when preparing the reference genome for mosdepth CRAM decoding.

---

### 4.4 IBM Aspera Connect

IBM. Aspera Connect — high-speed file transfer using the FASP protocol.
Download: <https://www.ibm.com/products/aspera/downloads>
Version bundled in container: **v4.2.12.780**

> **How used:** Optional (but recommended) high-speed download client in
> `example/1000G_highcov/01_download_and_mosdepth.sh` for retrieving CRAM files from the ENA
> FTP and the GRCh38 reference genome from the EBI 1000G FTP. Aspera FASP protocol typically
> achieves 10–100× faster transfers than plain FTP/HTTP. The pipeline falls back to `wget`
> automatically if `ascp` is not found.

---

### 4.5 Apptainer (formerly Singularity)

**Kurtzer, G. M., Sochat, V., & Bauer, M. W. (2017).** Singularity: Scientific containers for
mobility of compute. *PLOS ONE*, 12(5), e0177459.
<https://doi.org/10.1371/journal.pone.0177459>
Website: <https://apptainer.org/>

> **How used:** Recommended container runtime for HPC deployments. The pre-built container image
> (`ghcr.io/jlanej/ngs-pca:latest`) bundles NGS-PCA, mosdepth, and exclusion BED files so that
> no additional software installation is required. All four pipeline scripts in
> `example/1000G_highcov/` use Apptainer to invoke mosdepth and NGS-PCA.

---

### 4.6 Globus

Globus. Managed file transfer service.
Website: <https://www.globus.org/>

> **How used:** Mentioned in `example/1000G_highcov/README.md` as an alternative download method
> for sites where Aspera (FASP) ports are blocked. Requires manual setup outside the pipeline.

---

## 5. Java Libraries

### 5.1 Apache Commons Math

The Apache Software Foundation. Commons Math: The Apache Commons Mathematics Library.
<https://commons.apache.org/proper/commons-math/>
Version used: **3.6.1** (declared in `ngspca/pom.xml`)
License: **Apache License 2.0** — <https://github.com/apache/commons-math/blob/master/LICENSE.txt>
Full license text: <https://www.apache.org/licenses/LICENSE-2.0.txt>

> **How used:** Core numerical library throughout `ngspca/src/`. Provides:
> - `BlockRealMatrix` and `RealMatrix` — dense matrix representations for the coverage and SVD matrices
> - `SingularValueDecomposition` — standard SVD on the reduced matrix after random projection
> - `MersenneTwister` — pseudorandom number generator for the random test matrix (controlled by `-randomSeed`)
> - `Median` — robust median computation for normalization in `NormalizationOperations.java`
> - `QRDecomposition` — planned replacement for Jama QR (see `docs/future-work.md`)

---

### 5.2 EJML — Efficient Java Matrix Library

**Abeles, P. (2022).** Efficient Java Matrix Library (EJML).
<https://ejml.org/>
GitHub: <https://github.com/lessthanoptimal/ejml>
Version used: **0.30** (declared in `ngspca/pom.xml`)
License: **Apache License 2.0** — <https://github.com/lessthanoptimal/ejml/blob/master/LICENSE-2.0.txt>
Full license text: <https://www.apache.org/licenses/LICENSE-2.0.txt>

> **How used:** Dense matrix operations in the Java source. Declared as a dependency in `pom.xml`.

---

### 5.3 JAMA — Java Matrix Package

**Hicklin, J., Moler, C., Webb, P., Boisvert, R. F., Miller, B., Pozo, R., & Remington, K.**
JAMA: A Java Matrix Package.
<https://math.nist.gov/javanumerics/jama/>
Version used: **1.0.3** (declared in `ngspca/pom.xml`)
License: **Public Domain** — JAMA is a cooperative product of The MathWorks and the National
Institute of Standards and Technology (NIST) and is not subject to copyright protection in the
United States. See the project home page: <https://math.nist.gov/javanumerics/jama/#license>

> **How used:** Currently used for `QRDecomposition` inside the power iteration loop in
> `RandomizedSVD.java` (converting `RealMatrix` → `Matrix` → `QRDecomposition` → `RealMatrix`).
> Planned for removal in favor of Apache Commons Math's own `QRDecomposition` to eliminate
> redundant array copies (see `docs/future-work.md`, item 1).

---

### 5.4 Apache Commons CLI

The Apache Software Foundation. Commons CLI: Command Line Interface.
<https://commons.apache.org/proper/commons-cli/>
Version used: **1.4** (declared in `ngspca/pom.xml`)
License: **Apache License 2.0** — <https://github.com/apache/commons-cli/blob/master/LICENSE.txt>
Full license text: <https://www.apache.org/licenses/LICENSE-2.0.txt>

> **How used:** Parses command-line arguments in `CmdLine.java`, providing the `-input`,
> `-outputDir`, `-numPC`, `-bedExclude`, and all other user-facing parameters.

---

### 5.5 Apache Commons Lang

The Apache Software Foundation. Commons Lang.
<https://commons.apache.org/proper/commons-lang/>
Version used: **3.7** (declared in `ngspca/pom.xml`)
License: **Apache License 2.0** — <https://github.com/apache/commons-lang/blob/master/LICENSE.txt>
Full license text: <https://www.apache.org/licenses/LICENSE-2.0.txt>

> **How used:** General-purpose string and utility functions used across the Java source.

---

### 5.6 Apache Commons IO

The Apache Software Foundation. Commons IO.
<https://commons.apache.org/proper/commons-io/>
Version used: **2.6** (declared in `ngspca/pom.xml`)
License: **Apache License 2.0** — <https://github.com/apache/commons-io/blob/master/LICENSE.txt>
Full license text: <https://www.apache.org/licenses/LICENSE-2.0.txt>

> **How used:** File I/O utilities in `FileOps.java` (e.g., `FileUtils.writeStringToFile`).

---

### 5.7 HTSJDK

Broad Institute & collaborators. HTSJDK: A Java API for high-throughput sequencing data (HTS)
formats.
GitHub: <https://github.com/samtools/htsjdk>
Version used: **2.16.0** (declared in `ngspca/pom.xml`)
License: **MIT License** — <https://github.com/samtools/htsjdk/blob/master/LICENSE.txt>

> **How used:** Provides the `BEDFileReader` and `BEDFeature` classes used in
> `MosdepthUtils.java` and `BedUtils.java` to parse mosdepth BED output files and perform
> genomic-region overlap detection against the exclusion BED.

---

### 5.8 JUnit

JUnit Contributors. JUnit — Java Testing Framework.
<https://junit.org/>
Version used: **3.8.1** (declared in `ngspca/pom.xml`, test scope)
License: **IBM Common Public License v1.0 (CPL-1.0)** — used by JUnit 3.x releases.
Full license text: <https://opensource.org/licenses/cpl1.0.php>

> **How used:** Unit testing framework for `ngspca/src/test/java/org/pankratzlab/ngspca/AppTest.java`.

---

## 6. Build and Infrastructure

### 6.1 Apache Maven

The Apache Software Foundation. Apache Maven — Software Project Management Tool.
<https://maven.apache.org/>
Version required: **3.6+** (build stage base image: `maven:3.9-eclipse-temurin-11`)

> **How used:** Build system for the NGS-PCA Java project (`ngspca/pom.xml`). Produces the
> fat JAR via `maven-shade-plugin`. CI and release workflows invoke `mvn -B package`.

---

### 6.2 Eclipse Temurin (Java 11 JRE)

Eclipse Adoptium. Eclipse Temurin — OpenJDK Distribution.
<https://adoptium.net/>
Runtime image: `eclipse-temurin:11-jre-jammy`

> **How used:** Base runtime image in the multi-stage `Dockerfile`. Provides the Java 11 JRE
> used to run the NGS-PCA JAR inside the container.

---

### 6.3 GitHub Actions

GitHub. GitHub Actions — CI/CD platform.
<https://docs.github.com/en/actions>

> **How used:** Two workflows are defined in `.github/workflows/`:
> - `integration-test.yml` — runs the 18-sample 1000G chr1 example on every pull request and push to `main`/`master`, verifying output MD5 checksums.
> - `release.yml` — builds and publishes the fat JAR as a GitHub Release asset on every `v*` tag push.

---

### 6.4 GitHub Container Registry (GHCR)

GitHub. GitHub Packages — Container Registry.
<https://ghcr.io>
Image: `ghcr.io/jlanej/ngs-pca:latest`

> **How used:** Hosts the pre-built container image. Users pull the image via
> `apptainer pull ngs-pca.sif docker://ghcr.io/jlanej/ngs-pca:latest`.

---

### 6.5 SLURM Workload Manager

SchedMD. SLURM — Highly Scalable Cluster Management and Job Scheduling System.
<https://slurm.schedmd.com/>

> **How used:** The primary job scheduler targeted by the `example/1000G_highcov/` pipeline
> scripts. `01_download_and_mosdepth.sh` uses SLURM array jobs to parallelize download and
> coverage computation across 3,202 samples. `02_run_ngspca.sh` submits a single large-memory
> job for the PCA step. Documentation also covers adaptation for PBS/Torque, SGE, and LSF.

---

## 7. Needs Citation / Unverified

The following items are referenced in the codebase but lack a verified canonical citation or
publication. If you use these resources, please locate the current authoritative reference.

| Resource | Context | Status |
|---|---|---|
| Producer-consumer pattern reference (<https://dzone.com/articles/the-evolution-of-producer-consumer-problem-in-java>) | Code comment in `MosdepthUtils.java` for thread-pool design | Web article; no formal citation |
| 10x Genomics SV Blacklist | Exclusion BED layer; hosted at `http://cf.10xgenomics.com/supp/genome/` | No associated peer-reviewed publication identified; cite 10x Genomics directly if needed |
