[![Snakemake](https://img.shields.io/badge/snakemake-≥5.1.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# as_analysis
A Snakemake pipeline for detecting allelic imbalance from DNA and RNA seq reads
![Pipeline Skeleton](https://drive.google.com/uc?export=view&id=1xefUeBPLLKfFfn9vu_IY2PHUBfC8XmST)

# files

### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis. It uses DNA and RNA FASTQ files to generate a summary of allelic imbalance for each gene.

### config.yaml
Defines options and input for the Snakemake pipeline.

### run-all.bash
An example bash script for executing the entire pipeline on an SGE cluster using snakemake.

# execution
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See the [run-all.bash](https://github.com/aryam7/as_analysis/blob/master/run-all.bash) script for an example. Make sure to provide required input and options in the [config file](https://github.com/aryam7/as_analysis/blob/master/config.yaml) before executing. For more information about what is required in the config file, see the [READMEs for each portion of the pipeline](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.md).

The entire pipeline is made up of three different sections. We provide a single [Snakefile](https://github.com/aryam7/as_analysis/blob/master/Snakefile) to execute all of them at once, but you can also execute each of these sections on their own. For each section that you'd like to run separately, you must fill out a new config file. You can find more information about these individual portions of the pipeline and how to execute them in the [Snakefiles directory](https://github.com/aryam7/as_analysis/tree/master/Snakefiles).

Snakemake automatically detects files that it has already generated, so that it doesn't rerun steps it doesn't have to. You can take advantage of this to easily skip any of the steps of our pipeline. Just provide your files in place of the files that Snakemake would have generated. Make sure that you place them in the correct directories and name them appropriately.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# dependencies
If you have [Anaconda](https://conda.io/docs/user-guide/install/index.html) installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the following dependencies:
- [GATK](https://software.broadinstitute.org/gatk/gatk4) (v4.x)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [tabix](https://github.com/samtools/tabix)
- [bcftools](http://www.htslib.org/download/)
- [plyr](https://www.rdocumentation.org/packages/plyr)
- [rmutil](https://www.rdocumentation.org/packages/rmutil/)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [dplyr](https://www.rdocumentation.org/packages/dplyr)

Note that our pipeline uses [WASP](https://github.com/bmvdgeijn/WASP) (v3.x), which is not available from Anaconda. If the pipeline is unable to locate WASP, it will automatically install it. The `--use-conda` option will automatically handle WASP's dependencies, but you can see the [WASP README](https://github.com/bmvdgeijn/WASP/blob/master/README.md) if you'd like to install these dependencies manually.
