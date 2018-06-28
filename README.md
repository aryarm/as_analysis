[![Snakemake](https://img.shields.io/badge/snakemake-≥5.1.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# as_analysis
A pipeline for detecting allelic imbalance from DNA and RNA seq reads

# files

### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis. It uses DNA and RNA FASTQ files to generate a summary of allelic imbalance for each gene.

### config.yaml
Defines options and input for the Snakemake pipeline.

### run-all.bash
An example bash script for executing the entire pipeline on an SGE cluster using snakemake.

# how to execute the pipeline
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See `run-all.bash` for an example. Make sure to provide required input and options in the [config](https://github.com/aryam7/as_analysis/blob/master/config.yaml) file before executing.

The entire pipeline is made up of three different sections. Each of these sections can be executed on their own but require you to fill out separate config files for each of them. See the [Snakefiles README](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.md) for more information.

If you have [Anaconda](https://conda.io/docs/user-guide/install/index.html) installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the following dependencies:
- [pytables](https://www.pytables.org/) (v2)
- [pysam](https://github.com/pysam-developers/pysam)
- [GATK](https://software.broadinstitute.org/gatk/gatk4) (v4)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [tabix](https://github.com/samtools/tabix)

Note that our pipeline uses [WASP](https://www.biorxiv.org/content/early/2014/11/07/011221), which is not downloadable from Anaconda. You must [manually install it](https://github.com/bmvdgeijn/WASP).

Whether or not you use the `--use-conda` option, you must manually install a number of R packages. The full list of such packages is at the top of the [find_imbalance.r](https://github.com/aryam7/as_analysis/blob/master/scripts/find_imbalance.r) script.
