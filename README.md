[![Snakemake](https://img.shields.io/badge/snakemake-≥5.1.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# as_analysis
A Snakemake pipeline for detecting allelic imbalance from DNA and RNA seq reads

# files

### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis. It uses DNA and RNA FASTQ files to generate a summary of allelic imbalance for each gene.

### config.yaml
Defines options and input for the Snakemake pipeline.

### run-all.bash
An example bash script for executing the entire pipeline on an SGE cluster using snakemake.

# execution
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See the [run-all.bash](https://github.com/aryam7/as_analysis/blob/master/run-all.bash) script for an example. Make sure to provide required input and options in the [config file](https://github.com/aryam7/as_analysis/blob/master/config.yaml) before executing. For more information about what is required in the config file, see the [READMEs for each portion of the pipeline](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.md).

The entire pipeline is made up of three different sections. We provide a single Snakefile to execute all of them at once, but you can also execute each of these sections on their own. For each section that you'd like to run separately, you must fill out a new config file. You can find more information about these portions of the pipeline and how to execute them in the [Snakefiles directory](https://github.com/aryam7/as_analysis/tree/master/Snakefiles).

# dependencies
If you have [Anaconda](https://conda.io/docs/user-guide/install/index.html) installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the following dependencies:
- [GATK](https://software.broadinstitute.org/gatk/gatk4) (v4.x)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [tabix](https://github.com/samtools/tabix)
- [bcftools](http://www.htslib.org/download/)

Note that our pipeline uses [WASP](https://www.biorxiv.org/content/early/2014/11/07/011221) (v3.x), which is not available from Anaconda. You must [manually install it](https://github.com/bmvdgeijn/WASP). Our pipeline automatically handles WASP's dependencies, but you can see the [WASP README](https://github.com/bmvdgeijn/WASP/blob/master/README.md) if you'd like to install these dependencies manually.

Regardless of whether you use the `--use-conda` option, you must manually install a number of R packages. The full list of such packages is at the top of the [find_imbalance.r](https://github.com/aryam7/as_analysis/blob/master/scripts/find_imbalance.r) script.
