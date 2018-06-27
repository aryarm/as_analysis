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
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See `run-all.bash` for an example. Make sure to provide required input and options in the config file.

If you have anaconda installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the following dependencies:
- [pytables](https://www.pytables.org/) (v2)
- [pysam](https://github.com/pysam-developers/pysam)
- [GATK](https://software.broadinstitute.org/gatk/gatk4) (v4)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [tabix](https://github.com/samtools/tabix)

Whether or not you use the `--use-conda` option, you must manually install a number of R packages. The full list of such packages is at the top of the [find_imbalance.r](https://github.com/aryam7/as_analysis/blob/master/scripts/find_imbalance.r) script.
