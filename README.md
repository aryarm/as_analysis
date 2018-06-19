[![Snakemake](https://img.shields.io/badge/snakemake-≥5.1.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# as_analysis
A pipeline for detecting allelic imbalance from DNA and RNA seq reads

# files
### Snakefile-variant_calling
A snakemake pipeline for running the variant calling portion of the pipeline. It generates a filtered VCF of heterozygous sites from DNA seq reads.

### Snakefile-WASP
A snakemake pipeline for running [WASP](https://github.com/bmvdgeijn/WASP). It generates an un-biased, filtered BAM file from RNA seq reads.

### Snakefile-counts
A snakemake pipeline for generating and analyzing counts of heterozygous sites and detecting significant allelic imbalance from DNA and RNA BAM files.

### config-[*].yaml
Required input for each snakemake pipeline

### allele_imbalance.r
Defines a function for detecting allelic imbalance given count data. This code was adapted from [@Arkosen](https://github.com/Arkosen/Detecting-structural-variants-/blob/master/allele_imbalance.r)

### find_imbalance.r
An R script that sets up the necessary data for calling `allele_imbalance.r` and then calls it. This script is called, in turn, by Snakefile-counts.

### run-all.bash
A bash script for executing the entire pipeline on an SGE cluster using snakemake.

# how to execute the pipeline
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/).

Once you've installed snakemake and adapted the provided config files to your liking, you can execute each section of the pipeline individually using the `snakemake` command or all at once on an SGE cluster using the provided [run-all.bash](https://github.com/aryam7/as_analysis/blob/master/run-all.bash) script.

If you have anaconda installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the following dependencies:
- [pytables](https://www.pytables.org/) (v2)
- [pysam](https://github.com/pysam-developers/pysam)
- [GATK](https://software.broadinstitute.org/gatk/gatk4) (v4)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [tabix](https://github.com/samtools/tabix)

Whether or not you use the `--use-conda` option, you must manually install a number of R packages. The full list of such packages is at the top of the [find_imbalance.r](https://github.com/aryam7/as_analysis/blob/master/find_imbalance.r) script.
