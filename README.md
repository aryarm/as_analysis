# as_analysis
A pipeline for detecting allelic imbalance from DNA and RNA seq reads

# files
### final.bash
A pipeline for getting counts of heterozygous sites, written in bash

### Snakefile-variant_calling
A snakemake pipeline for running the variant calling portion of final.bash. It generates a filtered VCF of heterozygous sites from DNA seq reads.

### Snakefile-WASP
A snakemake pipeline for running [WASP](https://github.com/bmvdgeijn/WASP). It generates an un-biased, filtered BAM file from RNA seq reads.

### Snakefile-counts
A snakemake pipeline for generating and analyzing counts of heterozygous sites and detecting significant allelic imbalance from DNA and RNA BAM files.

### config.yaml
Required input for the snakemake pipeline

### allele_imbalance.r
Defines a function for detecting allelic imbalance given count data. This code was adapted from [@Arkosen](https://github.com/Arkosen/Detecting-structural-variants-/blob/master/allele_imbalance.r)

### find_imbalance.r
An R script that sets up the necessary data for calling allele_imbalance.r and calls it. This script is called, in turn, by Snakefile-counts.

### run-all.bash
A bash script for executing the entire pipeline on an SGE cluster using snakemake.
