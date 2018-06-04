# as_analysis
A pipeline for detecting allelic imbalance from DNA and RNA seq reads

# Files
### final.bash
A pipeline for getting counts of heterozygous sites, written in bash

### Snakefile
A snakemake pipeline for running final.bash and then analyzing counts of heterozygous sites and detecting significant allelic imbalance

### config.yaml
Required input for the snakemake pipeline

### allele_imbalance.r
Defines a function for detecting allelic imbalance given count data. This code was adapted from [@Arkosen](https://github.com/Arkosen/Detecting-structural-variants-/blob/master/allele_imbalance.r)

### stats.r
An R script that sets up the necessary data for calling allele_imbalance.r
