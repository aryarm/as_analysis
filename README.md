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
Defines a function for detecting allelic imbalance given count data

### stats.r
An R script that sets up the necessary data for calling allele_imbalance.r
