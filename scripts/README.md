Counts of reads that overlap heterozygous SNPs are analyzed in `Snakefile-counts`. This section of the pipeline requires several custom-made scripts.

### allele_imbalance.r
Defines a function for detecting allelic imbalance given count data. This code was adapted from [@Arkosen](https://github.com/Arkosen/Detecting-structural-variants-/blob/master/allele_imbalance.r)

### find_imbalance.r
An R script that calls `allele_imbalance.r`. This script is called, in turn, by Snakefile-counts.

### prepare_counts.r
This R script prepares counts for use by `find_imbalance.r`, incorporating GQ scores and gene information for each SNP.

### remove_indels.r
This R script extracts SNPs from a text file containing read counts. Its output is used by `prepare_counts.r`.
