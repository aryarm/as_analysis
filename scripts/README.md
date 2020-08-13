Counts of reads that overlap heterozygous SNPs are analyzed in `Snakefile-counts`. This section of the pipeline requires several custom-made scripts.

For those who don't have counts from DNA sequencing reads, we've included "RNA-only" versions of these scripts which have `-rna` appended to their file names.

### allele_imbalance.r
Defines a function for detecting allelic imbalance given count data. This code was adapted from [@Arkosen](https://github.com/Arkosen)

### benchmark.py
This python script can be used to summarize the runtime and memory usage of the pipeline. It is provided for the convenience of the user and is not officially part of the pipeline.

### find_imbalance.r
An R script that calls `allele_imbalance.r`. This script is called, in turn, by Snakefile-counts.

### prepare_counts.r
This R script prepares counts for use by `find_imbalance.r`, incorporating GQ scores and gene information for each SNP.

### remove_indels.r
This R script extracts SNPs from a text file containing read counts. Its output is used by `prepare_counts.r`.

### merge_results.py
This python script can be used to merge the final results of the pipeline into a single csv file. It is provided for the convenience of the user and is not officially part of the counts pipeline.
