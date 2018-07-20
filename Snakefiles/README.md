We provide a single [Snakefile](https://github.com/aryam7/as_analysis/blob/master/Snakefile) to execute the entire pipeline at once, but individual portions of the snakemake pipeline can also be run on their own. We have separate Snakefiles (and corresponding config files) for each of the following three portions of the pipeline:

### [Variant Calling](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/Snakefile-variant_calling)
A snakemake pipeline for running the variant calling portion of the pipeline. It generates a filtered VCF of heterozygous SNPs from DNA FASTQs. See our [variant calling README](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.variant_calling.md) for more information.

### [WASP](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/Snakefile-WASP)
A snakemake pipeline for running [WASP](https://github.com/bmvdgeijn/WASP). It generates an un-biased, filtered BAM file (containing only reads that overlap SNPs) from RNA FASTQs. See our [WASP README](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.WASP.md) for more information.

### [Analyzing Read Counts](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/Snakefile-counts)
A snakemake pipeline for generating and analyzing counts of heterozygous SNPs and detecting genes with significant allelic imbalance from DNA and RNA BAM files. See our [counts README](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.counts.md) for more information.
