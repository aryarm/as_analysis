Individual portions of the snakemake pipeline can also be run on their own. We have separate Snakefiles (and corresponding config files) for each of the following three portions of the pipeline:

### [Variant Calling](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.variant_calling.md)
A snakemake pipeline for running the variant calling portion of the pipeline. It generates a filtered VCF of heterozygous SNPs from DNA seq reads. See the variant calling README for more information.

### [WASP](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.WASP.md)
A snakemake pipeline for running [WASP](https://github.com/bmvdgeijn/WASP). It generates an un-biased, filtered BAM file from RNA seq reads. See the WASP README for more information.

### [Analyzing Read Counts](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.counts.md)
A snakemake pipeline for generating and analyzing counts of heterozygous SNPs and detecting genes with significant allelic imbalance from DNA and RNA BAM files. See the counts README for more information.