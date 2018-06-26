Individual portions of the snakemake pipeline can also be run on their own. We have separate Snakefiles (and corresponding config files) for each of the following three portion of the pipeline:
	
1) Snakefile-variant_calling - align DNA fastq's and generate a filtered VCF containing heterozygous SNPs

2) Snakefile-WASP - align RNA fastq's and filter using WASP to reduce mapping bias

3) Snakefile-counts - retrieve counts of reads overlapping SNPs from BAM files generated in each of the previous portions of the pipeline and use them to find genes which demonstrate allelic imbalance
