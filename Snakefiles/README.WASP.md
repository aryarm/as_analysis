# WASP Pipeline

If you'd like to run the WASP pipeline on its own, you should provide required input in `config-WASP.yaml`.

This pipeline is an adaption of the one by [Graham McVicker](https://github.com/gmcvicker). You can find the original in the WASP repository [here](https://github.com/bmvdgeijn/WASP/blob/master/mapping/Snakefile). Note that our pipeline differs from the original in several ways, the most notable being that reads that don't overlap a SNP are discarded in our pipeline.

## Inputs
 - FASTQ files containing RNA sequencing reads for each sample. Specify the location of these files in a tab delimited text file containing four columns: `vcf_sample_id | unique_sample_name | rna_fastq_path_1 | rna_fastq_path_2` where each row is a different sample.
 
     Sometimes, a VCF may contain unique identifiers for each sample that differ from the sample name. For this reason, we ask that you provide a `vcf_sample_id` column so the pipeline can map between the sample name and the VCF sample ID. If this is not your situation, you can simply duplicate the `unique_sample_name` column as the `vcf_sample_id`.
 - A VCF file containing SNPs for WASP. If you don't have this, you can generate it using the [variant calling pipeline](https://github.com/aryam7/as_analysis/blob/master/Snakefiles/README.variant_calling.md). The VCF can contain any type of variant, but if you are performing our allele specific analysis, it is much faster to provide a VCF containing only heterozygous SNPs, since all other variants will be ignored anyway.
 - A reference genome for RNA-seq alignment
 - A STAR index of the aforementioned reference genome

## Other Inputs and Options
 - A text file containing names and lengths of all chromosomes in the assembly. You can usually download these from the [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/) or use [the example that comes with WASP](https://github.com/bmvdgeijn/WASP/blob/master/examples/example_data/chromInfo.hg19.txt).
 - Your VCF must be converted to the HDF5 format before it can be used by WASP. You can optionally specify a directory to which you'd like these files written. Otherwise, the pipeline will default to `genotypes/snp_h5/`.

## Running the WASP pipeline on its own
Make sure to download [WASP](https://github.com/bmvdgeijn/WASP) before running the pipeline. In your config file, you must specify the path to the directory in which you downloaded it.

When calling [Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html), use options `-s` and `--configfile` to specify the location of the Snakefile and its corresponding config file. We also recommend using the `--use-conda` option to let Snakemake [handle all dependencies](http://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) of the pipeline.

    snakemake -s Snakefiles/Snakefile-WASP --configfile Snakefiles/config-WASP.yaml --use-conda

## Output
The WASP pipeline creates the following directories under the output directory specified in your config file. The `rmdup` folder will contain the final output, filtered BAM files containing RNA reads for each sample.
 - genotypes - your input VCF, split by chromosome
 - genotypes/snp_h5 - your input VCF in HDF5 format
 - map1 - results from first mapping of reads to genome
 - map1_sort - sorted BAM files from first mapping
 - find_intersecting_snps - results from find_intersecting_snps.py
 - map2 - results from second mapping of reads to genome
 - map2_sort - sorted BAM files from second mapping
 - filter_remapped_reads - results from filter_remapped_reads.py
 - rmdup - Final BAM files with duplicate reads removed. Only *sort* files need to be kept.
