# Pipeline for Analyzing Read Counts

If you'd like to run the counts pipeline on its own, you should provide required input in [configs/config-counts.yaml](/configs/config-counts.yaml).

## Inputs
 - BAM files containing DNA and RNA sequencing reads for each sample. Specify the location of these files in a tab delimited text file containing four columns: `vcf_sample_id | unique_sample_name | dna_bam_path | rna_bam_path` where each row is a different sample.
 
     If you don't have DNA sequencing reads, you can still run the pipeline using a more conservative version of our ASE test by setting the [`rna_only` config option](/configs/config-counts.yaml#L16) to true.
     
     Sometimes, a VCF may contain unique identifiers for each sample that differ from the sample name. For this reason, we ask that you provide a `vcf_sample_id` column so the pipeline can map between the sample name and the VCF sample ID. If this is not your situation, you can simply duplicate the `unique_sample_name` column as the `vcf_sample_id`.
 - A VCF file containing SNPs at which you'd like to get counts. If you don't have this, you can generate it using the [variant calling pipeline](/Snakefiles/README.variant_calling.md). The VCF can contain any type of variant, but if you are performing our allele specific analysis, it is much faster to provide a VCF containing only heterozygous SNPs, since all other variants will be ignored anyway. You should also make sure that your VCF has [GQ scores](https://gatk.broadinstitute.org/hc/en-us/articles/360035531392), as they are used in our test for allele specific expression.

## Other Inputs and Options
 - If your VCF is not in the hdf5 format, you must provide a text file containing names and lengths of all chromosomes in the assembly. You can usually download these from the [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/) or use [the example that comes with WASP](https://github.com/bmvdgeijn/WASP/blob/master/examples/example_data/chromInfo.hg19.txt).
 - If your VCF is not in the hdf5 format, your VCF must be converted to the HDF5 format before it can be used by WASP. You can optionally specify a directory to which you'd like these files written. Otherwise, the pipeline will default to `genotypes/snp_h5/`.
     
     If your VCF _has_ been converted to the hdf5 format, you should specify the directory of your hdf5 files as the `snp_h5_dir`. The `geno_probs.h5` files are not required. Note: you will still need to provide a VCF file.
 - Gene annotation information in the GTF format. These can ususally be obtained from GENCODE. For example, annotations for hg19 can be downloaded from [here](https://www.gencodegenes.org/releases/19.html).

## Running the counts pipeline on its own
Make sure to download [WASP](https://github.com/bmvdgeijn/WASP) before running the pipeline. In your config file, you must specify the path to the directory in which you downloaded it.

If you don't provide the `--use-conda` option, you must also download several R dependencies manually. The full list of these packages is at the top of the [find_imbalance.r](/scripts/find_imbalance.r) script.

When calling [Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html), use options `-s` and `--configfile` to specify the location of the Snakefile and its corresponding config file. We also recommend using the `--use-conda` option to let Snakemake [handle all dependencies](http://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) of the pipeline.

    snakemake -s Snakefiles/Snakefile-counts --configfile configs/config-counts.yaml --use-conda

## Output
The counts pipeline creates the following directories under the output directory specified in your config file. The `final` folder will contain the final output, a list of the resulting genes with allelic imbalance and their confidence.
 - genotypes - your input VCF, split by chromosome (if your VCF isn't already in the HDF5 format)
 - extract_gq - tables containing GQ scores for each of the samples in the input VCF
 - as_counts - allele-specific dna and rna read counts for all SNPs in the VCF from the output of [the `--txt_counts` argument of WASP's bam2h5.py script](https://github.com/bmvdgeijn/WASP/blob/master/CHT/README.bam2h5.md)
 - final - prepared counts and results of running the [allele_imbalance script](/scripts/allele_imbalance.r): genes with allelic imbalance and their FDR corrected p-value
