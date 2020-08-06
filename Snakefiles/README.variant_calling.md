# Variant Calling Pipeline

If you'd like to run the variant calling pipeline on its own, you should provide required input in [configs/config-variant_calling.yaml](/configs/config-variant_calling.yaml).

## Inputs
 - FASTQ files containing DNA sequencing reads for each sample. Specify the location of these files in a tab delimited text file containing three columns: `unique_sample_name | dna_fastq_path_1 | dna_fastq_path_2` where each row is a different sample.
 - A reference genome for DNA-seq alignment
 - A BWA index of the aforementioned reference genome
 
## Other Inputs and Options
- If you choose to use GATK's base quality ([BQSR](https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)) and variant quality ([VQSR](https://software.broadinstitute.org/gatk/documentation/article.php?id=39)) score recalibration, you must provide the following, as described in [this GATK article](https://software.broadinstitute.org/gatk/documentation/article.php?id=1259):

    - True sites training resource: HapMap
    - True sites training resource: Omni
    - Non-true sites training resource: 1000G
    - Known sites resource, not used in training: dbSNP
    
    These resources can usually be easily obtained from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle).

    You will also be required to specify a target sensitivity value, as described in a [previously mentioned GATK article](https://software.broadinstitute.org/gatk/documentation/article?id=39).
 
 - If you choose not to perform VQSR, the pipeline will default to [hard filtering](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set) your variants. You will need to provide a [GATK filter expression](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php#--filter-expression), as described in [this GATK article](https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). One example might be ```"QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"```.

## Running the variant calling pipeline on its own
When calling [Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html), use options `-s` and `--configfile` to specify the location of the Snakefile and its corresponding config file. We also recommend using the `--use-conda` option to let Snakemake [handle all dependencies](http://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) of the pipeline.

    snakemake -s Snakefiles/Snakefile-variant_calling --configfile configs/config-variant_calling.yaml --use-conda

## Output
The variant calling pipeline creates the following directories under the output directory specified in your config file. The `genotypes` folder will contain the final output, a filtered VCF containing heterozygous SNPs for all samples.
 - dna_align - output from BWA and samtools
 - base_recal - output from GATK's BQSR
 - haplotype - output from GATK's Haplotype Caller and a file `ALL.genotype.vcf.gz` containing genotyped variants for all samples
 - variant_filter - output from variant filtering (either VQSR or hard filtering) of SNPs
 - genotypes - heterozgyotes from the filtered VCF
