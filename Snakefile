import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.1.4")

configfile: "config.yaml"


def read_samples():
    """Function to get names and dna/rna fastq paths from a sample file
    specified in the configuration. Input file is expected to have five
    columns:
    <unique_sample_id> <fastq1_dna_path> <fastq2_dna_path> <fastq1_rna_path> <fastq2_rna_path>
    Modify this function as needed to provide the correct dictionary."""
    f = open(config['sample_file'], "r")
    samp_dict = {}
    for line in f:
        words = line.strip().split("\t")
        samp_dict[words[0]] = ((words[1], words[2]), (words[3], words[4]))
    return samp_dict
GLOBAL_SAMP = read_samples()

# the user can change config['SAMP_NAMES'] here (or define it in the config
# file) to contain whichever sample names they'd like to run the pipeline on
if 'SAMP_NAMES' not in config:
    config['SAMP_NAMES'] = list(GLOBAL_SAMP.keys())


rule all:
    input:
        # if you'd like to run the pipeline on only a subset of the samples,
        # you should specify them in the config['SAMP_NAMES'] variable above
        expand(config['output_dir'] + "/final/{sample}/result.csv.gz",
               sample=config['SAMP_NAMES'])


# variant calling pipeline
SAMP1 = {samp: GLOBAL_SAMP[samp][0] for samp in config['SAMP_NAMES']}
include: "Snakefiles/Snakefile-variant_calling"
config['vcf_file'] = rules.filter_hets.output

# WASP pipeline
SAMP2 = {samp: GLOBAL_SAMP[samp][1] for samp in config['SAMP_NAMES']}
SAMP_TO_VCF_ID = {samp: samp for samp, fastqs in GLOBAL_SAMP.items()}
config['snp_h5_dir'] = config['output_dir'] + "/genotypes/snp_h5"
include: "Snakefiles/Snakefile-WASP"

# counts analysis pipeline
SAMP3 = {}
for samp in config['SAMP_NAMES']:
    dna = rules.rm_dups.output.final_bam.format(sample=samp)
    rna = rules.rmdup_pe.output.sort.format(sample=samp)
    SAMP3[samp] = (dna, rna)
include: "Snakefiles/Snakefile-counts"
