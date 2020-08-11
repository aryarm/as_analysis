import warnings
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.18.0")

if 'imported' not in config:
    configfile: "config.yaml"


def check_config(value, default=False, place=config):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

no_variant_calling = check_config('vcf_file')

def read_samples():
    """Function to get names and dna/rna fastq paths from a sample file
    specified in the configuration. Input file is expected to have five
    columns:
    <unique_sample_id> <fastq1_dna_path> <fastq2_dna_path> <fastq1_rna_path> <fastq2_rna_path>
    Or if the vcf_file config option is provided, the input file is
    expected to have these five columns:
    <vcf_sample_id> <unique_sample_name> <dna_bam_path> <fastq1_rna_path> <fastq2_rna_path>
    or (optionally) just these four:
    <vcf_sample_id> <unique_sample_name> <fastq1_rna_path> <fastq2_rna_path>
    Modify this function as needed to provide the correct dictionary."""
    f = open(config['sample_file'], "r")
    samp_dict = {}
    for line in f:
        words = line.strip().split("\t")
        if no_variant_calling:
            if len(words) == 5:
                samp_dict[words[1]] = (words[2], (words[3], words[4]))
            elif len(words) == 4:
                # then the dna_bam_path hasn't been specified
                samp_dict[words[1]] = ((words[2], words[3]),)
                # sanity check to make sure that rna_only is set to true
                config['rna_only'] = True
            else:
                raise ValueError("The samples file is not formatted correctly for this line: "+line)
        else:
            samp_dict[words[0]] = ((words[1], words[2]), (words[3], words[4]))
    return samp_dict
GLOBAL_SAMP = read_samples()

# the user can define config['SAMP_NAMES'] to contain whichever sample names
# they'd like to run the pipeline on
if not check_config('SAMP_NAMES'):
    config['SAMP_NAMES'] = list(GLOBAL_SAMP.keys())
else:
    # double check that the user isn't asking for samples they haven't provided
    user_samps_len = len(config['SAMP_NAMES'])
    config['SAMP_NAMES'] = list(set(GLOBAL_SAMP.keys()).intersection(config['SAMP_NAMES']))
    if len(config['SAMP_NAMES']) != user_samps_len:
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")



rule all:
    input:
        # if you'd like to run the pipeline on only a subset of the samples,
        # you should specify them in the config['SAMP_NAMES'] variable above
        expand(config['output_dir'] + "/final/{sample}/result.csv.gz",
               sample=config['SAMP_NAMES'])

# variant calling pipeline
if not no_variant_calling:
    SAMP1 = {samp: GLOBAL_SAMP[samp][0] for samp in config['SAMP_NAMES']}
    include: "Snakefiles/Snakefile-variant_calling"
    config['vcf_file'] = rules.filter_hets.output.vcf

# WASP pipeline
if no_variant_calling:
    SAMP2 = {
        samp: GLOBAL_SAMP[samp][1] if len(GLOBAL_SAMP[samp]) == 2 else GLOBAL_SAMP[samp][0]
        for samp in config['SAMP_NAMES']
    }
    def read_vcf_samples():
        f = open(config['sample_file'], "r")
        samp_dict = {}
        for line in f:
            words = line.strip().split("\t")
            if words[1] in config['SAMP_NAMES']:
                samp_dict[words[1]] = words[0]
        return samp_dict
    SAMP_TO_VCF_ID = read_vcf_samples()
else:
    SAMP2 = {samp: GLOBAL_SAMP[samp][1] for samp in config['SAMP_NAMES']}
    SAMP_TO_VCF_ID = {samp: samp for samp, fastqs in GLOBAL_SAMP.items()}
include: "Snakefiles/Snakefile-WASP"

# counts analysis pipeline
SAMP3 = {}
for samp in config['SAMP_NAMES']:
    if no_variant_calling:
        if len(GLOBAL_SAMP[samp]) == 2:
            dna = GLOBAL_SAMP[samp][0]
        else:
            config['rna_only'] = True
    else:
        dna = rules.rm_dups.output.final_bam.format(sample=samp)
    rna = rules.rmdup.output.sort.format(sample=samp)
    if config['rna_only']:
        SAMP3[samp] = (rna,)
    else:
        SAMP3[samp] = (dna, rna)
include: "Snakefiles/Snakefile-counts"
