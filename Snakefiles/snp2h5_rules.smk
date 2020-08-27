""" rules that are common to both the counts and WASP pipeline """

# set the wasp_dir if the user hasn't specified one
if 'wasp_dir' not in config:
    config['wasp_dir'] = ".snakemake/WASP"
# set the default SNP H5 dir if the user hasn't specified one
if 'snp_h5_dir' not in config:
    config['snp_h5_dir'] = config['output_dir'] + "/genotypes/snp_h5"

rule get_WASP:
    """Download WASP if doesn't exist"""
    output:
        find_intersecting_snps_script = config['wasp_dir'] + "/mapping/find_intersecting_snps.py",
        filter_remapped_reads_script = config['wasp_dir'] + "/mapping/filter_remapped_reads.py",
        rmdup_script = config['wasp_dir'] + "/mapping/rmdup_pe.py",
        bam2h5_script = config['wasp_dir'] + "/CHT/bam2h5.py",
        makefile = config['wasp_dir'] + "/snp2h5/Makefile"
    conda: "../envs/default.yaml"
    benchmark: config['output_dir'] + "/benchmark/snp2h5/get_WASP/all.tsv"
    shell:
        "[ ! -d \"{config[wasp_dir]}\" ] && "
        "curl -Ls https://api.github.com/repos/bmvdgeijn/WASP/tarball | "
        "tar zxf - -C \"{config[wasp_dir]}\" --strip-components 1"

rule install_WASP:
    """Make WASP snp2h5 if it hasn't been compiled yet"""
    input:
        ancient(rules.get_WASP.output.makefile)
    output:
        snp2h5_script = config['wasp_dir'] + "/snp2h5/snp2h5"
    conda: "../envs/default.yaml"
    benchmark: config['output_dir'] + "/benchmark/snp2h5/install_WASP/all.tsv"
    shell:
        "conda_path=\"$CONDA_PREFIX\" && "
        "if [ -z $conda_path ]; then "
        "conda_path=\"$(dirname $(dirname $(which conda)))\"; "
        "fi && [ ! -z $conda_path ] && "
        "(cd \"{config[wasp_dir]}/snp2h5\" && "
        "make -s HDF_INSTALL=\"$conda_path\")"

checkpoint vcf_chroms:
    """get the chroms from a VCF"""
    input:
        vcf = config['vcf_file'],
        vcf_index = config['vcf_file']+".tbi"
    output:
        config['output_dir']+"/genotypes/chroms.txt"
    conda: "../envs/default.yaml"
    benchmark: config['output_dir'] + "/benchmark/snp2h5/vcf_chroms/all.tsv"
    resources:
        mem_mb=200
    shell:
        "tabix --list-chroms {input.vcf} > {output}"

rule split_vcf_by_chr:
    """Split the provided VCF file by chromosome and gzip it for WASP"""
    input:
        vcf = config['vcf_file'],
        vcf_index = config['vcf_file']+".tbi",
    output:
        temp(config['output_dir'] + "/genotypes/bychrom/ALL.{chr}.vcf.gz")
    conda: "../envs/default.yaml"
    benchmark: config['output_dir'] + "/benchmark/snp2h5/split_vcf_by_chr/{chr}.tsv"
    resources:
        mem_mb=200
    shell:
        "tabix -h {input.vcf} {wildcards.chr} | bgzip > {output}"

def get_split_vcf(wildcards):
    with checkpoints.vcf_chroms.get().output[0].open() as f:
        return expand(
            rules.split_vcf_by_chr.output,
            chr=filter(lambda x: len(x), f.read().split('\n'))
        )

rule vcf2h5:
    """Convert VCF data files to HDF5 format"""
    input:
        chrom = config['chrom_info'],
        vcfs = get_split_vcf,
        snp2h5_script = ancient(rules.install_WASP.output.snp2h5_script)
    output:
        snp_index = config['snp_h5_dir'] + "/snp_index.h5",
        snp_tab = config['snp_h5_dir'] + "/snp_tab.h5",
        haplotype = config['snp_h5_dir'] + "/haplotype.h5"
    conda: "../envs/default.yaml"
    benchmark: config['output_dir'] + "/benchmark/snp2h5/vcf2h5/all.tsv"
    resources:
        mem_mb=2000
    shell:
        "{input.snp2h5_script} "
            "--chrom {input.chrom} "
            "--format vcf "
            "--snp_index {output.snp_index} "
            "--snp_tab {output.snp_tab} "
            "--haplotype {output.haplotype} "
            "{input.vcfs}"
