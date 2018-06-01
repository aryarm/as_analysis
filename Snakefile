configfile: "config.yaml"

import glob


def read_samples():
    """Function to get names and fastq paths from a sample file specified
    in the configuration. Input file is expected to have 4 columns:
    <1000genomes_id> <unique_sample_id> <fastq1_path> <fastq2_path>. Modify
    this function as needed to provide a dictionary of sample_id keys and
    (fastq1, fastq1) values"""
    f = open(config['sample_file'], "r")
    samp_dict = {}
    for line in f:
        words = line.strip().split()
        samp_dict[words[1]] = ((words[2], words[3]), (words[4], words[5]))

    return samp_dict


def read_1kg_samples():
    f = open(config['sample_file'], "r")
    samp_dict = {}
    for line in f:
        words = line.strip().split()
        samp_dict[words[1]] = words[0]

    return samp_dict


SAMP_TO_1KG = read_1kg_samples()


def get_chromosomes():
    """Gets list of chromosomes from VCF files"""
    filenames = os.listdir(config['vcf_dir'])

    chr_names = set([])
    for filename in filenames:
        if filename.endswith(".vcf.gz"):
            m = re.match(".*(chr[0-9A-Z]+).*", filename)
            if m:
                chr_names.add(m.groups()[0])

    return chr_names


config['output_dir'] = "out"
config['num_threads'] = 24


def set_wasp_config():
    """WASP needs some config values. Let's load them, so the user doesn't have
    to."""
    config['vcf_dir'] = "genotypes"
    config['snp_h5_dir'] = "genotypes/snp_h5"




rule all:
    input:
        expand(config['output_dir'] + "/as_counts/{sample}.as_counts.txt.gz",
               sample=read_samples().keys())

rule align_dna:
    """Align DNA reads using BWA-MEM. Note that we use -R to specify read group
    info for haplotype caller."""
    input:
        fastq1 = lambda wildcards: read_samples()[wildcards.sample][0][0],
        fastq2 = lambda wildcards: read_samples()[wildcards.sample][0][1]
    output:
        pipe("{sample}.sam")
    threads: config['num_threads']
    shell:
        "bwa mem -M "
        "-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA' "
        "-t {threads} "
        "{input.fastq1} {input.fastq2} > {output}"

rule sam_to_bam:
    """Convert a SAM file to its more compressed counterpart. Note that -u to
    create an uncompressed bam file. We use -q to filter alignments with MAPQ
    scores less than 20."""
    input:
        rules.align_dna.output
    output:
        pipe("{sample}.bam")
    threads: config['num_threads']
    shell:
        "samtools view -u -b -F 4 -q 20 -@ {threads} <{input} >{output}"

rule sort_bam_by_name:
    """Sort the bam output by name (not by coordinates yet)"""
    input:
        rules.sam_to_bam.output
    output:
        pipe("{sample}.nameSort.bam")
    threads: config['num_threads']
    shell:
        "samtools sort -n -@ {threads} <{input} >{output}"

rule add_mate_info:
    """Use fixmate to fill in mate coordinates and mate related flags, since
    our data is pair-ended. We need the MC tags (included because we used the
    -m flag) that it creates for markdup"""
    input:
        rules.sort_bam_by_name.output
    output:
        pipe("{sample}.mate.nameSort.bam")
    threads: config['num_threads']
    shell:
        "samtools fixmate -m -@ {threads} {input} {output}"

rule sort_bam_by_coord:
    """Sort the bam output by coordinates. Needed for markdup use later on."""
    input:
        rules.add_mate_info.output
    output:
        pipe("{sample}.coordSort.mate.nameSort.bam")
    threads: config['num_threads']
    shell:
        "samtools sort -@ {threads} <{input} >{output}"

rule rm_dups:
    """Remove duplicates that may have occurred from PCR and index the
    resulting file."""
    input:
        rules.sort_bam_by_coord.output
    output:
        final_bam = config['output_dir'] + "/dna_align/{sample}.final.bam"
        final_bam_index = config['output_dir'] + "/dna_align/{sample}.final.bam.bai"
    threads: config['num_threads']
    shell:
        "mkdir -p {config[output_dir]}/dna_align;"
        "samtools markdup -@ {threads} {input} {output.final_bam}; "
        "samtools -b -@ {threads} {output.final_bam_index}"


rule base_recal:
    """Recalibrate the base quality scores. They might be biased"""
    input:
        ref = config['ref_genome']
        bam = rules.rm_dups.output.final_bam
        known_sites = config['dbSNP']
    output:
        config['output_dir'] + "/base_recal/{sample}.recal_data.table"
    shell:
        "mdkir -p {config[output_dir]}/base_recal; "
        "gatk BaseRecalibrator -R {input.ref} -I {input.bam} -known-sites {input.known_sites} -O {output}"

rule apply_base_recal:
    """Apply base quality score recalibration"""
    input:
        ref = config['ref_genome']
        bam = rules.rm_dups.output.final_bam
        recal_table = rules.base_recal.output
    output:
        config['output_dir'] + "/base_recal/{sample}.recal.final.bam"
    shell:
        "gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output}"


rule haplotype:
    """Make a file with annotated variants"""
    input:
        ref = config['ref_genome']
        bam = rules.apply_base_recal.output
    output:
        config['output_dir'] + "/haplotype/{sample}.snps.g.vcf.gz"
    shell:
        "gatk HaplotypeCaller "
        "-R {input.ref} -I {input.bam} -O {output} -ERC GVCF"
        "-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation"

# NOTE: combine and genotype are untested!

rule combine:
    """Combine the g.vcf files"""
    input:
        ref = config['ref_genome']
        vcf = ['-V '+file for file in glob.glob(rules.haplotype.output)]
    output:
        config['output_dir'] + "haplotype/ALL.g.vcf.gz"
    shell:
        "gatk CombineGVCFs -R {input.ref} -O {output} "
        "-G StandardAnnotation -G AS_StandardAnnotation "
        "{input.vcfs}"

rule genotype:
    """Perform joint genotyping on all of the samples"""
    input:
        ref = config['ref_genome']
        vcf = rules.combine.output
    output:
        config['output_dir'] + "/haplotype/ALL.genotype.vcf.gz"
    shell:
        "gatk VariantRecalibrator -R {input.ref} -V {input.vcf} "
        "-G StandardAnnotation -G AS_StandardAnnotation"

rule variant_filter:
    """Filter variants by QD, FS, MQ, MQRankSum, ReadPosRankSum, and SOR"""
    input:
        ref = config['ref_genome']
        vcf = rules.genotype.output
        hapmap = config['hapmap']
        omni = config['omni']
        1000G = config['1000G']
        dbsnp = config['dbSNP']
    output:
        recal = config['output_dir'] + "/variant_filter/snps.recal"
        tranches = config['output_dir'] + "/variant_filter/snps.tranches"
    shell:
        "mkdir -p {config[output_dir]}/variant_filter; "
        "gatk VariantRecalibrator -R {input.ref} -V {input.vcf} "
        "--resource hapmap,known=false,training=true,truth=true,prior=15.0:{input.hapmap} "
        "--resource omni,known=false,training=true,truth=false,prior=12.0:{input.omni} "
        "--resource 1000G,known=false,training=true,truth=false,prior=10.0:{input.1000G} "
        "--resource dbsnp,known=true,training=false,truth=false,prior=2.0:{input.dbsnp} "
        "-an QD -an FS -an MQ -an MQRankSum -an ReadPosRankSum -an SOR"
        "-mode SNP -O {output.recal} --tranches-file {output.tranches}"

rule apply_variant_filter:
    """Create a file with only variants that have passed filtration"""


rule vcf2h5:
    """Convert VCF data files to HDF5 format"""
    input:
        chrom = config['chrom_info'],
        vcfs = glob.glob(config['vcf_dir'] + "/*chr*.vcf.gz")
    output:
        snp_index = config['snp_h5_dir'] + "/snp_index.h5",
        snp_tab = config['snp_h5_dir'] + "/snp_tab.h5",
        haplotype = config['snp_h5_dir'] + "/haplotype.h5"
    shell:
        "mkdir -p {config[snp_h5_dir]}; "
        "{config[wasp_dir]}/snp2h5/snp2h5 "
        "  --chrom {input.chrom} "
        "  --format vcf "
        "  --snp_index {output.snp_index} "
        "  --snp_tab {output.snp_tab} "
        "  --haplotype {output.haplotype} "
        "  {input.vcfs}"


rule find_intersecting_snps_paired_end:
    """find intersecting SNPs using WASP script"""
    input:
        bam = config["output_dir"] + "/map1_sort/{sample}.bam",
        snp_index = config["snp_h5_dir"] + "/snp_index.h5",
        snp_tab = config["snp_h5_dir"] + "/snp_tab.h5",
        haplotype = config['snp_h5_dir'] + "/haplotype.h5"
    output:
        fastq1 = config["output_dir"] + "/find_intersecting_snps/{sample}.remap.fq1.gz",
        fastq2 = config["output_dir"] + "/find_intersecting_snps/{sample}.remap.fq2.gz",
        keep_bam = config["output_dir"] + "/find_intersecting_snps/{sample}.keep.bam",
        remap_bam = config["output_dir"] + "/find_intersecting_snps/{sample}.to.remap.bam"
    shell:
        "mkdir -p {config[output_dir]}/find_intersecting_snps ; "
        "{config[py2]} {config[wasp_dir]}/mapping/find_intersecting_snps.py "
        "    --is_paired_end "
        "    --is_sorted "
        "    --output_dir {config[output_dir]}/find_intersecting_snps "
        "    --snp_tab {input.snp_tab} "
        "    --snp_index {input.snp_index} "
        "    --haplotype {input.haplotype} "
        "    --samples {config[sample_file]} "
        "    {input.bam}"


rule map_bowtie2_paired_end1:
    """map reads using bowtie2"""
    input:
        fastq1 = lambda wildcards: read_samples()[wildcards.sample][1][0],
        fastq2 = lambda wildcards: read_samples()[wildcards.sample][1][1]
    output:
        config["output_dir"] + "/map1/{sample}.bam"
    shell:
        "mkdir -p " + config["output_dir"] + "/map1 ; "
        "{config[bowtie2]} -x {config[bowtie2_index]} -1 {input.fastq1} -2 {input.fastq2} -p 12"
        "| {config[samtools]} view -b -q 10 - > {output} "


rule sort_and_index_bam1:
    """sort and index bam generated by first mapping step"""
    input:
        config["output_dir"] + "/map1/{sample}.bam"
    output:
        config["output_dir"] + "/map1_sort/{sample}.bam",
        config["output_dir"] + "/map1_sort/{sample}.bam.bai"
    shell:
        "mkdir -p {config[output_dir]}/map1_sort ; "
        "{config[samtools]} sort -o {output[0]} {input}; "
        "{config[samtools]} index {output[0]}"


rule map_bowtie2_paired_end2:
    """map reads a second time using bowtie2"""
    input:
        fastq1 = config['output_dir'] + "/find_intersecting_snps/{sample}.remap.fq1.gz",
        fastq2 = config['output_dir'] + "/find_intersecting_snps/{sample}.remap.fq2.gz"
    output:
        config["output_dir"] + "/map2/{sample}.bam"
    shell:
        "mkdir -p " + config["output_dir"] + "/map2 ; "
        "{config[bowtie2]} -x {config[bowtie2_index]} -1 {input.fastq1} -2 {input.fastq2} -p 12"
        "| {config[samtools]} view -b -q 10 - > {output}"


rule sort_and_index_bam2:
    """sort and index bam generated by second mapping step"""
    input:
        config["output_dir"] + "/map2/{sample}.bam"
    output:
        config["output_dir"] + "/map2_sort/{sample}.bam",
        config["output_dir"] + "/map2_sort/{sample}.bam.bai"
    shell:
        "mkdir -p {config[output_dir]}/map2_sort ; "
        "{config[samtools]} sort -o {output[0]} {input} ; "
        "{config[samtools]} index {output[0]}"


rule filter_remapped_reads:
    """filter reads from second mapping step"""
    input:
        to_remap_bam = config['output_dir'] + "/find_intersecting_snps/{sample}.to.remap.bam",
        remap_bam = config['output_dir'] + "/map2_sort/{sample}.bam",
    output:
        keep_bam = config['output_dir'] + "/filter_remapped_reads/{sample}.keep.bam"
    shell:
        "mkdir -p {config[output_dir]}/filter_remapped_reads ; "
        "{config[py2]} {config[wasp_dir]}/mapping/filter_remapped_reads.py "
        "  {input.to_remap_bam} {input.remap_bam} {output.keep_bam}"


rule merge_bams:
    """merge 'keep' BAM files from mapping steps 1 and 2, then sort and index"""
    input:
        keep1 = config['output_dir'] + "/find_intersecting_snps/{sample}.keep.bam",
        keep2 = config['output_dir'] + "/filter_remapped_reads/{sample}.keep.bam"
    output:
        merge = config['output_dir'] + "/merge/{sample}.keep.merge.bam",
        sort = config['output_dir'] + "/merge/{sample}.keep.merge.sort.bam"
    shell:
        "mkdir -p {config[output_dir]}/merge ; "
        "{config[samtools]} merge {output.merge} {input.keep1} {input.keep2}; "
        "{config[samtools]} sort -o {output.sort} {output.merge}; "
        "{config[samtools]} index {output.sort}"


rule rmdup_pe:
    """remove duplicate read pairs"""
    input:
        config['output_dir'] + "/merge/{sample}.keep.merge.sort.bam"
    output:
        rmdup = config['output_dir'] + "/rmdup/{sample}.keep.merge.rmdup.bam",
        sort = config['output_dir'] + "/rmdup/{sample}.keep.merge.rmdup.sort.bam"
    shell:
        "mkdir -p {config[output_dir]}/rmdup ; "
        "{config[py2]} {config[wasp_dir]}/mapping/rmdup_pe.py {input} {output.rmdup} ;"
        "{config[samtools]} sort -o {output.sort} {output.rmdup}; "
        "{config[samtools]} index {output.sort}"


rule get_as_counts:
    """get allele-specific read counts for SNPs"""
    input:
        bam = config['output_dir'] + "/rmdup/{sample}.keep.merge.rmdup.sort.bam",
        snp_index = config["snp_h5_dir"] + "/snp_index.h5",
        snp_tab = config["snp_h5_dir"] + "/snp_tab.h5",
        haplotype = config['snp_h5_dir'] + "/haplotype.h5",
    params:
        samp1kg = lambda wildcards: SAMP_TO_1KG[wildcards.sample]
    output:
        config['output_dir'] + "/as_counts/{sample}.as_counts.txt.gz"
    shell:
        "mkdir -p {config[output_dir]}/as_counts ; "
        "{config[py2]} {config[wasp_dir]}/mapping/get_as_counts.py "
        "  --snp_tab {input.snp_tab} "
        "  --snp_index {input.snp_index} "
        "  --haplotype {input.haplotype} "
        "  --samples {config[sample_file]} "
        "  --genotype_sample {params.samp1kg} "
        "  {input.bam} | gzip > {output}"
