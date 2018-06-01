#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd
#$ -o /iblm/netapp/home/amassarat/allele_specific_analysis/final/out

main="/iblm/netapp/home/amassarat/allele_specific_analysis"
wasp_dir="/iblm/netapp/home/amassarat/bin/WASP"
reads="${main}/test_data"
ref_data="${main}/ref_data"
other_data="${main}/other_data"
example_data="${main}/example_data"
main="${main}/final"
num_threads=20

# make directory for BWA and samtools
mkdir -p $main/dna_align
cd $main/dna_align

# align using BWA-MEM
printf "________________________ aligning using BWA-MEM ________________________\n"
# note that we use the "-M" flag for picard compatibility
# we use -R to specify read group info so that haplotype caller doesn't scream at us
# and we use -t to specify the number of threads to run on
bwa mem -M -R "@RG\tID:group1\tSM:Jurkat\tPL:ILLUMINA\tLB:lib1\tPU:unit1" -t $num_threads $ref_data/BWAIndex/genome.fa $reads/DNA_R1.fq.gz $reads/DNA_R2.fq.gz |
	{
		# convert sam to bam (so that we use less space on the cluster) as specified by the -b flag
		# we use -u to create an uncompressed bam file
		# we use -q to filter alignments with MAPQ scores less than 20
		# I'm not sure what the -F option was for. you should probably figure this out
		samtools view -u -b -F 4 -q 20 -@ $num_threads
	} | {
		# and sort the bam output by name
		# this is needed for fixmate to run properly downstream of us
		samtools sort -n -@ $num_threads > $main/dna_align/aln.bam
	}

# use fixmate to fill in mate coordinates and mate related flags, since our data is pair-ended.
# we need the MC tags (included because we used the -m flag) that it creates for markdup
printf "\n________________________ filling in mate coordinates and mate related flags ________________________\n"
samtools fixmate -m $main/dna_align/aln.bam $main/dna_align/aln.mated.bam

# sort the file (not by name this time)
# this is needed for markdup to run properly downstream of us
printf "\n________________________ sorting by coordinates ________________________\n"
samtools sort -@ $num_threads -o $main/dna_align/aln.mated.sorted.bam $main/dna_align/aln.mated.bam

# remove duplicates that occurred from PCR
printf "\n________________________ removing duplicates ________________________\n"
samtools markdup -@ $num_threads $main/dna_align/aln.mated.sorted.bam $main/dna_align/aln.rmdup.bam

# index the bam output so its easier/faster to use
printf "\n________________________ indexing ________________________\n"
samtools index -b -@ $num_threads $main/dna_align/aln.rmdup.bam

# did everything work?
if [ "$?" -ne "0" ]; then
	echo "samtools failed" 1>&\n2
	exit 1
fi

# make a directory for base recalibration
mkdir -p $main/base_recal
cd $main/base_recal

# use gatk's base quality score recalibration to make the base quality scores less biased
printf "\n________________________ recalibrating base quality scores ________________________\n"
gatk BaseRecalibrator -R $ref_data/WholeGenomeFasta/genome.fa -I $main/dna_align/aln.rmdup.bam -known-sites $other_data/dbsnp.vcf -O $main/base_recal/recal_data.table
# use the quality report file (recal_data.table) to create a bam file with the corrected quality scores
gatk ApplyBQSR -R $ref_data/WholeGenomeFasta/genome.fa -I $main/dna_align/aln.rmdup.bam --bqsr-recal-file $main/base_recal/recal_data.table -O $main/base_recal/aln_final.bam

# make a directory for haplotype calling
mkdir -p $main/haplotype
cd $main/haplotype

# make VCF using GATK
printf "\n________________________ creating annotated VCF file ________________________\n"
# note that we use the -G options to tell gatk to include allele-specific versions of the standard annotations
# also note that we must use GVCF mode if we want it to output annotations
gatk HaplotypeCaller -R $ref_data/WholeGenomeFasta/genome.fa -I $main/base_recal/aln_final.bam -O $main/haplotype/snps.noannotate.vcf.gz
#	-ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
# we also want to annotate our VCFs. note that this step doesn't need to happen when we have multiple samples because haplotypecaller can do it on the fly
gatk VariantAnnotator -R $ref_data/WholeGenomeFasta/genome.fa -V $main/haplotype/snps.noannotate.vcf.gz -O $main/haplotype/snps.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation

# FOR MANY SAMPLES:
# combine and genotype the g.vcf files...
# printf "\n________________________ combining and genotyping the g.vcf files ________________________\n"
# gatk CombineGVCFs -R $ref_data/WholeGenomeFasta/genome.fa -V $main/haplotype/snps.vcf.gz -O $main/haplotype/output-genotyped.vcf -G StandardAnnotation -G AS_StandardAnnotation
# gatk GenotypeGVCFs -R $ref_data/WholeGenomeFasta/genome.fa -V $main/haplotype/snps.vcf.gz -O $main/haplotype/output-genotyped.vcf -G StandardAnnotation -G AS_StandardAnnotation

# make a directory for filtering variants
mkdir -p $main/variant_filter
cd $main/variant_filter

# Filter variants by a couple different parameters (QD, FS, MQ, MQRankSum, ReadPosRankSum, and SOR - as described below).
# Note that we use the VariantRecalibrator to decide these parameters based on some pre-determined data (using machine learning).
# QD - QualByDepth (the variant confidence divided by the unfiltered depth of non-reference samples) <-- I used a more stringent QD because apparenly depth is important
# FS - FisherStrand (the Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads)
# MQ - RMSMappingQuality (the Root Mean Square of the mapping quality of the reads across all samples)
# MQRankSum - MappingQualityRankSumTest (the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities; this is only applied to heterozygous calls)
# ReadPosRankSum - ReadPosRankSumTest (the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele; this is only applied to heterozygous calls)
# SOR - StrandOddsRatio (aims to evaluate whether there is strand bias in the data)
printf "\n________________________ filtering variants ________________________\n"
gatk VariantRecalibrator -R $ref_data/WholeGenomeFasta/genome.fa -V $main/haplotype/snps.vcf.gz \
	--resource hapmap,known=false,training=true,truth=true,prior=15.0:$other_data/hapmap.vcf \
	--resource omni,known=false,training=true,truth=false,prior=12.0:$other_data/omni.vcf \
	--resource 1000G,known=false,training=true,truth=false,prior=10.0:$other_data/1000G.vcf \
	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$other_data/dbsnp.vcf \
	-an QD -an FS -an MQ -an MQRankSum -an ReadPosRankSum -an SOR \
	-mode SNP \
	-O $main/variant_filter/snps.recal --tranches-file $main/variant_filter/snps.tranches
gatk ApplyVQSR -R $ref_data/WholeGenomeFasta/genome.fa -V $main/haplotype/snps.vcf.gz -mode SNP --ts-filter-level 97.5 --recal-file $main/variant_filter/snps.recal --tranches-file $main/variant_filter/snps.tranches -O $main/variant_filter/snps.filter.vcf.gz

# # filter variants by a couple parameters (almost the same as the defaults listed on the GATK website)
# printf "\n________________________ filtering variants ________________________\n"
# gatk VariantFiltration -R $ref_data/WholeGenomeFasta/genome.fa -V $main/snps.vcf.gz -O $main/snps.filter.vcf.gz \
# 	-filter "QD < 12.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_filter"
# 	# -G-filter "isHet == 1" -G-filter-name "hets"

# first, make a directory to store the final VCF data
mkdir -p $main/genotypes
cd $main/genotypes

printf "\n________________________ subsetting by heterozygotes and filters ________________________\n"
zcat $main/variant_filter/snps.filter.vcf.gz | SnpSift filter "(countHet() > 0) && (FILTER == 'PASS')" > $main/genotypes/ALL

# split the VCF file that we got by chromosome, since that's what WASP expects
printf "\n________________________ creating chromosome VCF files ________________________\n"
# now break down the vcf file, and gzip all of the vcf files that are created
SnpSift split ALL
gzip *.vcf
# switch back to the main directory
cd $main

# TODO: consider using STAR instead of manually calling the WASP pipeline?
# STAR --runThreadN $num_threads --genomeDir $ref_data/StarIndex --readFilesIn $reads/RNA_R1.fq.gz $reads/RNA_R2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $main/ --waspOutputMode SAMtag --varVCFfile $main/dna-snps.vcf

# switch into the right conda environments to run WASP
source activate py35
source activate wasp

# use snakemake to carry out WASP
printf "\n________________________ initiating WASP pipeline ________________________\n"
snakemake -s Snakefile --cluster "qsub -t 1 -V -q iblm.q -j y -o $main/out" -j 24 --rerun-incomplete --configfile $main/snake_conf.yaml

# Get rid of many of the files that WASP created.
# but first switch into the directory where WASP output is defined
cd $main/WASP
# remove unsorted files from rmdup dir:
ls $PWD/rmdup/* | grep -v sort | xargs rm
# remove intermediate files and directories:
rm -r map1 map1_sort find_intersecting_snps map2 map2_sort filter_remapped_reads merge
# switch back to the previous directory
cd -

# switch back out of the WASP-specific conda environments
source deactivate wasp
source deactivate py35

# lastly, create a table with GQ scores from the filtered VCF for later analysis
# the table will be tab separated and can be easily imported to R
gatk VariantsToTable -R $ref_data/WholeGenomeFasta/genome.fa -V $main/variant_filter/snps.filter.vcf.gz -O gq_subset.table \
	-F CHROM -F POS -F REF -F ALT -GF GQ --moltenize

printf "\n________________________ DONE! ________________________\n"
echo "The final file(s) with allele specific counts is in ${main}/WASP/as_counts"
