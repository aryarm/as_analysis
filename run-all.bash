#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd

# An example bash script demonstrating how to run the entire snakemake pipeline
# on an SGE cluster
# This script creates two separate log files:
# 	1) log - the basic snakemake log of completed rules
# 	2) qlog - a more detailed log of the progress of each rule and any errors

# you can specify a directory for all output here:
out_path="out"
mkdir -p "$out_path"

# Before running this snakemake pipeline, remember to verify that the config.yaml
# file has been appropriately completed with the required input info. In
# particular, make sure that you have created a samples.tsv file specifying
# paths to the fastq files for each of your samples.

snakemake \
--cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
-j 24 \
--config output_dir="${out_path}" \
--latency-wait 60 \
--use-conda \
-k \
--restart-times 2 \
&>"${out_path}/log"







# Individual portions of the snakemake pipeline can also be run on their own.
# See the "Snakefiles/" directory for Snakefiles of each of the following
# three portions of the pipeline (and the "configs/" directory for their
# corresponding config files):
# 
# 
# 1) Snakefile-variant_calling - align DNA fastq's and generate a filtered
#    VCF containing heterozygous SNPs
# 
# snakemake \
# -s Snakefiles/Snakefile-variant_calling \
# --configfile configs/config-variant_calling.yaml \
# --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
# -j 24 \
# --config output_dir=${out_path} \
# --latency-wait 60 \
# --use-conda \
# -k \
# --restart-times 2 \
# >>${out_path}/log 2>&1
# 
# 
# 2) Snakefile-WASP - align RNA fastq's and filter using WASP to reduce
#    mapping bias
# 
# snakemake \
# -s Snakefiles/Snakefile-WASP \
# --configfile configs/config-WASP.yaml \
# --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
# -j 24 \
# --config output_dir=${out_path} \
# --latency-wait 60 \
# --use-conda \
# -k \
# --restart-times 2 \
# >>${out_path}/log 2>&1
# 
# 
# 3) Snakefile-counts - retrieve counts of reads overlapping SNPs from BAM
#    files generated in each of the previous portions of the pipeline and use
#    them to find genes which demonstrate allelic imbalance
# 
# snakemake \
# -s Snakefiles/Snakefile-counts \
# --configfile configs/config-counts.yaml \
# --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
# -j 24 \
# --config output_dir=${out_path} \
# --latency-wait 60 \
# --use-conda \
# --restart-times 2 \
# -k \
# >>${out_path}/log 2>&1
