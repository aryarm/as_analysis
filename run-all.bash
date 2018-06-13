#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd

out_path="out"
mkdir -p $out_path

snakemake -s Snakefile-variant_calling --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qout" -j 24 --config output_dir=${out_path} --latency-wait 60 --use-conda >>${out_path}/out 2>&1

snakemake -s Snakefile-WASP --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qout" -j 24 --config output_dir=${out_path} --latency-wait 60 --use-conda >>${out_path}/out 2>&1

snakemake -s Snakefile-counts --cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qout" -j 24 --config output_dir=${out_path} --use-conda >>${out_path}/out 2>&1
