#!/bin/bash
#BSUB -W 120:00
#BSUB -o /path/to/log.%I.%J.out
#BSUB -e /path/to/log.%I.%J.out
#BSUB -q cgel
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -J "MK545-A"

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
# 1) Create a pair manifest (tab-delimited)
#printf "pair_id\ttumor_id\ttumor_bam\tnormal_id\tnormal_bam\tgender\nMK545-A_pair\tMK545-A\t/rsrch6/home/genetics/zzhang18/projects/WGS_pipeline/scripts/pipelines/tests/data/MK545-A.bam\tMK545-Control\t/rsrch6/home/genetics/zzhang18/projects/WGS_pipeline/scripts/pipelines/tests/data/MK545-Control.bam\tmale\n" \
#  > /rsrch6/home/genetics/zzhang18/projects/WGS_pipeline/scripts/pipelines/tests/samplesheet.tsv
ml nextflow
Dir="/dir/to/pipeline"
# 2) Run locally

cd ${Dir}/tests

nextflow run ${Dir}/main.nf \
  -resume \
  -c ${Dir}/conf/nextflow.config \
  --pairs_tsv /path/to/samplesheet.tsv \
  --outdir /path/to/nf_out \
  --ref hg38 \
  --pipeline ${Dir} \
  --refdir /Dir/to/refs \
  --threads 10 
