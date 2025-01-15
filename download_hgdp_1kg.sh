#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N download_hgdp_genomes
#$ -o /data/scratch/hmy117
#$ -t 1:22

~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${SGE_TASK_ID}\.vcf.bgz \
/data/scratch/hmy117/hgdp_1kg_genomes

