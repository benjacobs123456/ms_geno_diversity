#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N make_kg_reference_for_imputation
#$ -o /data/scratch/hmy117
#$ -t 1:22



cd /data/Wolfson-UKBB-Dobson/1kg_reference/raw_files
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# make reference 
~/Minimac4/bin/minimac4 \
--compress-reference ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
--output /data/Wolfson-UKBB-Dobson/1kg_reference/msav_files/compressed_chr${SGE_TASK_ID}\.msav 

