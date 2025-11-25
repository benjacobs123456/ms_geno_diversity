#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N bcf_filter
#$ -o /data/scratch/hmy117
#$ -t 1:22

# filter VCFs
cd /data/scratch/hmy117/adams_imputed/

module load bcftools
bcftools filter -i "R2>0.7 & MAF>0.01" chr${SGE_TASK_ID}\.dose.vcf.gz -Oz -o chr${SGE_TASK_ID}\_filtered.vcf.gz
