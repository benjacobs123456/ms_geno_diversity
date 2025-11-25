#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_qc_ukb
#$ -o /data/scratch/hmy117
#$ -t 1:22

# set wd
cd /data/scratch/hmy117

# Filter UKB SNPs to ADAMS SNPs

~/plink2 \
--pfile /data/Wolfson-PNU-dementia/UKB/imputed_genotypes/plink2_files/chr_${SGE_TASK_ID} \
--extract /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_keep.tsv \
--make-bed \
--rm-dup exclude-all \
--out /data/scratch/hmy117/ukb_chr${SGE_TASK_ID}\_filtered_tmp
