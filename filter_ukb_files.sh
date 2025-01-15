#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N filter_ukb
#$ -o /data/scratch/hmy117
#$ -t 1:22

~/plink2 --pfile /data/Wolfson-PNU-dementia/UKB/imputed_genotypes/plink2_files/chr_${SGE_TASK_ID} \
--extract /data/scratch/hmy117/ukb_snps_to_keep.tsv \
--chr ${SGE_TASK_ID} \
--make-bed \
--rm-dup exclude-all \
--out /data/scratch/hmy117/post_qc_ukb_tmp_chr${SGE_TASK_ID}

~/plink2 --bfile /data/scratch/hmy117/post_qc_ukb_tmp_chr${SGE_TASK_ID} \
--update-name /data/scratch/hmy117/ukb_snps_to_update.tsv \
--out /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID} \
--make-bed 

rm /data/scratch/hmy117/post_qc_ukb_tmp_chr${SGE_TASK_ID}\.bed
rm /data/scratch/hmy117/post_qc_ukb_tmp_chr${SGE_TASK_ID}\.bim
rm /data/scratch/hmy117/post_qc_ukb_tmp_chr${SGE_TASK_ID}\.fam
