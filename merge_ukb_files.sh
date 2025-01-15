#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N merge_ukb
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/


# split adams chroms & filter to UKB variants
~/plink2 --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed \
--chr ${SGE_TASK_ID} \
--extract /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID}\.bim \
--make-bed \
--out /data/scratch/hmy117/adams_for_merge_chr${SGE_TASK_ID} \
--rm-dup exclude-all 

~/plink2 --bfile /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID} \
--chr ${SGE_TASK_ID} \
--extract /data/scratch/hmy117/adams_for_merge_chr${SGE_TASK_ID}\.bim \
--make-bed \
--out /data/scratch/hmy117/ukb_for_merge_chr${SGE_TASK_ID} \
--rm-dup exclude-all 


~/plink --bfile /data/scratch/hmy117/adams_for_merge_chr${SGE_TASK_ID} \
--bmerge /data/scratch/hmy117/ukb_for_merge_chr${SGE_TASK_ID} \
--out /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}


rm /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID}\.bed;
rm /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID}\.bim;
rm /data/scratch/hmy117/post_qc_ukb_chr${SGE_TASK_ID}\.fam;

rm /data/scratch/hmy117/ukb_for_merge_chr${SGE_TASK_ID}\.bed; 
rm /data/scratch/hmy117/ukb_for_merge_chr${SGE_TASK_ID}\.bim;
rm /data/scratch/hmy117/ukb_for_merge_chr${SGE_TASK_ID}\.fam;
