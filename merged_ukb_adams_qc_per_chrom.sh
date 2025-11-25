#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N filter_out_palindromes_and_do_qc_ukb
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/

# exclude palindromes

awk '{if( ($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="G" && $6=="C") || ($5=="C" && $6=="G")) print $2}' /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\.bim > /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_palindromic_snps_to_exclude

~/plink --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID} \
--maf 0.01 \
--exclude /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_palindromic_snps_to_exclude \
--hwe 1e-10 \
--geno 0.1 \
--allow-no-sex \
--out /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_no_palindromes \
--make-bed
