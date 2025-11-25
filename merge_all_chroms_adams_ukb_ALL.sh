#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=25G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N merge_ukb
#$ -o /data/scratch/hmy117


cd /data/scratch/hmy117/
rm filelist_for_merge
for i in {2..22}; do echo /data/scratch/hmy117/ukb_adams_merged_genotypes_chr$i >> filelist_for_merge; done
~/plink --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_chr1 \
--merge-list filelist_for_merge \
--make-bed \
--out /data/scratch/hmy117/ukb_adams_merged_genotypes_all_chroms_all

