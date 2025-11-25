#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_qc_ukb2
#$ -o /data/scratch/hmy117
#$ -t 1:22

# set wd
cd /data/scratch/hmy117

# Filter UKB SNPs to ADAMS SNPs


~/plink2 \
--bfile  /data/scratch/hmy117/ukb_chr${SGE_TASK_ID}\_filtered_tmp \
--update-name /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_update_names.tsv \
--out /data/scratch/hmy117/filtered_ukb_chr${SGE_TASK_ID}\_hg19_cpra \
--make-bed

wait

echo "All done"
