#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_qc
#$ -o /data/scratch/hmy117
#$ -t 1:22

# set wd
cd ~/ADAMS/


# filter by sample
~/plink2 --pfile /data/Wolfson-PNU-dementia/UKB/imputed_genotypes/plink2_files/chr_${SGE_TASK_ID} \
--keep ~/ADAMS/genotypes/QMUL_Aug_23/outputs/sas_ukb_ids.tsv \
--out /data/scratch/hmy117/ukb_geno_sas_chr${SGE_TASK_ID} \
--make-pgen \
--set-all-var-ids @:#:\$1:\$2 \
--snps-only

wait

# filter
~/plink2 --pfile /data/scratch/hmy117/ukb_geno_sas_chr${SGE_TASK_ID} \
--maf 0.01 \
--hwe 1e-5 \
--geno 0.1 \
--out /data/scratch/hmy117/qcd_ukb_geno_sas_chr${SGE_TASK_ID} \
--make-bed
