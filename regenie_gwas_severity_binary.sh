#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_sev_bt
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr eur)
ancestry=${ancestries[$index]}

echo doing $ancestry


# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# step 1 REGENIE
	~/regenie \
--step 1 \
--extract /data/scratch/hmy117/adams_imputed/snps_to_keep_step1_regenie.tsv \
--bed ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile ./pheno/$ancestry\_covars_with_pcs.tsv \
--phenoFile ./pheno/$ancestry\_pheno.tsv \
--bsize 1000 \
--bt --lowmem --firth \
--out bt_step1_$ancestry\_fit \
--catCovarList Sex,Site \
--covarColList ageatedss,agesq,agesex,PC{1:4} \
--phenoCol edss6


# step 2 REGENIE
	~/regenie \
--step 2 \
--bed ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile ./pheno/$ancestry\_covars_with_pcs.tsv \
--phenoFile ./pheno/$ancestry\_pheno.tsv \
--bsize 1000 \
--bt --lowmem --firth \
--pred bt_step1_$ancestry\_fit_pred.list \
--catCovarList Sex,Site \
--covarColList ageatedss,agesq,agesex,PC{1:4} \
--out ./outputs/sev_gwas_$ancestry \
--phenoCol edss6
