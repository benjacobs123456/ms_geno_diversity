#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_sev
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr eur)
ancestry=${ancestries[$index]}

echo doing $ancestry


# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/Cambridge/

# step 1 REGENIE
	~/regenie \
--step 1 \
--extract /data/scratch/hmy117/adams_imputed/snps_to_keep_step1_regenie.tsv \
--bed  /data/scratch/hmy117/cambridge_genos/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile $ancestry\_covars_with_pcs.tsv \
--phenoFile $ancestry\_pheno.tsv \
--bsize 1000 \
--qt --lowmem --apply-rint \
--out step1_$ancestry\_fit \
--catCovarList Sex,study \
--covarColList age_at_recruitment,PC{1:4}


# step 2 REGENIE
	~/regenie \
--step 2 \
--bed  /data/scratch/hmy117/cambridge_genos/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile $ancestry\_covars_with_pcs.tsv \
--phenoFile $ancestry\_pheno.tsv \
--bsize 1000 \
--qt --lowmem --apply-rint \
--pred step1_$ancestry\_fit_pred.list \
--out sev_gwas_$ancestry \
--catCovarList Sex,study \
--covarColList age_at_recruitment,PC{1:4}
