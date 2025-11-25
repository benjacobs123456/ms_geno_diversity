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
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# step 2 REGENIE
	~/regenie \
--step 2 \
--bed ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile ./pheno/$ancestry\_covars_with_pcs.tsv \
--phenoFile ./pheno/$ancestry\_pheno.tsv \
--bsize 1000 \
--qt --lowmem --apply-rint \
--pred step1_$ancestry\_fit_pred.list \
--catCovarList Sex \
--out ./outputs/sev_recessive_gwas_$ancestry \
--catCovarList Sex,Site \
--covarColList ageatedss,agesq,agesex,PC{1:4} \
--phenoExcludeList edss6 \
--test recessive

