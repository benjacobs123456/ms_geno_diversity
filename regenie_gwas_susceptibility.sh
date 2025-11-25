#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(CSA AFR EUR)
ancestries_2=(sas afr eur)
ancestry=${ancestries[$index]}
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry


# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# run in PLINK
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--extract /data/scratch/hmy117/snps_to_keep_for_gwas$ancestry_out\.tsv \
--covar ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--pheno ./pheno/susceptibility_$ancestry_out\_pheno.tsv \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--covar-name age,sex,PC1,PC2,PC3,PC4 \
--glm hide-covar \
--mac 20 \
--out ./outputs/plink_risk_gwas_$ancestry_out


# step 1 REGENIE
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--extract /data/scratch/hmy117/snps_to_keep_for_gwas$ancestry_out\.tsv \
--bed /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--covarFile ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_$ancestry_out\_pheno.tsv \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--bsize 1000 \
--bt --lowmem --firth approx \
--force-step1 \
--pThresh 0.001 \
--out ./outputs/case_control_step1_$ancestry_out\_fit \
--catCovarList sex \
--covarColList age,PC{1:4}

# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--extract /data/scratch/hmy117/snps_to_keep_for_gwas$ancestry_out\.tsv \
--bed /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--covarFile ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_$ancestry_out\_pheno.tsv \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--bsize 1000 \
--bt --lowmem --firth approx \
--pThresh 0.001 \
--pred ./outputs/case_control_step1_$ancestry_out\_fit_pred.list \
--catCovarList sex \
--covarColList age,PC{1:4} \
--minMAC 20 \
--af-cc \
--out ./outputs/case_control_gwas_$ancestry_out
