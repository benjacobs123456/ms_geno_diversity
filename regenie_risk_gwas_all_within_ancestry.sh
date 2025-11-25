#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_all_ancestry
#$ -o /data/scratch/hmy117
#$ -t 1:2


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(CSA AFR EUR)
ancestries_2=(sas afr eur)
ancestry=${ancestries[$index]}
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# fisher's exact test 
~/plink --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--pheno ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--assoc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_fisher$ancestry_out \
--allow-no-sex 

# step 1 REGENIE
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--keep ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--covarFile ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--firth --approx \
--force-step1 \
--out ./outputs/case_control_step1_ALL_$ancestry_out\_fit \
--catCovarList sex \
--niter 1000 \
--covarColList age,PC{1:10} \
--loocv


# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--bed /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--keep ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--covarFile ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--loocv \
--bt --lowmem \
--pred ./outputs/case_control_step1_ALL_$ancestry_out\_fit_pred.list \
--catCovarList sex \
--covarColList age,PC{1:10} \
--minMAC 20 \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_$ancestry_out


# PLINK glm  
~/plink2 --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--pheno ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--covar ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--glm hide-covar \
--covar-variance-standardize \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_plink_glm_$ancestry_out

