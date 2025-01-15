#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
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

# filter to common vars


# get step 1 snps
~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--exclude mhc_snps_to_exclude_pca_$ancestry_out \
--extract /data/scratch/hmy117/info_stats_for_step1_$ancestry_out\.tsv \
--keep ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--out /data/scratch/hmy117/reimputed_step1_snps_$ancestry_out \
--threads $NSLOTS \
--make-bed

~/plink --bfile /data/scratch/hmy117/reimputed_step1_snps_$ancestry_out \
--indep-pairwise 1000 100 0.9 \
--out /data/scratch/hmy117/reimputed_step1_snps_pruned_$ancestry_out \
--threads $NSLOTS


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--extract /data/scratch/hmy117/reimputed_step1_snps_pruned_$ancestry_out\.prune.in \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--force-step1 \
--out /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit \
--catCovarList sex \
--covarColList PC{1:10} \
--loocv


# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--loocv \
--bt --lowmem \
--pred /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_pred.list \
--catCovarList sex \
--covarColList PC{1:10} \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_$ancestry_out

################################
# Sensitivity analyses

#################################
# Age 
#################################

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--extract /data/scratch/hmy117/reimputed_step1_snps_pruned_$ancestry_out\.prune.in \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--force-step1 \
--out /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_age \
--catCovarList sex \
--covarColList age,PC{1:10} \
--loocv


# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--loocv \
--bt --lowmem \
--pred /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_age_pred.list \
--catCovarList sex \
--covarColList age,PC{1:10} \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_age_$ancestry_out

#################################
# PC 1 - 2 
#################################


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--extract /data/scratch/hmy117/reimputed_step1_snps_pruned_$ancestry_out\.prune.in \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--force-step1 \
--out /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_pcs1_2 \
--catCovarList sex \
--covarColList age,PC{1:2} \
--loocv


# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--loocv \
--bt --lowmem \
--pred /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_pcs1_2_pred.list \
--catCovarList sex \
--covarColList age,PC{1:2} \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_pcs_1_2_$ancestry_out

#################################
# Use all SNPs for step 1 SNPs 
#################################

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--force-step1 \
--out /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_allsnps \
--catCovarList sex \
--covarColList PC{1:10} \
--loocv


# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--keep ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--bed /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--covarFile ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--loocv \
--bt --lowmem \
--pred /data/scratch/hmy117/case_control_step1_ALL_$ancestry_out\_fit_allsnps_pred.list \
--catCovarList sex \
--covarColList PC{1:10} \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_allsnps_step1$ancestry_out
