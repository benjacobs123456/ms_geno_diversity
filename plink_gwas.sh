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


# glm
~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--pheno ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--covar ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--glm hide-covar \
--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm$ancestry_out \
--threads $NSLOTS

# glm no covar
~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--pheno ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--glm hide-covar allow-no-covars \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm_nocov$ancestry_out \
--threads $NSLOTS

# glm with age
~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--pheno ./pheno/reimputed_$ancestry_out\_pheno.tsv \
--covar ./pheno/reimputed_$ancestry_out\_covars_with_pcs.tsv \
--glm hide-covar \
--covar-name age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm_age$ancestry_out \
--threads $NSLOTS

