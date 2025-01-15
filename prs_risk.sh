#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N prs_risk
#$ -o /data/scratch/hmy117
#$ -t 1:2


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr)
ancestry=${ancestries[$index]}

echo doing $ancestry


# run PRS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# filter in PLINK
~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry\_ALL_CHRS \
--keep ./pheno/reimputed_$ancestry\_covars_with_pcs.tsv \
--maf 0.05 \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_prs_tmp_$ancestry \
--make-bed

~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_prs_tmp_$ancestry \
--set-all-var-ids chr@:# \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_prs_$ancestry \
--make-bed


# NO MHC
for r2 in 0.001 0.01 0.1 0.2 0.4 0.6 0.8;
    do
        # run PRSice
        ~/PRSice_linux \
        --base /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/snps_for_risk_prs_all_snps \
        --all-score \
        --target /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_prs_$ancestry \
        --thread 1 \
        --clump-kb 1000 \
        --clump-r2 $r2 \
        --stat OR \
        --score std \
        --or \
        --binary-target T \
        --cov ./pheno/reimputed_$ancestry\_covars_with_pcs.tsv \
        --pheno ./pheno/reimputed_$ancestry\_pheno.tsv \
        --cov-factor sex \
        --cov-col sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --bar-levels 5e-8,5e-7,5e-6,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
        --fastscore \
        --out ./outputs/prs_risk_nomhc_$r2\_$ancestry \
        --perm 10000 \
        --print-snp \
        --ultra \
        --x-range chr6:25000000-35000000

        # WITH MHC
        # run PRSice
        ~/PRSice_linux \
        --base /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/snps_for_risk_prs_all_snps \
        --all-score \
        --target /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_prs_$ancestry \
        --thread 1 \
        --clump-kb 1000 \
        --clump-r2 $r2 \
        --stat OR \
        --score std \
        --or \
        --binary-target T \
        --cov ./pheno/reimputed_$ancestry\_covars_with_pcs.tsv \
        --pheno ./pheno/reimputed_$ancestry\_pheno.tsv \
        --cov-factor sex \
        --cov-col sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --bar-levels 5e-8,5e-7,5e-6,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
        --fastscore \
        --out ./outputs/prs_risk_$r2\_$ancestry \
        --perm 10000 \
        --print-snp \
        --ultra
    done