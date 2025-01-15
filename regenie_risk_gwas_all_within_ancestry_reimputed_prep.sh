#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_all_ancestry_prep_imputed
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

# find mhc snps
awk '{if($1==6 && $4>25000000 && $4<35000000) print $2}' /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS.bim > mhc_snps_to_exclude_pca_$ancestry_out

# LD prune, remove MHC SNPs, and restrict to good SNPs for step 1
~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--indep-pairwise 1000 500 0.1 \
--exclude mhc_snps_to_exclude_pca_$ancestry_out \
--maf 0.1 \
--geno 0.01 \
--out /data/scratch/hmy117/reimputed_risk_gwas_genotypes_pca_snps_$ancestry_out \
--threads $NSLOTS

# restrict to these snps 
~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
--extract /data/scratch/hmy117/reimputed_risk_gwas_genotypes_pca_snps_$ancestry_out\.prune.in \
--make-bed \
--out /data/scratch/hmy117/reimputed_risk_gwas_genotypes_for_pca_$ancestry_out \
--threads $NSLOTS

# run PCA
 ~/flashpca_x86-64 --bfile /data/scratch/hmy117/reimputed_risk_gwas_genotypes_for_pca_$ancestry_out \
--ndim 10 \
--outpc /data/scratch/hmy117/pcs_$ancestry_out \
--outvec /data/scratch/hmy117/eigenvecs_$ancestry_out \
--outval /data/scratch/hmy117/eigenvals_$ancestry_out \
--outpve /data/scratch/hmy117/pve_$ancestry_out \
--outload /data/scratch/hmy117/snp_loadings_$ancestry_out

module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/pca_outlier_removal_reimputed.R $ancestry $ancestry_out





