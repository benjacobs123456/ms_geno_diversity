#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N pca_per_ancestry
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(CSA AFR EUR)
ancestries_2=(sas afr eur)
ancestry=${ancestries[$index]}
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/


# find mhc snps
awk '{if($1==6 && $4>25000000 && $4<35000000) print $2}' /data/scratch/hmy117/matched_case_control_$ancestry_out\.bim > mhc_snps_to_exclude_pca_$ancestry_out

# LD prune, remove MHC SNPs, and restrict to good SNPs for step 1
~/plink --bfile /data/scratch/hmy117/matched_case_control_$ancestry_out \
--indep-pairwise 1000 500 0.05 \
--exclude mhc_snps_to_exclude_pca_$ancestry_out \
--maf 0.1 \
--geno 0.01 \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_independent_snps_$ancestry_out

# restrict to these snps 
~/plink --bfile /data/scratch/hmy117/matched_case_control_$ancestry_out \
--extract /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_independent_snps_$ancestry_out\.prune.in \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_pca_qcd_$ancestry_out

# run PCA
 ~/flashpca_x86-64 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_pca_qcd_$ancestry_out \
--ndim 10 \
--outpc /data/scratch/hmy117/pcs_$ancestry_out \
--outvec /data/scratch/hmy117/eigenvecs_$ancestry_out \
--outval /data/scratch/hmy117/eigenvals_$ancestry_out \
--outpve /data/scratch/hmy117/pve_$ancestry_out \


module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/pca_match.R $ancestry $ancestry_out
