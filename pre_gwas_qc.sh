#!/bin/bash
#$ -pe smp 4
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

# restrict to case-control cohort
# remove duplicates
# restrict to good snps
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$ancestry_out \
--keep ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--extract /data/scratch/hmy117/snps_to_keep_for_gwas$ancestry_out\.tsv \
--maf 0.05 \
--hwe 1e-5 \
--geno 0.1 \
--out /data/scratch/hmy117/matched_case_control_$ancestry_out \
--make-bed 

# find differential freqs 
~/plink -bfile /data/scratch/hmy117/matched_case_control_$ancestry_out \
--freq case-control \
--allow-no-sex \
--pheno ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--out /data/scratch/hmy117/case_control_diff_freqs_$ancestry_out\.tsv

awk '{if($5<0.05 || $6 <0.05) print}' /data/scratch/hmy117/case_control_diff_freqs_$ancestry_out\.tsv.frq.cc > rare_snps_to_exclude$ancestry_out

# diff miss 
~/plink -bfile /data/scratch/hmy117/matched_case_control_$ancestry_out \
--test-missing \
--allow-no-sex \
--pheno ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--out /data/scratch/hmy117/case_control_diff_miss_$ancestry_out\.tsv

awk '{if($3<0.1 && $4 <0.1 && $5 > 1e-5) print}' /data/scratch/hmy117/case_control_diff_miss_$ancestry_out\.tsv.missing > missing_snps_to_keep$ancestry_out

# filter on these snps 
~/plink -bfile /data/scratch/hmy117/matched_case_control_$ancestry_out \
--extract  missing_snps_to_keep$ancestry_out \
--exclude rare_snps_to_exclude$ancestry_out \
--out /data/scratch/hmy117/case_control_postqc_$ancestry_out \
--make-bed 

# prune for pca 
~/plink2 --bfile /data/scratch/hmy117/case_control_postqc_$ancestry_out \
--extract /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_independent_snps_$ancestry_out\.prune.in \
--out /data/scratch/hmy117/matched_case_control_for_pca_$ancestry_out \
--make-bed

# repeat pca 
~/flashpca_x86-64 --bfile /data/scratch/hmy117/matched_case_control_for_pca_$ancestry_out \
--ndim 2 \
--outpc /data/scratch/hmy117/matched_case_control_pcs_$ancestry_out

# remake covar & pheno files 
module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/make_matched_covar_file.R $ancestry $ancestry_out

# lightly LD prune for step 1
~/plink --bfile /data/scratch/hmy117/case_control_postqc_$ancestry_out \
--indep-pairwise 1000 500 0.8 \
--maf 0.1 \
--geno 0.01 \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_independent_snps_matched_$ancestry_out
