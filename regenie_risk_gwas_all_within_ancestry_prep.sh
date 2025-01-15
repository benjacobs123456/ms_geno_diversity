#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_all_ancestry
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

# check for duplicates
~/king -b /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$ancestry_out\.bed \
--duplicate \
--prefix /data/scratch/hmy117/dups_$ancestry_out

awk 'NR>1{print $1,$2}' /data/scratch/hmy117/dups_$ancestry_out\.con > /data/scratch/hmy117/dups_$ancestry_out\_to_remove



# case-case gwas
# exclude SNPs from case-case GWAS
awk '{if($3=="MS") print $1,$1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv > ms_cases_$ancestry_out
awk '{if($3=="MS" && $2=="UKB") print $1,$1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv > ms_cases_ukb_$ancestry_out

# filter plink file 
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$ancestry_out \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--keep ms_cases_$ancestry_out \
--out just_ms_$ancestry_out \
--make-bed \
--make-pheno ms_cases_ukb_$ancestry_out \* \
--allow-no-sex

# fisher's test 
~/plink --bfile just_ms_$ancestry_out \
--assoc \
--allow-no-sex \
--out fisher_case_case_$ancestry_out

# with counts
~/plink --bfile just_ms_$ancestry_out \
--assoc counts \
--allow-no-sex \
--out fisher_case_case_counts_$ancestry_out

# explore case-case gwas 

# exclude SNPs which 
## show association with cohort 
## MAF < 5% in either cohort
module load R/4.2.2
Rscript ./scripts/exclude_snps_case_case.gwas.R $ancestry_out

echo "Keeping this many SNPs from case-case GWAS:"
awk 'END{print NR}' snps_to_keep_case_case_gwas_$ancestry_out 

# exclude bad SNPs
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$ancestry_out \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_$ancestry_out \
--extract snps_to_keep_case_case_gwas_$ancestry_out \
--remove /data/scratch/hmy117/dups_$ancestry_out\_to_remove \
--make-bed 

# differential missingness 
awk '{if($2=="ADAMS") print $1,$1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv > adams_cases_$ancestry_out

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_$ancestry_out \
--make-bed \
--make-pheno adams_cases_$ancestry_out \* \
--allow-no-sex \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_cohort_pheno_$ancestry_out 

# calculate missingess 
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_cohort_pheno_$ancestry_out \
--test-missing \
--allow-no-sex \
--out missingness_cohort_$ancestry_out 

awk '{if($3 > 0.1 || $4 > 0.1) print}' missingness_cohort_$ancestry_out\.missing > missing_snps_$ancestry_out 

# frequencies 
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_cohort_pheno_$ancestry_out \
--freq case-control \
--allow-no-sex \
--out freqs_cohort_$ancestry_out 

awk '{if($5 > 0.01 && $6 > 0.01) print}' freqs_cohort_$ancestry_out\.frq.cc > common_snps_$ancestry_out 

# filter to these snps 
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_postqc_$ancestry_out \
--extract common_snps_$ancestry_out \
--exclude missing_snps_$ancestry_out \
--allow-no-sex \
--out /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--make-bed


# find mhc snps
awk '{if($1==6 && $4>25000000 && $4<35000000) print $2}' /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out\.bim > mhc_snps_to_exclude_pca_$ancestry_out

# LD prune, remove MHC SNPs, and restrict to good SNPs for step 1
~/plink --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--indep-pairwise 1000 500 0.1 \
--exclude mhc_snps_to_exclude_pca_$ancestry_out \
--maf 0.1 \
--geno 0.01 \
--out /data/scratch/hmy117/risk_gwas_genotypes_pca_snps_$ancestry_out

# restrict to these snps 
~/plink --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--extract /data/scratch/hmy117/risk_gwas_genotypes_pca_snps_$ancestry_out\.prune.in \
--make-bed \
--out /data/scratch/hmy117/risk_gwas_genotypes_for_pca_$ancestry_out

# run PCA
 ~/flashpca_x86-64 --bfile /data/scratch/hmy117/risk_gwas_genotypes_for_pca_$ancestry_out \
--ndim 10 \
--outpc /data/scratch/hmy117/pcs_$ancestry_out \
--outvec /data/scratch/hmy117/eigenvecs_$ancestry_out \
--outval /data/scratch/hmy117/eigenvals_$ancestry_out \
--outpve /data/scratch/hmy117/pve_$ancestry_out \
--outload /data/scratch/hmy117/snp_loadings_$ancestry_out

module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/pca_outlier_removal.R $ancestry $ancestry_out




