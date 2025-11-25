#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr eur)
ancestry=${ancestries[$index]}
echo doing $ancestry 

# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# remove PCA outliers and filter to MAF 0.05 within each cluster
# filter in PLINK
~/plink2 --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_all_chroms \
--keep ./outputs/ukb_pca_non_outliers_$ancestry\.tsv \
--maf 0.05 \
--hwe 1e-5 \
--mind 0.1 \
--geno 0.1 \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19

# find freq differences 
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19 \
--test-missing \
--allow-no-sex \
--pheno ./pheno/susceptibility_$ancestry\_pheno.tsv \
--out $ancestry\_diff_missing
awk 'NR>1{if($5<0.05) print $2}' $ancestry\_diff_missing.missing > /data/scratch/hmy117/diff_missing_snps_$ancestry

# find diff missing SNPs
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19 \
--test-missing \
--allow-no-sex \
--pheno ./pheno/susceptibility_$ancestry\_pheno.tsv \
--out $ancestry\_diff_missing
awk 'NR>1{if($5<0.05) print $2}' $ancestry\_diff_missing.missing > /data/scratch/hmy117/diff_missing_snps_$ancestry


# filter out these snps
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19 \
--exclude /data/scratch/hmy117/diff_missing_snps_$ancestry \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19_qc \
--make-bed

