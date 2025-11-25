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
ancestries=(CSA AFR EUR)
ancestries_2=(sas afr eur)
ancestry=${ancestries[$index]}
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# remove PCA outliers and filter to MAF 0.05 within each cluster
# filter in PLINK

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19 \
--indep-pairwise 1000 500 0.01 \
--allow-no-sex \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_$ancestry_out


# restrict to pruned snps
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19 \
--extract  /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_$ancestry_out.prune.in \
--make-bed \
--allow-no-sex \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_for_pca_$ancestry_out

# calculate PCs
~/plink2 --bfile /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_for_pca_$ancestry_out \
--pca approx 20 \
--out /data/scratch/hmy117/adams_ukb_pca_projections_$ancestry_out

