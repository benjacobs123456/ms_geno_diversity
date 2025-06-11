#!/bin/bash
#$ -pe smp 16
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N split_ancestries
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

# find people of inferred ancestry x
awk -v anc=$ancestry 'NR>1{if($3==anc) print $1,$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_all_ukb_adams.tsv > /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_pca_non_outliers_$ancestry\.tsv

# remove PCA outliers and filter to MAF 0.01 within each cluster
# filter in PLINK
for i in {22..1};
do
echo doing chrom $i
~/plink2 --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_chr$i \
--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_pca_non_outliers_$ancestry\.tsv \
--maf 0.01 \
--hwe 1e-20 \
--geno 0.1 \
--make-bed \
--chr $i \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19_chr$i
done

# merge all chroms for gwas
# merge all chroms
echo "Merging all chroms"

cd /data/scratch/hmy117/
rm filelist_for_gwas_merge_$ancestry_out
for i in {2..22}; do echo /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19_chr$i >> filelist_for_gwas_merge_$ancestry_out; done
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19_chr1 \
--merge-list filelist_for_gwas_merge_$ancestry_out \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$ancestry_out
