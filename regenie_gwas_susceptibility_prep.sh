#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N regenie_gwas_prep
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(CSA AFR EUR)
ancestries_2=(sas afr eur)
ancestry=${ancestries[$index]}
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry


# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# remove PCA outliers and filter to MAF 0.05 within each cluster
# filter in PLINK
for i in {22..1};
do
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry_out\_hg19_chr$i \
--keep ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--make-bed \
--maf 0.05 \
--hwe 1e-5 \
--geno 0.1 \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19_qc_chr$i
done

# merge all chroms for gwas
# merge all chroms
echo "Merging all chroms"

cd /data/scratch/hmy117/
rm filelist_for_gwas_merge_$ancestry_out
for i in {2..22}; do echo /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19_qc_chr$i >> filelist_for_gwas_merge_$ancestry_out; done
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_$ancestry\_hg19_qc_chr1 \
--merge-list filelist_for_gwas_merge_$ancestry_out \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out

# check for duplicates
~/king -b /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out\.bed \
--duplicate \
--prefix /data/scratch/hmy117/dups_$ancestry_out

awk 'NR>1{print $1,$2}' /data/scratch/hmy117/dups_$ancestry_out\.con > /data/scratch/hmy117/dups_$ancestry_out\_to_remove

# calculate freqs
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry_out\_pheno.tsv \
--allow-no-sex \
--freq case-control \
--out /data/scratch/hmy117/freqs_$ancestry_out

# calculate differential missingness
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry_out\_pheno.tsv \
--allow-no-sex \
--test-missing \
--out /data/scratch/hmy117/missingness_$ancestry_out

# find good snps
module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/find_good_snps_for_gwas.R $ancestry_out
