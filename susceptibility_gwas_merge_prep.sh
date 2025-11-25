#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N regenie_gwas
#$ -o /data/scratch/hmy117
#$ -t 1:2


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


# step 1 REGENIE
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
	~/regenie \
--step 1 \
--cc12 \
--bed /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--covarFile ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem --spa \
--force-step1 \
--pThresh 0.05 \
--out ./outputs/case_control_step1_$ancestry_out\_fit \
--catCovarList sex \
--covarColList age,PC{1:50}

# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--bed /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--covarFile ./pheno/susceptibility_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem --spa \
--pThresh 0.05 \
--pred ./outputs/case_control_step1_$ancestry_out\_fit_pred.list \
--catCovarList Sex \
--catCovarList sex \
--covarColList age,PC{1:50}
--minMAC 20 \
--af-cc
