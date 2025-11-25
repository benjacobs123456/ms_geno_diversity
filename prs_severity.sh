#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N prs_sev
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr eur)
ancestry=${ancestries[$index]}

echo doing $ancestry


# run PRS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/


# filter in PLINK
~/plink2 --bfile ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--set-all-var-ids chr@:# \
--out ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38_chrpos \
--make-bed

for r2 in 0.01 0.1 0.2 0.4 0.6 0.8;
do
# run PRSice
~/PRSice_linux \
--base snps_for_severity_prs_all_snps \
--all-score \
--target ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38_chrpos \
--thread 1 \
--clump-kb 1000 \
--clump-r2 $r2 \
--stat BETA \
--score std \
--beta \
--binary-target F \
--pheno ./pheno/$ancestry\_prs_pheno.tsv \
--cov ./pheno/$ancestry\_covars_with_pcs.tsv \
--cov-factor Sex,Site \
--pheno-col rint_gARMSS \
--bar-levels 5e-8,5e-7,5e-6,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
--fastscore \
--out ./outputs/prs_$r2\_$ancestry \
--perm 10000 \
--print-snp \
--ultra
done

# run for other phenotypes
~/PRSice_linux \
--base snps_for_severity_prs_all_snps \
--all-score \
--target ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38_chrpos \
--thread 1 \
--clump-kb 1000 \
--clump-r2 0.1 \
--stat BETA \
--score std \
--beta \
--pheno ./pheno/$ancestry\_prs_pheno.tsv \
--pheno-col rint_gARMSS,rint_edss,rint_msis_physical_normalised,rint_eq5d_vas,rint_age_at_dx \
--cov ./pheno/$ancestry\_covars_with_pcs.tsv \
--cov-factor Sex,Site \
--bar-levels 5e-8,5e-7,5e-6,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
--fastscore \
--out ./outputs/prs_all_phenos_$ancestry \
--perm 10000 \
--print-snp \
--ultra
