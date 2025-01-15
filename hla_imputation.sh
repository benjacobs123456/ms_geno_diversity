#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N hla_imp
#$ -o /data/scratch/hmy117
#$ -t 1:2

# initialise
index=$((SGE_TASK_ID-1))
ancestries=(sas afr)
ancestry=${ancestries[$index]}

echo doing $ancestry

cd /data/home/hmy117/HLA-TAPAS

# rename SNPs using reference marker names
module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/get_ref_names_for_hla_tapas.R \
/data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry\_chr6.bim \


# update names
~/plink2 --bfile /data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry\_chr6 \
--chr 6 \
--from-bp 25000000 \
--to-bp 35000000 \
--maf 0.05 \
--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry\_pheno.tsv \
--update-name /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry\_ALL_CHRS.bim_snp_name_updates.tsv \
--out /data/scratch/hmy117/risk_gwas_genotypes_$ancestry\_for_hla_imp \
--make-bed \
--threads $NSLOTS

# restrict to SNPs in ref panel 
~/plink2 --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry\_for_hla_imp \
--extract "/data/home/hmy117/HLA-TAPAS/resources/1000G.bglv4.bim" \
--out /data/scratch/hmy117/risk_gwas_genotypes_$ancestry\_for_hla_imp_just_ref_snps \
--make-bed 

# load java 
module load java


# impute HLA
python3 -m SNP2HLA \
--target /data/scratch/hmy117/risk_gwas_genotypes_$ancestry\_for_hla_imp_just_ref_snps \
--out /data/scratch/hmy117/imputed_hla_$ancestry \
--reference resources/1000G.bglv4 \
--nthreads 1 \
--mem 64g \
--niterations 1


