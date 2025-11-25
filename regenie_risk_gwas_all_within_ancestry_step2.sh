#!/bin/bash
#$ -pe smp 4
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

#wd 
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/



# step 2 REGENIE
	~/regenie \
--step 2 \
--cc12 \
--bed /data/scratch/hmy117/matched_case_control_$ancestry_out \
--covarFile ./pheno/susceptibility_ALL_$ancestry_out\_covars_with_pcs.tsv \
--phenoFile ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--bsize 1000 \
--bt --lowmem \
--pred ./outputs/case_control_step1_ALL_$ancestry_out\_fit_pred.list \
--catCovarList sex \
--covarColList age,PC{1:10} \
--minMAC 20 \
--af-cc \
--out ./outputs/case_control_gwas_ALL_ANCESTRY_$ancestry_out


