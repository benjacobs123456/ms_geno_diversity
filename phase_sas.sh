#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N phase
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

ancestry_out="sas"
# filter plink file
~/plink --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
--keep ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
--out /data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry_out\_chr${SGE_TASK_ID} \
--chr ${SGE_TASK_ID} \
--make-bed \
--threads $NSLOTS

# phasing

~/Eagle_v2.4.1/eagle \
--bfile=/data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry_out\_chr${SGE_TASK_ID} \
--geneticMapFile=/data/home/hmy117/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--outPrefix=/data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID} \
--pbwtIters=1 \
--numThreads=$NSLOTS
