#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N hibag
#$ -o /data/scratch/hmy117
#$ -t 1:6




module load R/3.6.1
Rscript "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hibag_imputation.R" ${SGE_TASK_ID} sas
