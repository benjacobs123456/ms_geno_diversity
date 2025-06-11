#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N reimp_metrics
#$ -o /data/scratch/hmy117
#$ -t 1:3


# initialise
index=$((SGE_TASK_ID-1))
ancestries_2=(sas afr eur)
ancestry_out=${ancestries_2[$index]}

echo doing $ancestry

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
module load R
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/reimputation_metrics.R $ancestry_out
