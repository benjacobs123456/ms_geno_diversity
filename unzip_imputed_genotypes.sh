#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N unzip
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/adams_imputed
unzip -P joGUiWvc=?C57B -o chr_${SGE_TASK_ID}\.zip
