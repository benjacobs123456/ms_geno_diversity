#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N unzip
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/adams_imputed/

unzip -P wSqJCdOGLh87xf chr_${SGE_TASK_ID}\.zip
