#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N admixture
#$ -o /data/scratch/hmy117

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/
~/dist/admixture_linux-1.3.0/admixture \
ADAMS_geno_qc.bed \
5

