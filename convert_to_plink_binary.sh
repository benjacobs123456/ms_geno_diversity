#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N convert_to_plink1
#$ -o /data/scratch/hmy117


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
module load plink/1.9-170906


grep Control ./raw/adams_250823.ped | cut -f1,2 > ./raw/controls_to_remove

plink --file /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/raw/adams_250823 \
--make-bed \
--chr 1-22 \
--geno 0.1 \
--mind 0.1 \
--maf 0.0001 \
--hwe 1e-10 \
--remove ./raw/controls_to_remove \
--out ./outputs/ADAMS_geno
