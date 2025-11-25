#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N convert_1kg
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117
module load plink/1.9-170906

# download 
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

wait 

# convert to plink 
plink \
--vcf ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
--make-bed \
--out 1kg_chr${SGE_TASK_ID}

