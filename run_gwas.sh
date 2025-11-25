# merge chroms
#!/bin/bash
#$ -pe smp 25
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N gwas
#$ -o /data/scratch/hmy117

# set wd
cd ~/data/scratch/hmy117/

module load plink/2.0-20220603

plink2 --bfile /data/scratch/hmy117/ms_ukb_merged_qc3_sas \
--glm hide-covar \
--covar covars.txt \
--covar-variance-standardize



