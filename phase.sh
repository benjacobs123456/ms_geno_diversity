#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N phase
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

for ancestry_out in sas afr;
    do
        # filter plink file
        ~/plink --bfile /data/scratch/hmy117/risk_gwas_genotypes_$ancestry_out \
        --keep ./pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
        --out /data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry_out\_chr${SGE_TASK_ID} \
        --chr ${SGE_TASK_ID} \
        --make-bed

        # phasing 

        ~/Eagle_v2.4.1/eagle \
        --bfile=/data/scratch/hmy117/risk_gwas_genotypes_for_phase_$ancestry_out\_chr${SGE_TASK_ID} \
        --geneticMapFile=/data/home/hmy117/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
        --outPrefix=/data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID} \
        --pbwtIters=1 \
        --numThreads=10
    done 
    
# make reference 
~/Minimac4/bin/minimac4 \
--compress-reference /data/Wolfson-UKBB-Dobson/1kg_reference/vcf_files/1kg_chr${SGE_TASK_ID}\.recode.vcf \
--output /data/Wolfson-UKBB-Dobson/1kg_reference/compressed_chr${SGE_TASK_ID}\.msav 

