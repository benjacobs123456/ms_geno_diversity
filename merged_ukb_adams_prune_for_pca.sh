#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N prune_for_pca
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/

# prune
~/plink --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID} \
--allow-no-sex \
--extract /data/scratch/hmy117/hgdp_1kg_genomes/combined_hgdp_1kg_unrelated_hg19_cpra.bim \
--indep-pairwise 1000 500 0.01 \
--out /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_pruned

# filter to pruned markers
~/plink --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID} \
--allow-no-sex \
--extract /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_pruned.prune.in \
--out /data/scratch/hmy117/ukb_adams_merged_genotypes_chr${SGE_TASK_ID}\_pruned_for_pca \
--make-bed
