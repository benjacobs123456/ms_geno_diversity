#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N download_hgdp_genomes
#$ -o /data/scratch/hmy117
#$ -t 1:22

module load plink/1.9-170906
cd /data/scratch/hmy117/hgdp_1kg_genomes/

plink --vcf gnomad.genomes.v3.1.2.hgdp_tgp.chr${SGE_TASK_ID}\.vcf.bgz \
--double-id \
--extract adams_snps_for_filtering.tsv \
--make-bed \
--out filtered_1kg_hgdp_1kg_chr${SGE_TASK_ID}
