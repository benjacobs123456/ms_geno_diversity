#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N download_hgdp_genomes
#$ -o /data/scratch/hmy117
#$ -t 1:22

module load plink/1.9-beta6.27-gcc-12.2.0

cd /data/scratch/hmy117/hgdp_1kg_genomes/

plink --vcf gnomad.genomes.v3.1.2.hgdp_tgp.chr${SGE_TASK_ID}\.vcf.bgz \
--double-id \
--extract /data/scratch/hmy117/ukb_snps_to_keep.tsv \
--make-bed \
--out filtered_1kg_hgdp_1kg_chr${SGE_TASK_ID}
