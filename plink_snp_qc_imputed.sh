#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_qc
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/scratch/hmy117/adams_imputed/

~/plink2 \
--vcf chr${SGE_TASK_ID}\.dose.vcf.gz dosage=HDS \
--geno 0.1 \
--double-id \
--mac 1 \
--rm-dup exclude-all \
--hwe 1e-10 \
--extract /data/scratch/hmy117/adams_snps_to_keep.tsv \
--make-bed \
--out ADAMS_imputed_tmp_qc_chr${SGE_TASK_ID} 

~/plink2 \
--bfile ADAMS_imputed_tmp_qc_chr${SGE_TASK_ID} \
--new-id-max-allele-len 9999 \
--set-all-var-ids @:#:\$r:\$a \
--out ADAMS_imputed_qc_chr${SGE_TASK_ID} \
--make-bed