#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N plink_qc_step2
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/scratch/hmy117/adams_imputed/
for anc in sas afr eur;
    do 
        # do QC
        ~/plink2 \
        --vcf /data/scratch/hmy117/imputed_risk_gwas_genotypes_$anc\_chr${SGE_TASK_ID}\.vcf.gz dosage=HDS \
	--id-delim \
	--hard-call-threshold 0.49 \
        --maf 0.05 \
        --rm-dup exclude-all \
        --make-bed \
        --threads 1 \
        --extract-if-info R2 '>'= 0.7 \
        --out /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr${SGE_TASK_ID}

    done
