#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N get_imputation_metrics
#$ -o /data/scratch/hmy117
#$ -t 1:22

module load bcftools
cd /data/scratch/hmy117/adams_imputed/
for anc in sas afr;
    do 

        # save imputation QC metrics
        bcftools query -f "%CHROM %POS %REF %ALT %INFO/MAF %INFO/R2" \
                /data/scratch/hmy117/imputed_risk_gwas_genotypes_$anc\_chr${SGE_TASK_ID}\.vcf.gz > \
                /data/scratch/hmy117/snp_qc_chr${SGE_TASK_ID}\_anc_$anc 
    done
