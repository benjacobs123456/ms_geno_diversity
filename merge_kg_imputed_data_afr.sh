#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N merge_imputed_chroms
#$ -o /data/scratch/hmy117


cd /data/scratch/hmy117/adams_imputed/

for anc in afr;
    do 
        # update bim to cpra 
        for i in {1..22};
            do 

                ~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i \
                --new-id-max-allele-len 9999 \
                --set-all-var-ids @:#:\$r:\$a \
                --make-just-bim \
                --out /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i\_cpra

                # replace bim 
                rm /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i\.bim
                mv /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i\_cpra\.bim /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i\.bim
            done
        # make merge list 
        rm merge_list_$anc
        for i in {2..22};
            do 
                echo /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr$i >> merge_list_$anc
            done
        ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_chr1 \
        --merge-list merge_list_$anc \
        --make-bed \
        --out /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
        --threads 1 

    done
