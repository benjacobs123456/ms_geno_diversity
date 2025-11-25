#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N hla_assoc
#$ -o /data/scratch/hmy117
#$ -t 1:2


# initialise
index=$((SGE_TASK_ID-1))
ancestries_2=(sas afr eur)
ancestry=${ancestries_2[$index]}

echo doing $ancestry


cd /data/home/hmy117/HLA-TAPAS
module load R/3.6.1
module load python/3.8


# write pheno files
awk 'NR>1{print 0,$2,$3-1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry\_pheno.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry\_pheno_for_omnibus.tsv

awk 'NR==1{print "FID","IID","MS_status"};NR>1{print 0,$2,$3}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry\_pheno.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv

# write covar files
awk 'NR==1{print "FID",$2,$3,$4,$5,$6,$7,$8};NR>1{print 0,$2,$3,$4,$5,$6,$7,$8}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry\_covars_with_pcs.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv

awk 'NR==1{print "FID",$2,$3,$4,$5,$6,$7,$8};NR>1{print 0,$2,$3,$4,$5,$6,$7,$8}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry\_covars_with_pcs.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_logistic.tsv

# make plink file with phased vcf
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--make-bed \
--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_logistic.tsv \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus \
--threads $NSLOTS

# run omnibus logistic model
python3 -m HLAassoc OMNIBUS_LOGISTIC \
--vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--bim /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus.bim \
--fam  /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus.fam \
--covars /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_omnibus_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_omnibus.tsv \
--aa-only \
--maf-threshold 0.01

# standard analysis
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_logistic.tsv \
--glm  hide-covar \
--maf 0.01 \
--threads $NSLOTS
