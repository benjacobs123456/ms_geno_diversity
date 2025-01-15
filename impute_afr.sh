#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N phase
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

module load bcftools 
module load htslib 

ancestry_out="afr"

# convert haps to vcf 
~/plink2 --haps /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.haps.gz \
--sample /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.sample \
--export vcf \
--out /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID} \
--threads $NSLOTS

# index 
bgzip -c /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.vcf > /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.vcf.gz
tabix -p vcf /data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.vcf.gz

# imputation 
~/Minimac4/bin/minimac4 \
/data/Wolfson-UKBB-Dobson/1kg_reference/msav_files/compressed_chr${SGE_TASK_ID}\.msav \
/data/scratch/hmy117/phased_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.vcf.gz \
-o /data/scratch/hmy117/imputed_risk_gwas_genotypes_$ancestry_out\_chr${SGE_TASK_ID}\.vcf.gz


