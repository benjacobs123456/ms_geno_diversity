#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N run_vep
#$ -o /data/scratch/hmy117


cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

module load R

# liftover
for ancestry in sas afr eur
    do
        # run liftover
        Rscript ./scripts/liftover_hg19_to_hg38.R ./outputs/case_control_gwas_ALL_ANCESTRY_$ancestry\_MS_status.regenie
        Rscript ./scripts/liftover_plink_to_hg38.R ./outputs/case_control_gwas_ALL_ANCESTRY_glm$ancestry\.MS_status.glm.logistic.hybrid
    done

module unload R
module load miniforge
mamba activate vep_env

for ancestry in sas afr eur
do
    # prepare for vep
	  awk 'NR>1{if($12<1e-3) print $1,$15,$15,$3"/"$4,"+",$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm$ancestry\.MS_status.glm.logistic.hybrid_hg38 > ./outputs/susceptibility_vep_input_$ancestry
    #awk 'NR>1{if($15>3) print $1,$18,$18,$3"/"$4,"+",$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_$ancestry\_MS_status.regenie_hg38 > ./outputs/susceptibility_vep_input_$ancestry

    # annotate with vep
    vep -i ./outputs/susceptibility_vep_input_$ancestry -o ./outputs/snp_annotations_nearest_susceptibility_$ancestry \
    --cache \
    --dir_cache /data/scratch/hmy117/.vep \
    --force_overwrite \
    --nearest symbol \
    --tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence"


done

for ancestry in sas afr eur
do
        # frequencies
    vep -i ./outputs/susceptibility_vep_input_$ancestry -o ./outputs/snp_annotations_freqs_susceptibility_$ancestry \
    --cache \
    --dir_cache /data/scratch/hmy117/.vep \
    --force_overwrite \
    --af_gnomadg \
    --af \
    --tab
done
