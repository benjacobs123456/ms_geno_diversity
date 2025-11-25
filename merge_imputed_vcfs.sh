# merge chroms
#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N merge_vcfs
#$ -o /data/scratch/hmy117

cd /data/scratch/hmy117/adams_imputed/

module load bcftools
rm filelist
for i in {22..1}; do echo chr$i\_filtered.vcf.gz >> filelist; done
bcftools concat --file-list filelist -Ov -o merged_all_chroms
