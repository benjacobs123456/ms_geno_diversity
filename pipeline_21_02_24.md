
# Genetic analysis of Multiple Sclerosis in diverse ancestries
### A genetic _A_ssociation Study of individuals of _D_iverse _A_ncestral backgrounds with _M_ultiple _S_clerosis (ADAMS)<a href =https://app.mantal.co.uk/adams><p><center> <img src="https://app.mantal.co.uk/files/db2440daa093e6b2d688b8a02b82cff2.png" width=20% height=40%></a></p></center>


This repo contains code used to analyze data from the ADAMS project - a genotype-phenotype cohort of people with Multiple Sclerosis focussed on those from diverse ancestral backgrounds. The baseline cohort description is [here](https://bmjopen.bmj.com/content/13/5/e071656.full).

# Code
Prior to imputation genotypes were called using Illumina Genome Studio v2.0. Raw PLINK binary files (.map and .ped) were exported from Genome Studio with calls coerced to the forward strand.

## Imputation
### Initial genotype quality control
First the .ped and .map files are converted to PLINK1 binary format with some light QC.
- Conversion to PLINK1 binary
- Basic QC pre-imputation: autosomes, missingess for genotypes and people <10%, MAC > 10, HWE deviation

````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23

~/plink --file /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/raw/ADAMS_16_01_24 \
--make-bed \
--chr 1-22 \
--geno 0.1 \
--hwe 1e-10 \
--mac 20 \
--mind 0.01 \
--out ./outputs/ADAMS_geno
````

#### Update FIDs to IIDs
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
awk '{print $1,$2,$2,$2}'  ./outputs/ADAMS_geno.fam > ./outputs/new_ids
~/plink --bfile ./outputs/ADAMS_geno \
--update-ids ./outputs/new_ids \
--make-bed \
--out ./outputs/ADAMS_geno_fid_iid
````

#### Remove people missing in phenotype data
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid \
--keep ./pheno/adams_pheno.tsv \
--make-bed \
--out ./outputs/ADAMS_geno_fid_iid_inpheno
````


### Kinship / duplicates
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/
~/king -b ADAMS_geno_fid_iid_inpheno.bed --duplicate
````

#### Download king.con & run
#### explore_duplicates.R on QM server
#### upload ids_to_exclude.tsv
#### Remove duplicates
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid_inpheno \
--remove ./outputs/ids_to_exclude.tsv \
--make-bed \
--out ./outputs/ADAMS_geno_fid_iid_inpheno_nodups

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/
~/king -b ADAMS_geno_fid_iid_inpheno_nodups.bed --duplicate
~/king -b ADAMS_geno_fid_iid_inpheno_nodups.bed --related --degree 3
````

#### Missingness
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid_inpheno_nodups \
--missing \
--out miss_check

````
#### Heterozygosity
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid_inpheno_nodups \
--het small-sample \
--out het_check
````

````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23")
het = read_table("het_check.het")


p=ggplot(het,aes(F))+
	geom_histogram()+
	theme_minimal()+
	labs(x="Heterozygosity statistic (F)",y="N")+
	geom_vline(xintercept = mean(het$F)-5*sd(het$F),linetype="dashed",color="blue")+
	geom_vline(xintercept = mean(het$F)+5*sd(het$F),linetype="dashed",color="blue")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/het_stats.png",res=900,units="in",width=4,height=4)
p
dev.off()

````
#### Sex check
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23")

fam_file = read_table("./outputs/ADAMS_geno_fid_iid_inpheno_nodups.fam",col_names=F)
covars = read_tsv("./pheno/adams_covars.tsv")
fam_file = fam_file %>%
	dplyr::rename("IID" = X1) %>%
	left_join(covars,by="IID")

fam_file %>%
	group_by(Sex) %>%
	dplyr::count(X5) %>%
	mutate(prop = n/sum(n))

# plot
ggplot(fam_file,aes(Sex,fill=factor(X5)))+
	geom_bar(position="fill",color="black")

````

### Export VCFs

````unix
# flip strand
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
for i in {1..22};
do
	~/plink --bfile ./outputs/ADAMS_geno_fid_iid_inpheno_nodups \
	--chr $i \
	--recode vcf \
	--output-chr chrMT \
	--out ./imputation_raw_files/chr$i &
done
````

### Sort with bcf tools
````unix
module load bcftools
for i in {1..22};
do
  bcftools sort ./imputation_raw_files/chr$i\.vcf \
  -Oz -o ./imputation_raw_files/sorted_chr$i\.vcf.gz &
done
````

### First pass imputation
- Attempt imputation at TOPMED server (TOPMED-R3)
- This will fail due to strand flips
- Download excluded snps and use strand flips to flip

````unix
# wget https://imputation.biodatacatalyst.nhlbi.nih.gov/results/job-20230901-122606-755/statisticDir/snps-excluded.txt

awk 'FS="\t"{if($2=="Strand flip") print $1}' ./imputation_raw_files/snps-excluded.txt > ./imputation_raw_files/snps_to_flip
awk 'FS="\t"{if($2!="Strand flip") print $1}' ./imputation_raw_files/snps-excluded.txt | uniq > ./imputation_raw_files/snps_to_exclude

# set var IDs to chr:pos:ref:alt
for i in {1..22};
do
	~/plink2 --vcf ./imputation_raw_files/chr$i\.vcf \
	--set-all-var-ids chr@:#:\$r\:\$a \
	--make-bed \
	--out ./imputation_raw_files/updated_ids_chr$i &
done

# recode VCFs
for i in {1..22};
do
	~/plink --bfile ./imputation_raw_files/updated_ids_chr$i \
	--exclude ./imputation_raw_files/snps_to_exclude \
	--out ./imputation_raw_files/flipped_chr$i \
	--output-chr chrMT \
	--recode vcf &
done

# flip strands
for i in {1..22};
do
	~/plink --vcf ./imputation_raw_files/flipped_chr$i\.vcf --double-id \
	--flip ./imputation_raw_files/snps_to_flip \
	--out ./imputation_raw_files/dbl_flipped_chr$i \
	--output-chr chrMT \
	--recode vcf &
done

	# sort with bcf tools
module load bcftools
for i in {1..22};
do
	bcftools sort ./imputation_raw_files/dbl_flipped_chr$i\.vcf \
	 -Oz -o ./imputation_raw_files/sorted_chr$i\.vcf.gz &
done
````

###########################
### Download imputed data
###########################

### 2nd pass imputation
- Performed imputation again
- Download imputed data
- Rsq filter 0.001
- TOPMED-r3

````unix
cd /data/scratch/hmy117/adams_imputed

wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/d7a25dc866a65bd63cb4a4067159d85198dc6023b01480600e7d23aede14fb16/chr_1.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/eb5b07fca297cd98b43e9b7dc83718f34e636e3ca1fa3b0abb98a74406bdbf44/chr_10.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/a1dbea2988d4dee952d906bb15f94ab0d66cbc0a05029cd367bfac3172dd6741/chr_11.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/677c338bbe78f55279fdcb089797d7f27cbe16ad75cf0f6d01559ee7f200a1fe/chr_12.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/6f4c37f14c4fbb9b25372448d2d57055b05a9ad10e8b04c848c28ca253eb1e38/chr_13.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/8c99bdabec3599885d19c4035dcaafafcade23a9f0d321c79e45283ad9f97ccb/chr_14.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/539b1d02c3ea59873a448142451eef07cfe0e57d45e463175667ca459fed37c4/chr_15.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/4f55f533b68bfd22c89ba5f2a1addfe6ea7caaff7b60189780bfdfc9d2f57268/chr_16.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/128d7c3f5489f1efd906d94f05cf1ef6220b912e62ad1f0e10428e73a6905287/chr_17.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/95a4d57f8283707acb0b6d8b85a4e38a2a48513e46a0f0795475b71079f86d88/chr_18.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/06794fb7a71606baa8ef54ad1e8af1a1a80a1adab02668a04448ebe6c6719523/chr_19.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/dfb1c399b286aeb2d40b34b78b0e550cf27796f89ac59bcf5d9c888fec356ee3/chr_2.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/592e31c53ca50a9d41cd2bb24047f9793dd79ca477fd7e5dbeed3d9590f9e01c/chr_20.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/df298e8b3eab1e8f71682280635e61c756ada0f4ce05708d46e695f8ca5bebf0/chr_21.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/0e20bb858e52b73e0ac920fe4ee7999bce67d4986d1ef6d7735bf2caab3fd562/chr_22.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/bb8981e8025ead8cd7da6779fb363aab0e762b2c9fff2a2974c8cb81f3064e08/chr_3.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/af927d519118be1fdfbb2df2bada73ac7b58bfa38b696b08b30bd4f25900afaa/chr_4.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/d3c79e79f5db2febdeb5e69310791bbef8625e16bc2b14edfee305ed09cc8567/chr_5.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/cabd159975ce321d6da95280eb0d4838b331f72bb0040cc766c28bf9a039d523/chr_6.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/282f0d22969ffe45b340cf316f75975e05b3bdce7f93bee431e8ab4087f0a399/chr_7.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/317e84d729b7f3f44f5c5f33280c00ff4798d7cc5004fca268cd9c758c4ce462/chr_8.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/eb0994cd8288bef8513d2781ae8359381ebeb8086fe93b29599512a2e483227c/chr_9.zip &

# unzip
for i in {1..22};
do
  unzip -o -P rZvZSBy4p0NiXw chr_$i\.zip.1 &
done

````

### HLA imputation
- HLA imputation run via MIS https://imputationserver.sph.umich.edu/index.html#!run/imputationserver-hla
- Note that this output is _hg19_
````unix
cd /data/scratch/hmy117/adams_imputed
wget https://imputationserver.sph.umich.edu/share/results/503cd6ac93b7dd8e208374b21b169d58fd9a1d6065f811492afed12a0e8ac2c1/chr_6.zip
unzip -P "pV)PSgZfUt4pu2" chr_6.zip

````
## Ancestry inference

### Liftover (for bigSNPR ancestry)
````unix

cd /data/scratch/hmy117/hgdp_1kg_genomes/
# lift over to hg19
awk '{print "chr"$1,$4-1,$4,$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid_inpheno_nodups.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$1":"$3":"$4}' hg19_bedfile > hg19_snps
awk '{print $1":"$3":"$4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions and IDs
~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid_inpheno_nodups \
--update-name hg19_snps \
--make-bed \
--out adams_temp_hg19_for_bigsnpr

~/plink --bfile adams_temp_hg19_for_bigsnpr \
--update-map hg19_snp_positions \
--make-bed \
--out adams_hg19_for_bigsnpr_cpra
````

### Download reference data
````unix
# run with qsub ./scripts/download_hgdp_1kg.sh
#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N download_hgdp_genomes
#$ -o /data/scratch/hmy117
#$ -t 1:22

~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${SGE_TASK_ID}\.vcf.bgz \
/data/scratch/hmy117/hgdp_1kg_genomes
````

### Filter to ADAMS variants
````R
library(tidyverse)
df = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid_inpheno_nodups.bim",col_names=F) %>%
  dplyr::select(2)
message(nrow(df)," SNPs")
df$X2 = str_remove_all(df$X2,"GSA-")
df = df %>%
  filter(grepl("^rs",X2))
message(nrow(df)," SNPs with rsids")

write_tsv(df,"/data/scratch/hmy117/hgdp_1kg_genomes/adams_snps_for_filtering.tsv")
````
### Filter HGDP genomes to ADAMS genotyped variants
````unix
#qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_hgdp_to_adams_vars.sh
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
````

### Download ancestry calls
````unix
~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/data_intersection/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv \
/data/scratch/hmy117/hgdp_1kg_genomes/
````

### Merge hgdp-1kg
````unix
cd /data/scratch/hmy117/hgdp_1kg_genomes/
for i in {2..22};
  do
    echo filtered_1kg_hgdp_1kg_chr$i >> merge_filelist
  done

cd /data/scratch/hmy117/hgdp_1kg_genomes/

# try merge
~/plink --bfile filtered_1kg_hgdp_1kg_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg

# remove duplicates & do QC
for i in {1..22};
  do
    ~/plink --bfile filtered_1kg_hgdp_1kg_chr$i \
    --exclude combined_hgdp_1kg-merge.missnp \
    --make-bed \
    --out hgdp_1kg_nodups_chr$i \
    --maf 0.05 \
    --geno 0.01 \
    --hwe 1e-10
  done

# try merge again
rm merge_filelist
for i in {2..22};
  do
    echo hgdp_1kg_nodups_chr$i >> merge_filelist
  done
~/plink --bfile hgdp_1kg_nodups_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg_filtered

````

### Do same QC for adams
````unix

cd /data/scratch/hmy117/hgdp_1kg_genomes/
~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid_inpheno_nodups \
--extract combined_hgdp_1kg.bim \
--out ADAMS_qc_genos_for_pcs \
--make-bed
````

### Missingness filter
````unix
~/plink --bfile combined_hgdp_1kg_filtered \
--mind 0.1 \
--out combined_hgdp_1kg_nonmissing \
--make-bed
````

### Get intersection of nonpalindromic, compatible snps
````R
library(tidyverse)
adams_snps = read_table("ADAMS_qc_genos_for_pcs.bim",col_names=F)
kg_snps = read_table("combined_hgdp_1kg_nonmissing.bim",col_names=F)

# filter
adams_snps = adams_snps %>% filter(X2 %in% kg_snps$X2)
kg_snps = kg_snps %>% filter(X2 %in% adams_snps$X2)

# merge
combo = adams_snps %>% left_join(kg_snps,by="X2")

# filter to compatible alleles & filter out palindromes
combo = combo %>%
  filter(
      (X5.x == X5.y & X6.x == X6.y ) |
      (X5.x == X6.y & X6.x == X5.y )
      ) %>%
    filter(
      ! ( X5.x == "G" & X5.y == "C") &
      ! ( X5.x == "C" & X5.y == "G") &
      ! ( X5.x == "A" & X5.y == "T") &
      ! ( X5.x == "T" & X5.y == "A")
      )

combo = combo %>% dplyr::select(X2)
write_tsv(combo,"snps_to_keep.tsv",col_names=F)
````

### Project PCs
````unix

# rename as cpra
~/plink2 --bfile combined_hgdp_1kg_nonmissing \
--extract snps_to_keep.tsv \
--make-bed \
--out combined_hgdp_1kg_hg38

~/plink2 --bfile combined_hgdp_1kg_hg38 \
--set-all-var-ids @:#:\$r\:\$a \
--make-bed \
--out combined_hgdp_1kg_hg38_cpra

# filter to compatible snps
# remove related indivs
~/king -b combined_hgdp_1kg_hg38_cpra.bed --degree 3
awk 'NR>1{print $1,$2}' king.kin0 > indivs_to_remove
~/plink --bfile combined_hgdp_1kg_hg38_cpra \
--remove indivs_to_remove \
--out combined_hgdp_1kg_unrelated \
--make-bed

# prune 1kg
~/plink --bfile combined_hgdp_1kg_unrelated \
--indep-pairwise 1000 100 0.1 \
--out pruned_snps_hgdp_1kg

# filter
~/plink --bfile combined_hgdp_1kg_unrelated \
--out pruned_hgdp_1kg \
--extract pruned_snps_hgdp_1kg.prune.in \
--make-bed

# calculate PCs
~/plink2 --bfile pruned_hgdp_1kg \
--pca allele-wts 50 \
--freq \
--out hgdp_kg_pcs

# rename ADAMS as cpra
~/plink2 --bfile ADAMS_qc_genos_for_pcs \
--set-all-var-ids @:#:\$r\:\$a \
--make-bed \
--out adams_geno_for_pcs_hg38_cpra

# project ADAMS samples
~/plink2 --bfile adams_geno_for_pcs_hg38_cpra \
--read-freq hgdp_kg_pcs.afreq \
--score hgdp_kg_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize list-variants \
--score-col-nums 6-55 \
--out adams_pcs

# project original dataset
~/plink2 --bfile pruned_hgdp_1kg \
--read-freq hgdp_kg_pcs.afreq \
--extract adams_pcs.sscore.vars \
--score hgdp_kg_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize \
--score-col-nums 6-55 \
--out hgdp_kg_pcs_rescored

````

### Download metadata
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv


### Ancestry inference
````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")

# read in data
adams = read_table("adams_pcs.sscore")
kg_hgdp = read_table("hgdp_kg_pcs_rescored.sscore")
meta = read_tsv("gnomad_meta_v1.tsv") %>%
  dplyr::select(2,hgdp_tgp_meta.Population,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("IID","pop","superpop")
kg_hgdp = kg_hgdp %>%
left_join(meta,by="IID")

# build RF on HGDP data
# filter out those with missing labels
kg_hgdp = kg_hgdp %>%
  filter(!is.na(superpop)) %>%
  dplyr::select(superpop,contains("PC"))
library(caret)

rf_fit = train(superpop ~ .,
                      data=kg_hgdp,
                      method='rf',
                      metric='Accuracy')


# save rf
saveRDS(rf_fit,"rf_fit.rds")

# can reload here
# rf_fit = readRDS("rf_fit.rds")

# predict
adams$predicted_ancestry = predict(rf_fit,adams)

# count
counts = adams %>% dplyr::count(predicted_ancestry) %>% mutate(pct = n/sum(n))
counts = counts %>% arrange(desc(n))
counts$predicted_ancestry= factor(counts$predicted_ancestry,levels = counts$predicted_ancestry,ordered=T)
p=ggplot(counts,aes(n,predicted_ancestry,fill=predicted_ancestry,
	label=paste0(round(pct*100,2),"%")))+
	geom_col(color="black")+
	scale_fill_brewer(palette="Set1")+
	theme_minimal()+
	geom_text(hjust=0,position=position_nudge(x=10))+
	scale_x_continuous(limits = c(0,400))+
	labs(x="N",y="Inferred genetic ancestry")+
	theme(legend.position="none")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_counts.png",res=900,units="in",width=4,height=3)
p
dev.off()

adams %>%
	dplyr::count(predicted_ancestry) %>%
	mutate(total = sum(n), prop = n/total)

p=ggplot(adams,aes(PC1_AVG,PC2_AVG,col=predicted_ancestry))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Inferred genetic ancestry")

p2=ggplot(adams,aes(PC1_AVG,PC2_AVG,col=predicted_ancestry))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Genetic ancestry")+
  ggtitle("ADAMS")+
	scale_x_continuous(limits = c(-0.2,0.2))+
	scale_y_continuous(limits = c(-0.15,0.15))

p3 = ggplot(kg_hgdp %>% na.omit(),aes(PC1_AVG,PC2_AVG,col=superpop))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Genetic ancestry")+
  ggtitle("1kg + HGDP reference")+
	scale_x_continuous(limits = c(-0.2,0.2))+
	scale_y_continuous(limits = c(-0.15,0.15))

p4=ggplot(adams,aes(PC3_AVG,PC4_AVG,col=predicted_ancestry))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC3",y="PC4",col="Genetic ancestry")+
  ggtitle("ADAMS")+
	scale_y_continuous(limits = c(-0.25,0.05))+
	scale_x_continuous(limits = c(-0.15,0.15))
p5 = ggplot(kg_hgdp %>% na.omit(),aes(PC3_AVG,PC4_AVG,col=superpop))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC3",y="PC4",col="Genetic ancestry")+
  ggtitle("1kg + HGDP reference")+
	scale_y_continuous(limits = c(-0.25,0.05))+
	scale_x_continuous(limits = c(-0.15,0.15))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pc_plot_with_ref_pc1_4.png",res=900,units="in",width=10,height=8)
cowplot::plot_grid(p2,p4,p3,p5,align="v",ncol=2)
dev.off()

# write all ancestry calls
anc_calls = adams %>%
  dplyr::select(1,2,predicted_ancestry)
write_tsv(anc_calls,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls.tsv")

# combine with self-reported ethnicity
covars = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv")
adams = adams %>%
	left_join(covars,by="IID")

# plot

ethnicity_ancestry_counts = adams %>% group_by(ethnicity_clean) %>%
dplyr::count(predicted_ancestry) %>%
mutate(prop = n/sum(n),total = sum(n))
ethnicity_ancestry_counts
write_csv(ethnicity_ancestry_counts,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/anc_vs_ethnicity.csv")
p = ggplot(ethnicity_ancestry_counts %>% filter(!is.na(ethnicity_clean)),
aes(prop,ethnicity_clean,fill=predicted_ancestry,label = paste0("N = ",total)))+
geom_col(color="black")+
geom_text(data = ethnicity_ancestry_counts %>% distinct(ethnicity_clean,total,.keep_all=T),
mapping = aes(x=1.1))+
scale_fill_brewer(palette="Set1")+
theme_bw()+
theme(legend.position="top")+
labs(y="Self-reported ethnicity",fill="Genetic ancestry",x="Proportion")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/anc_vs_ethnicity.png",res=900,units="in",width=10,height=4)
p
dev.off()


# repeat with bigsnpr
library(bigsnpr)
library(tidyverse)

# download reference files
DIR = "/data/scratch/hmy117"
all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620968",
                         dir = DIR, fname = "ref_freqs.csv.gz"))
projection <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620953",
                         dir = DIR, fname = "projection.csv.gz"))

# match adams to reference

path_to_geno = "/data/scratch/hmy117/hgdp_1kg_genomes/adams_hg19_for_bigsnpr_cpra"
adams_snps = bigreadr::fread2(paste0(path_to_geno,".bim"),
  select = c(1, 4:6),
                   col.names = c("chr", "pos", "a1", "a0")) %>%
  mutate(beta = 1) %>%
  snp_match(all_freq[1:5]) %>%
  print()

# read matched SNPs
rds <- snp_readBed2(paste0(path_to_geno,".bed"),
 ind.col = adams_snps$`_NUM_ID_.ss`)
obj.bigsnp <- snp_attach(rds)
G <- obj.bigsnp$genotypes
G <- snp_fastImputeSimple(G)

# PCA projection
# project individuals (divided by 2) onto the PC space
PROJ <- as.matrix(projection[adams_snps$`_NUM_ID_`, -(1:5)])

correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
all_proj <- big_prodMat(G, sweep(PROJ, 2, correction / 2, '*'),
                        # scaling to get G if beta = 1 and (2 - G) if beta = -1
                        center = 1 - adams_snps$beta,
                        scale = adams_snps$beta)


X <- crossprod(PROJ,
               as.matrix(all_freq[adams_snps$`_NUM_ID_`, -(1:5)]))
               cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
Amat <- cbind(1, diag(ncol(X)))
bvec <- c(1, rep(0, ncol(X)))

# define groups
group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"
grp_fct <- factor(group, unique(group))


# assign to one group
all_centers <- t(X)
all_sq_dist <- apply(all_centers, 1, function(one_center) {
  rowSums(sweep(all_proj, 2, one_center, '-')^2)
})

THR <- 0.05  # you can adjust this threshold
thr_sq_dist <- max(dist(all_centers)^2) * THR / 0.16

cluster <- group[
  apply(all_sq_dist, 1, function(x) {
    ind <- which.min(x)
    if (isTRUE(x[ind] < thr_sq_dist)) ind else NA
  })
]

# read ancestry calls
ancestry_calls_hgdp_1kg = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls.tsv")
ancestry_calls_hgdp_1kg$bigsnpr_ancestry = cluster
write_tsv(ancestry_calls_hgdp_1kg,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_detailed.tsv")

# plot
ncol = ancestry_calls_hgdp_1kg$bigsnpr_ancestry %>% unique %>% length
pal = colorRampPalette(RColorBrewer::brewer.pal(ncol, "Set1"))(ncol)
p=ggplot(ancestry_calls_hgdp_1kg,aes(predicted_ancestry,fill=bigsnpr_ancestry))+
geom_bar(position="fill",color="black")+
theme_minimal()+
scale_fill_manual(values = pal)+
labs(y="Proportion of samples",x="Inferred genetic ancestry \nfrom HGDP-1KG reference \nusing random forest classifier",fill="Inferred genetic ancestry\nwith UKB reference\nand bigsnpr method")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_vs_bigsnpr.png",res=900,units="in",width=7,height=6)
p
dev.off()

# counts
ancestry_calls_hgdp_1kg %>%
  group_by(predicted_ancestry) %>%
  dplyr::count(bigsnpr_ancestry) %>%
  mutate(pct = 100*n/sum(n)) %>%
  print(n=100)
````

## QC of imputed data

### Explore imputation quality
- Check imputation quality
- Examine snps in R
- Plot info in MAF bins

#### Combine stats from all chroms using bcftools
````unix
cd /data/scratch/hmy117/adams_imputed/
module load bcftools
rm info_stats_all_snps
touch info_stats_all_snps
for i in {1..22};
do
	echo doing chrom $i
	bcftools query -f '%CHROM %POS %INFO/MAF %REF %ALT %INFO/R2 %INFO/ER2\n' chr$i\.info.gz >> info_stats_all_snps
done
````

#### Inpsect in R & get SNPs with high INFO & MAF
````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")
snps = read_table("info_stats_all_snps",col_names = F)
colnames(snps) = c("CHR","BP","MAF","REF", "ALT","INFO","ER2")

# cut into MAF bins
snps$maf_bin = Hmisc::cut2(snps$MAF,cuts = c(0,0.01,0.05,0.1,0.5))

# write genotyped snps for step1 regenie
snps_qcd_regenie = snps %>%
	filter(INFO >= 0.99 & MAF >= 0.05)


# find high-quality SNPs for downstream analyses
snps_qcd = snps %>%
	filter(INFO >= 0.7 & MAF >= 0.01)

# make new snp name
snps_qcd = snps_qcd %>%
	mutate(snp_name = paste0(CHR,":",BP,":",REF,":",ALT)) %>%
	dplyr::select(snp_name)
write_tsv(snps_qcd,"snps_to_keep_maf1e-2_info_0.7.tsv",col_names=F)

message("Before filtering N:")
print(nrow(snps))
message("After filtering N:")
print(nrow(snps_qcd))

snps_qcd_regenie = snps_qcd_regenie %>%
	mutate(snp_name = paste0(CHR,":",BP,":",REF,":",ALT)) %>%
	dplyr::select(snp_name)
write_tsv(snps_qcd_regenie,"snps_to_keep_step1_regenie.tsv",col_names=F)


# summarise per maf bin
summary_info = snps %>%
	group_by(maf_bin) %>%
	summarise(mean_info = mean(INFO), sd = sd(INFO))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_maf.png",res = 300, units="in",width=4, height=4)
ggplot(summary_info,aes(maf_bin,mean_info))+
	geom_point(size=3)+
	geom_errorbar(mapping = aes(x=maf_bin,ymin = mean_info - sd,ymax = mean_info + sd),width=0.3)+
	theme_minimal()+
	labs(x="Minor allele frequency",y="Imputation quality (INFO score)")+
	scale_x_discrete(labels = c("<0.01","<0.05","<0.1","<0.5"))+
	scale_fill_brewer(palette="Set1")+
	theme(legend.position = "none")
dev.off()

# genotyped snps
genotyped = snps %>% filter(ER2 !=".")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_r.png",res = 300, units="in",width=4, height=4)
ggplot(genotyped,aes(maf_bin,as.numeric(ER2),fill=maf_bin))+
geom_boxplot()+
theme_minimal()+
labs(x="Minor allele frequency",y="Imputation-Genotyping correlation (R)")+
scale_x_discrete(labels = c("<0.01","<0.05","<0.1","<0.5"))+
scale_fill_brewer(palette="Set1")+
theme(legend.position = "none")
dev.off()

p = ggplot(snps_qcd,aes(INFO,fill=maf_bin))+
geom_histogram(alpha=0.5,color="black")+
theme_minimal()+
scale_fill_brewer(palette="Set1")+
labs(x="Imputation quality score (INFO)",y="N SNPs",fill="MAF bin")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_final_snps.png",res = 300, units="in",width=6, height=4)
p
dev.off()

````

### SNP QC
- SNP QC on individual VCFs
- Filter to INFO > 0.7 & MAF > 0.01 <10% and HWE P>1e-5
- Conversion back to plink
````unix
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_snp_qc_imputed.sh

````

### Update fam files
- Modify IDs in R
````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")

for(i in c(1:22)){
df = read_table(paste0("ADAMS_imputed_qc_chr",i,".fam"),col_names=F) %>%
	mutate(X1 = str_remove(X1,"^0_")) %>%
	tidyr::separate(X1,sep = "_",into = c("part1","oragene")) %>%
	dplyr::select(X2,oragene) %>%
	mutate(oldfid = X2, oldiid = X2, newfid = oragene, newiid = oragene) %>%
	dplyr::select(oldfid,oldiid,newfid,newiid)
write_tsv(df,paste0("chr",i,"_newids.tsv"),col_names = F)
}

````

#### Update IDs in PLINK
- Merge across chromosomes in PLINK
- Rename sample IDs
````unix
# rename IDs
cd /data/scratch/hmy117/adams_imputed/
for i in {22..1}; do ~/plink --bfile ADAMS_imputed_qc_chr$i --update-ids chr$i\_newids.tsv --out chr$i\_combined_adams_imputed_newids --make-bed; done
````

### Merge chromosomes
````unix
cd /data/scratch/hmy117/adams_imputed/
rm filelist_for_merge
for i in {2..22}; do echo chr$i\_combined_adams_imputed_newids >> filelist_for_merge; done
~/plink --bfile chr1_combined_adams_imputed_newids \
--merge-list filelist_for_merge \
--out combined_adams_imputed \
--make-bed

````


### Further QC
````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")
snps = read_table("combined_adams_imputed.bim",col_names = F)

message(nrow(snps))

# exclude indels
snps = snps %>%
	filter(nchar(X5)==1 & nchar(X6)==1)

# exclude duplicate positions (i.e. non-biallelics)
dups = snps %>% dplyr::count(X1,X4)
dups =  dups %>% filter(n>1)

snps = snps %>% left_join(dups,by=c("X1","X4"))
snps = snps %>% filter(is.na(n)) %>% dplyr::select(-n)

message(nrow(snps))

write_tsv(snps,"non_duplicated_no_indel_snps.tsv",col_names=F)
````

### QC
````unix
cd /data/scratch/hmy117/adams_imputed/

~/plink --bfile combined_adams_imputed \
--extract non_duplicated_no_indel_snps.tsv \
--make-bed \
--maf 0.01 \
--geno 0.1 \
--out combined_adams_imputed_qc \
--chr 1-22
````

### Individual QC
#### Missingness
````unix
cd /data/scratch/hmy117/adams_imputed/
~/plink --bfile combined_adams_imputed_qc \
--mind 0.1 \
--make-bed \
--out combined_adams_imputed_qc_mind_0.1
````
#### Heterozygosity
````unix
cd /data/scratch/hmy117/adams_imputed/
~/plink --bfile combined_adams_imputed_qc \
--indep-pairwise 1000 500 0.1 \
--out pruned_for_het_and_kinship

~/plink --bfile combined_adams_imputed_qc \
--extract pruned_for_het_and_kinship.prune.in \
--make-bed \
--out pruned_for_het_and_kinship_genotypes

~/plink --bfile pruned_for_het_and_kinship_genotypes \
--het small-sample \
--out het_check
````

````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")
het = read_table("het_check.het")


p=ggplot(het,aes(F))+
	geom_histogram()+
	theme_minimal()+
	labs(x="Heterozygosity statistic (F)",y="N")+
	geom_vline(xintercept = mean(het$F)-5*sd(het$F),linetype="dashed",color="blue")+
	geom_vline(xintercept = mean(het$F)+5*sd(het$F),linetype="dashed",color="blue")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/het_stats.png",res=900,units="in",width=4,height=4)
p
dev.off()

````
### Kinship
#### Inference
````unix
cd /data/scratch/hmy117/adams_imputed/
~/king -b combined_adams_imputed_qc.bed --duplicate
~/king -b combined_adams_imputed_qc.bed --related --degree 3
~/plink --bfile combined_adams_imputed_qc --missing --out missingness_report
````

#### Plot missingness
````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")
miss = read_table("missingness_report.imiss")

miss$F_MISS %>% summary()

p=ggplot(miss,aes(F_MISS))+
	geom_histogram()+
	theme_minimal()+
	labs(x="Proportion of missing genotypes\n(per individual)",y="N")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/missing_stats.png",res=900,units="in",width=4,height=4)
p
dev.off()

````


# copy to home folder
````unix
cp combined_adams_imputed_qc* ~/ADAMS/genotypes/QMUL_Aug_23/outputs/
````

## Severity GWAS

### PCA
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# split ancestry groups
awk 'NR>1{if($3=="CSA") print $2,$2}' ./outputs/ancestry_calls.tsv > ./outputs/sas_iids
awk 'NR>1{if($3=="AFR") print $2,$2}' ./outputs/ancestry_calls.tsv > ./outputs/afr_iids
awk 'NR>1{if($3=="EUR") print $2,$2}' ./outputs/ancestry_calls.tsv > ./outputs/eur_iids

# recompute PCs within each cluster
for ancestry in sas eur afr
do
	~/plink2 \
	--bfile ./outputs/combined_adams_imputed_qc \
	--keep ./outputs/$ancestry\_iids \
	--maf 0.05 \
	--geno 0.01 \
	--indep-pairwise 1000 100 0.1 \
	--out ./outputs/$ancestry\_pruned

	# filter
	~/plink2 \
	--bfile ./outputs/combined_adams_imputed_qc \
	--keep ./outputs/$ancestry\_iids \
	--extract ./outputs/$ancestry\_pruned.prune.in \
	--pca 4 \
	--out ./outputs/$ancestry\_pcs
done
````

#### Add genetic sex to covar file
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
covars = read_tsv("./pheno/adams_covars.tsv")
fam = read_table("./outputs/ADAMS_geno_fid_iid_inpheno_nodups.fam",col_names = F)
covars = covars %>%
filter(IID %in% fam$X1) %>%
left_join(fam %>%
		dplyr::select(X1,X5) %>%
		dplyr::rename("IID" = X1, "genetic_sex" = X5),
		by="IID") %>%
		dplyr::select(-Sex) %>%
		dplyr::rename("Sex" = genetic_sex)
write_tsv(covars,"./pheno/adams_covars_fixed.tsv")
````

#### Make covar & pheno files
````unix
Rscript ./scripts/prepare_pheno_and_covar_files.R
````


#### Run GWAS
````unix
# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# remove PCA outliers and filter to MAC 20 within each cluster
for ancestry in sas eur afr;
do
	# filter in PLINK
	~/plink2 --bfile ./outputs/combined_adams_imputed_qc \
	--keep ./outputs/$ancestry\_iids \
	--remove ./outputs/pca_outliers_$ancestry\.tsv \
	--make-bed \
	--out ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
  --mac 20 \
	--hwe 1e-5
done
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_severity.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_severity_binary.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_severity_recessive.sh
````

#### Process results for FUMA
````unix
# lift over to hg19
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/

awk '{print "chr"$1,$4-1,$4,$2}' combined_adams_imputed_qc.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped
````

````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/")
hg19_file = read_table("hg19_bedfile",col_names=F)
for(ancestry in c("sas","afr","eur")){
  dat = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/sev_gwas_",ancestry,"_gARMSS.regenie"))
  dat = dat %>%
    dplyr::select(-GENPOS,-EXTRA,-TEST,-CHISQ) %>%
    left_join(hg19_file %>%
      dplyr::select(X4,X3) %>%
      dplyr::rename("ID"=X4),
      by="ID") %>%
      dplyr::rename("BP_hg19"=X3) %>%
      mutate(P = 10^-LOG10P) %>%
      dplyr::select(-LOG10P)
  dat = dat %>% filter(!is.na(BP_hg19))

   write_tsv(dat,paste0("/data/scratch/hmy117/",ancestry,"_sev_gwas.tsv.gz"))  
  }

````

#### For SAS - import to GH
````unix
cp ./outputs/filtered_genotypes_for_sev_gwas_sas_hg38* ./outputs/sas_for_gh
cp ./pheno/sas_* ./outputs/sas_for_gh
tar czvf ./outputs/sas_for_gh.tar.gz ./outputs/sas_for_gh
````

#### Power calc
````unix
Rscript ./scripts/power_calcs.R
````

#### Clump GWAS results
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

rm ./outputs/vep_input
for ancestry in sas eur afr
do
~/plink --bfile ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--clump /data/scratch/hmy117/$ancestry\_sev_gwas.tsv.gz \
--clump-p1 1e-5 \
--clump-r2 0.1 \
--clump-snp-field ID \
--out ./outputs/$ancestry\_clumped
done

````

#### VEP
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# print hits at P < 1e-5
rm ./outputs/vep_input
for ancestry in sas eur afr
do
	awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_gARMSS.regenie >> ./outputs/vep_input
	awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_edss.regenie >> ./outputs/vep_input
	awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_eq5d_vas.regenie >> ./outputs/vep_input
	awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_msis_physical_normalised.regenie >> ./outputs/vep_input
	awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_age_at_dx.regenie >> ./outputs/vep_input
  awk 'NR>1{if($12>5) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_edss6.regenie >> ./outputs/vep_input

done

# annotate all hits with P <1e-5
module load ensembl-vep

~/ensembl-vep/vep -i ./outputs/vep_input -o ./outputs/snp_annotations_nearest \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--check_existing \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence,Existing_variation"

````

#### Combine with IMSGC severity GWAS
##### Liftover
````unix

# lift over to hg38
cd /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS
awk 'NR>1{print "chr"$1,$3-1,$3,$2}' imsgc_mssev_discovery.tsv > hg19_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg19_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
hg38_bedfile \
unmapped
````

#### Compare freqs in IMSGC vs ADAMS-EUR
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
~/plink2 --bfile ./outputs/filtered_genotypes_for_sev_gwas_eur_hg38 \
--freq \
--out ./outputs/eur_sev_gwas_freqs
````

##### Liftover
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS")
imsgc_sev = read_tsv("imsgc_mssev_discovery.tsv")
hg38_positions = read_table("hg38_bedfile",col_names=F)

# prepare for merge
hg38_positions = hg38_positions %>%
			dplyr::select(X3,X4) %>%
			dplyr::rename("SNP" = X4,"BP_hg38" = X3)

# find duplicate positions
dups = hg38_positions %>%
	dplyr::count(SNP) %>%
	filter(n>1)

# filter to biallelics
hg38_positions = hg38_positions %>%
	filter(!SNP %in% dups$SNP)

# filter
imsgc_sev = imsgc_sev %>%
	filter(SNP %in% hg38_positions$SNP)


# merge
imsgc_sev = imsgc_sev %>%
	left_join(hg38_positions,by="SNP")

# find compatible SNPs with ADAMS
adams_geno = read_table("../QMUL_Aug_23/outputs/combined_adams_imputed_qc.bim",col_names=F)

# make chr:pos in hg38
imsgc_sev = imsgc_sev %>%
	dplyr::select(-POS) %>%
	dplyr::rename("POS" = BP_hg38) %>%
	mutate(chrpos = paste0("chr",CHR,":",POS))

# chrpos for adams
adams_geno = adams_geno %>%
	mutate(chrpos = paste0("chr",X1,":",X4))

# find intersecting positions
intersecting_positions = imsgc_sev %>%
	inner_join(adams_geno,by="chrpos")

# find compatible alleles
intersecting_positions = intersecting_positions %>%
	filter(  (A1 == X5 & A2 == X6) | (A1 == X6 & A2 == X5)   )

# exclude palindromes
intersecting_positions = intersecting_positions %>%
filter(
	!(A1 == "C" & A2 == "G" ) &
	!(A1 == "G" & A2 == "C" ) &
	!(A1 == "A" & A2 == "T" ) &
	!(A1 == "T" & A2 == "A" )
	)
# compare with ADAMS-EUR freqs
adams_eur_freqs = read_table("../QMUL_Aug_23/outputs/eur_sev_gwas_freqs.afreq") %>%
	separate(ID,sep=":",into=c("chr","pos","a1","a2")) %>%
	mutate(chrpos = paste0(chr,":",pos))

intersecting_positions = intersecting_positions %>%
left_join(adams_eur_freqs %>% dplyr::select(-c(1:5)),
by="chrpos")

# find delta freq
intersecting_positions = intersecting_positions %>%
	mutate(delta_freq = ifelse(
		A1 == ALT, abs(AF1 - ALT_FREQS), abs( (1 - AF1) - ALT_FREQS)
		))

# remove SNPs with AF diff >10%
intersecting_positions = intersecting_positions %>%
	filter(delta_freq < 0.1)
# write to file
write_tsv(intersecting_positions %>%
	dplyr::select(c(1:11)),"imsgc_mssev_discovery_hg38.tsv")
````


#### PRS
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# PRSice
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
awk 'NR==1{print "SNP","A1","A2","P","BETA"};NR>1{print "chr"$1":"$10,$3,$4,$8,$6}' /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_mssev_discovery_hg38.tsv > snps_for_severity_prs_all_snps

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/prs_severity.sh
````

#### Explore PRS
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

# read summary of all prs at r2 of 0.1
all_prs = data.frame()
for(ancestry in c("sas","eur","afr")){
	prs = read_table(paste0("./outputs/prs_0.1_",ancestry,".prsice")) %>%
	mutate(anc=ancestry)
	all_prs <<- bind_rows(all_prs,prs)
}

# plot
p = ggplot(all_prs,
  aes(factor(Threshold),Coefficient,fill=toupper(anc),label=paste0("SNPs: ",Num_SNP,"\nP: ",round(P,2))))+
  geom_errorbar(mapping = aes(x = factor(Threshold),
  ymin = Coefficient - 1.96 * `Standard.Error`,
  ymax = Coefficient + 1.96 * `Standard.Error`),
  width=0.1)+
  geom_point(size=3,shape=21,color="black")+
  facet_wrap(~toupper(anc),nrow=3)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  geom_text(mapping = aes(y = Coefficient + 1.96 * `Standard.Error` + 0.1),size=3)+
  geom_hline(yintercept=0,linetype="dashed",color="pink")+
  scale_y_continuous(limits = c(-0.3,0.5))+
  theme(legend.position="none")+
  labs(x="P value threshold for PRS",y="Beta coefficient\n(per-SD effect of PRS on gARMSS)")


png("./outputs/prs_plots_summ.png",res=900,units="in",width=10,height=8)
p
dev.off()

# read in all prs
all_prs = data.frame()
for(r2 in c(0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","eur","afr")){
	prs = read_table(paste0("./outputs/prs_",r2,"_",ancestry,".prsice")) %>%
	mutate(anc=ancestry,clump_r2 = r2)
	all_prs <<- bind_rows(all_prs,prs)
}
}

# plot
p = ggplot(all_prs,
  aes(factor(Threshold),-log10(P),fill=toupper(clump_r2)))+
  geom_col(color="black",position=position_dodge())+
  facet_wrap(~toupper(anc),nrow=3)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
  labs(x="P value threshold for PRS",
  y=bquote(-log[10]~P),
  fill=bquote(Clumping~R^2))
png("./outputs/prs_plots_all_clumps_summ.png",res=900,units="in",width=8,height=6)
p
dev.off()

# read in scores for all phenos
all_prs = data.frame()
for(ancestry in c("sas","eur","afr")){
	prs = read_table(paste0("./outputs/prs_all_phenos_",ancestry,".prsice")) %>%
	mutate(anc=ancestry)
	all_prs <<- bind_rows(all_prs,prs)
}

# fix pheno names
all_prs$Pheno =  str_remove_all(all_prs$Pheno,"rint_")

all_prs = all_prs %>%
mutate(Pheno = case_when(
  Pheno == "gARMSS" ~ "gARMSS",
  Pheno == "edss" ~ "EDSS",
  Pheno == "msis_physical_normalised" ~ "MSIS",
  Pheno == "eq5d_vas" ~ "EQ5D",
  Pheno == "age_at_dx" ~ "Age at diagnosis"))

# plot
p = ggplot(all_prs,
  aes(factor(Threshold),-log10(P),fill=Pheno))+
  geom_col(color="black",position=position_dodge())+
  facet_wrap(~toupper(anc),nrow=3)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  geom_hline(yintercept=-log10(0.05/5),linetype="dashed",color="pink")+
  labs(x="P value threshold for PRS",
  y=bquote(-log[10]~P),
  fill="Phenotype")
png("./outputs/prs_plots_all_phenos_summ.png",res=900,units="in",width=8,height=6)
p
dev.off()

# look at coeffs
# plot
p = ggplot(all_prs,
  aes(factor(Threshold),Coefficient,fill=toupper(anc)))+
  geom_errorbar(mapping = aes(x = factor(Threshold),
  ymin = Coefficient - 1.96 * `Standard.Error`,
  ymax = Coefficient + 1.96 * `Standard.Error`),
  width=0.1)+
  geom_point(size=3,shape=21,color="black")+
  facet_grid(toupper(anc)~Pheno)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  geom_hline(yintercept=0,linetype="dashed",color="pink")+
  theme(legend.position="none",axis.text.x = element_text(angle=90))+
  labs(x="P value threshold for PRS",y="Beta coefficient\n(per-SD effect of PRS on phenotype)")

png("./outputs/prs_plots_all_phenos_coefs.png",res=900,units="in",width=10,height=8)
p
dev.off()


````

#### Explore GWAS results
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
# read in GWAS
read_regenie_gwas = function(x){
  read_table(x, col_types = cols_only(
		  CHROM = col_double(),
		  GENPOS = col_double(),
		  ID = col_character(),
		  ALLELE0 = col_character(),
		  ALLELE1 = col_character(),
		  A1FREQ = col_double(),
		  N = col_double(),
		  TEST = col_character(),
		  BETA = col_double(),
		  SE = col_double(),
		  CHISQ = col_double(),
		  LOG10P = col_double())) %>%
  filter(TEST == "ADD" & !is.na(LOG10P)) %>%
  dplyr::rename("CHR" = CHROM, "BP" = GENPOS, "SNP" = ID) %>%
  dplyr::select(-TEST) %>%
	mutate(P = 10^-LOG10P)
}

# define ancestries
ancestries=c("afr","sas","eur")

# define phenotypes
phenos = c("gARMSS","edss","edss6","msis_physical_normalised","eq5d_vas","age_at_dx")

# read in all gwas
all_res = list()
for(pheno in phenos){
	for(ancestry in ancestries){

		gwas_res = read_regenie_gwas(
			paste0("./outputs/sev_gwas_",ancestry,"_",pheno,".regenie")) %>%
		mutate(anc = ancestry,phenotype = pheno)
		all_res[[length(all_res)+1]] = gwas_res
	}
}
# combine
all_res = do.call("bind_rows",all_res)

saveRDS(all_res,"/data/scratch/hmy117/all_gwas_res_ms_sev.rds")
````

#### QQ plots and manhattan plots
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
all_res = readRDS("/data/scratch/hmy117/all_gwas_res_ms_sev.rds")

# get N
counts = all_res %>%
	group_by(anc,phenotype) %>%
	summarise(n = median(N)) %>%
	ungroup() %>%
	pivot_wider(id_cols = phenotype,values_from=n,names_from=anc)
write_csv(counts,"./outputs/sev_gwas_pheno_counts.csv")

# read annotations
anno = read_table("./outputs/snp_annotations_nearest",skip=31,col_types = cols(.default="c")) %>%
	dplyr::select(1,NEAREST,Consequence,Existing_variation)
colnames(anno) = c("SNP","Gene","Effect","rsid")
anno = anno %>% distinct(SNP,Gene,rsid,.keep_all=T)


# define ancestries
ancestries=c("afr","sas","eur")

# define phenotypes
phenos = c("gARMSS","edss","edss6","msis_physical_normalised","eq5d_vas","age_at_dx")
phenolabs = c("gARMSS","EDSS","EDSS6","MSIS","EQ5D","Age at diagnosis")

# inflation
qqplots = list()
make_qq_plot = function(ancestry,pheno){
	gwas_res = all_res %>% filter(anc == ancestry & phenotype==pheno)
	lambda = median(qchisq(gwas_res$P,df=1,lower.tail=F)) /  qchisq(0.5,lower.tail=F,df=1)
	message(ancestry)
	message(lambda)
	pval_dat = data.frame(observed = -log10(gwas_res$P)) %>%
		arrange(observed)
	simulated_pvals = sort(runif(n = nrow(pval_dat)),decreasing=T)
	pval_dat = pval_dat %>% mutate(exp_pval = simulated_pvals)
	pval_dat = pval_dat %>% mutate(expected = -log10(exp_pval))

 # sample
 pval_dat = pval_dat %>% sample_n(size=500000)
	# plot
	p = ggplot(pval_dat,aes(expected,observed))+
		geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+
		geom_point()+
		labs(x="Expected -log10(P)",y="Observed -log10(P)")+
		ggtitle(paste0("Ancestry: ",toupper(ancestry),"\nLambda =  ",round(lambda,2)))+
		theme_bw()

	qqplots[[length(qqplots)+1]] <<- p
	png(paste0("./outputs/qq_plot_",ancestry,"_",pheno,".png"),res=900,units="in",width=4,height=4)
	print(p)
	dev.off()

}

for(pheno in phenos){
	for(ancestry in ancestries){
		make_qq_plot(ancestry,pheno)
	}
}

# print all together
library(gridExtra)
png("./outputs/all_qq_plots_garmss.png",res=900,units="in",width=8,height=4)
cowplot::plot_grid(plotlist = qqplots[c(1:3)],align="v",ncol=3)
dev.off()


# manhattans
make_manhattans = function(ancestry,pheno,phenolab){
	gwas_res = all_res %>% filter(anc == ancestry, phenotype == pheno)

	# get coords
	coords = gwas_res %>%
		group_by(CHR) %>%
		summarise(min_bp = min(BP),max_bp=max(BP),median_bp = median(BP))

	new_coords = coords %>%
		mutate(CHR = CHR+1) %>%
		mutate(cumbp = cumsum(max_bp))

	midpoints = coords$median_bp+c(0,new_coords$cumbp[1:21])
	gwas_res = gwas_res %>%
			left_join(new_coords,by="CHR") %>%
			mutate(cum_bp = ifelse(is.na(cumbp),BP,BP+cumbp))

	# colors
	gwas_res = gwas_res %>%
		mutate(significance = case_when(
			P < 1e-5 ~ "sig",
			CHR %% 2 == 0 ~ "even",
			CHR %% 2 != 0 ~ "odd"
			))

	col_pal = c("lavenderblush1","lavenderblush2","orange")
	names(col_pal) = c("even","odd","sig")

	plot_dat = gwas_res %>% filter(P < 0.05)
	p=ggplot(plot_dat,aes(cum_bp,-log10(P),col=significance))+
		geom_hline(yintercept = 5,color="blue",linetype="dashed",alpha=0.5)+
		geom_point()+
		scale_color_manual(values = col_pal)+
		scale_x_continuous(breaks = midpoints,labels = c(1:22))+
		ggrepel::geom_label_repel(plot_dat %>%
			filter(P < 1e-5 & SNP %in% anno$SNP) %>%
			left_join(anno,by="SNP") %>%
			distinct(Gene,.keep_all=T),
		mapping = aes(label = Gene),min.segment.length = 0,color="black",nudge_y=1,direction = "x")+
		labs(x="Genomic position (hg38)",y="-log10(P)")+
		ggtitle(paste0("GWAS of ",phenolab," in ",toupper(ancestry)))+
		theme_bw()+
		theme(legend.position="none")

	manhattans[[length(manhattans)+1]] <<- p

	png(paste0("./outputs/manhattan_plot_",pheno,"_",ancestry,".png"),res=900,units="in",width=10,height=4)
	print(p)
	dev.off()
}

# do all phenos
for(i in c(length(phenos):1)){  
manhattans = list()
for(ancestry in ancestries){
  message(ancestry)
	make_manhattans(ancestry,pheno=phenos[i],phenolab = phenolabs[i])
}
# print all together
png(paste0("./outputs/all_manhattans_",phenos[i],"_plots.png"),res=900,units="in",width=8,height=10)
print(cowplot::plot_grid(plotlist = manhattans,align="v",ncol=1))
dev.off()
}




# compare with IMSGC GWAS
imsgc_hg38 = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_mssev_discovery_hg38.tsv")

# read in clump results
all_clump_res = list()
for(ancestry in ancestries){
  all_clump_res[[length(all_clump_res)+1]] = read_table(paste0("./outputs/",ancestry,"_clumped.clumped")) %>%
  mutate(anc = ancestry)  %>%
  dplyr::select(SNP,anc) %>%
  mutate(topsnp = "topsnp")
}
all_clump_res = do.call("bind_rows",all_clump_res)

# print garmss sig hits
garmss_sig = all_res %>%
	filter(phenotype=="gARMSS" & P < 1e-5) %>%
  left_join(all_clump_res,by=c("SNP","anc")) %>%
	left_join(anno,by="SNP") %>%
  mutate(chrpos = paste0("chr",CHR,":",BP)) %>%
  left_join(imsgc_hg38 %>% dplyr::rename("rsid" = SNP),
  by=c("chrpos")) %>%
  filter(topsnp == "topsnp")

write_csv(garmss_sig,"./outputs/garmss_sig_res.csv")

# print msis sig hits

# print garmss sig hits
msis_sig = all_res %>%
	filter(phenotype=="msis_physical_normalised" & P < 1e-5) %>%
  left_join(all_clump_res,by=c("SNP","anc")) %>%
	left_join(anno,by="SNP") %>%
  mutate(chrpos = paste0("chr",CHR,":",BP)) %>%
  left_join(imsgc_hg38 %>% dplyr::rename("rsid" = SNP),
  by=c("chrpos"))
write_csv(msis_sig,"./outputs/msis_sig_res.csv")

age_at_dx = all_res %>%
	filter(phenotype=="age_at_dx" & P < 1e-5) %>%
  left_join(all_clump_res,by=c("SNP","anc")) %>%
	left_join(anno,by="SNP") %>%
  mutate(chrpos = paste0("chr",CHR,":",BP)) %>%
  left_join(imsgc_hg38 %>% dplyr::rename("rsid" = SNP),
  by=c("chrpos"))
write_csv(age_at_dx,"./outputs/age_at_dx_res.csv")

# combine garmss GWAS and imsgc join on chr:pos
all_res_with_imsgc = imsgc_hg38 %>%
		dplyr::rename("BP"=POS,"ALLELE1" = A1,"ALLELE0" = A2,"A1FREQ" = AF1) %>%
		mutate(chrpos = paste0(CHR,":",BP)) %>%
		mutate(anc = "IMSGC") %>%
		bind_rows(all_res %>%
      filter(phenotype=="gARMSS") %>%
      mutate(chrpos = paste0(CHR,":",BP)))

# save to pick up here
saveRDS(all_res_with_imsgc,"/data/scratch/hmy117/all_sev_res_gwas_with_imsgc.rds",compress=F)
# all_res_with_imsgc = readRDS("/data/scratch/hmy117/all_sev_res_gwas_with_imsgc.rds")

# search sig hits for nearby suggestive hits
imsgc_hg38 %>%
filter(P < 1e-5)

# forest
make_forest = function(snp){
  plot_dat = all_res_with_imsgc %>% filter(chrpos == snp)

  # Add study label
  plot_dat = plot_dat %>%
  	mutate(anc = ifelse(anc == "IMSGC","IMSGC-EUR",paste0("ADAMS-",toupper(anc))))

  # flip beta if necessary
  plot_dat = plot_dat %>%
  	mutate(BETA = ifelse(ALLELE1 == plot_dat$ALLELE1[1],BETA,BETA*-1))

  print(plot_dat)

  # refactor
  plot_dat$anc = factor(plot_dat$anc,levels = c("IMSGC-EUR","ADAMS-EUR","ADAMS-SAS","ADAMS-AFR"),ordered=T)
  pictured_allele = paste0(snp,"-",plot_dat$ALLELE1[1])
  ggplot(plot_dat,aes(BETA,anc,fill=anc))+
  	geom_errorbarh(mapping = aes(xmin = BETA - 1.96*SE,xmax = BETA + 1.96*SE,y=anc),height=0.3)+
  	geom_point(shape=21,color="black",size=3)+
  	geom_vline(xintercept=0,linetype="dashed")+
  	theme_bw()+
  	labs(x=paste0("Effect of ",pictured_allele," on gARMSS"),y="Study & ancestry")+
  	theme(legend.position="none")
}

png("./outputs/dysf_snp.png",res=900,units="in",width=6,height=4)
	make_forest("2:71449869")
dev.off()

png("./outputs/armc_forest_snp.png",res=900,units="in",width=6,height=4)
	make_forest("10:22886937")
dev.off()

png("./outputs/nrxn3_forest_snp.png",res=900,units="in",width=6,height=4)
	make_forest("14:79410873")
dev.off()

png("./outputs/trem_forest_snp.png",res=900,units="in",width=6,height=4)
	make_forest("6:41320217")
dev.off()


# abbreviated list of IMSGC sig hits (from ST9 of Nature GWAS)
top_hits = c("rs10191329",
  "rs149097173",
  "rs2876767",
  "rs181310516",
  "rs147933117",
  "rs194722",
  "rs12494504",
  "rs112663015",
  "rs9397000",
  "rs61215450",
  "rs4251626",
  "rs115687581")

# get sig snps in IMSGC
sig_imsgc_variants = all_res_with_imsgc %>%
  filter(anc=="IMSGC") %>%
  filter(SNP %in% top_hits)

# loop through all IMSGC sig hits
plots = list()
for(i in c(1:nrow(sig_imsgc_variants))){
  plots[[i]] = make_forest(sig_imsgc_variants$chrpos[i])
}

png("./outputs/all_imsgc_snps.png",res=900,units="in",width=8,height=10)
cowplot::plot_grid(plotlist = plots,ncol = 3)
dev.off()

png("./outputs/dysf_snps.png",res=900,units="in",width=8,height=6)
cowplot::plot_grid(
	make_forest("2:71449869"),
	make_forest("2:71450816"),
	make_forest("2:71452575"),
	align="h",
	nrow=3
	)
dev.off()


# region plot
chr=6
pos=41320217
window=5e5
gene_name="TREM1"
make_locus_plot = function(chr,pos,window=1e6,gene_name="TREM1"){
  plot_dat = all_res_with_imsgc %>% filter(CHR == chr & BP > pos - window/2 & BP < pos + window/2)


	# calculate LD
	for(study in unique(plot_dat$anc) ){
		# get top SNP for this ancestry

		snp = plot_dat %>% filter(anc == study)  %>%
		slice_min(P)
    ancestry = ifelse(study=="IMSGC","eur",snp$anc[1])
		this_snp = ifelse(study=="IMSGC",
    paste0("chr",snp$chrpos[1],":",snp$ALLELE0[1],":",snp$ALLELE1[1]),
    snp$SNP[1])

		message(this_snp)

		cmd = paste0("~/plink --bfile ./outputs/filtered_genotypes_for_sev_gwas_",ancestry,"_hg38 ",
							 "--ld-snp ",this_snp," --r2 --ld-window-r2 0 --ld-window 9999999 --ld-window-kb 999999 --out ./outputs/ld_",study)
	system(cmd)
	}

	# read in LD
	ld = purrr::map(unique(plot_dat$anc),function(x){
		y = paste0("./outputs/ld_",x,".ld")
		read_table(y) %>%
			mutate(anc = x)
		})
	ld = do.call("bind_rows",ld) %>%
    mutate(top_snp = ifelse(SNP_A == SNP_B,"top_snp"," "))
	ld = ld %>%
		mutate(chrpos = paste0(CHR_B,":",BP_B)) %>%
		dplyr::select(chrpos,anc,R2,top_snp)

	# Add study label
  plot_dat = plot_dat %>%
		mutate(ancestry = ifelse(anc=="IMSGC","eur",anc)) %>%
  	mutate(anc = ifelse(anc == "IMSGC","IMSGC-EUR",paste0("ADAMS-",toupper(anc)))) 	
  ld = ld %>%
		mutate(ancestry = ifelse(anc=="IMSGC","eur",anc)) %>%
  	mutate(anc = ifelse(anc == "IMSGC","IMSGC-EUR",paste0("ADAMS-",toupper(anc)))) 	

	# add ld
	plot_dat = plot_dat %>%
    left_join(ld %>% dplyr::select(-ancestry),
    by=c("chrpos","anc"))


  plot_dat_wide = plot_dat %>%
    dplyr::select(CHR,BP,ALLELE0,ALLELE1,A1FREQ,BETA,SE,P,R2,top_snp,anc) %>%
    pivot_wider(id_cols = c(CHR,BP),
    names_from = anc,
    values_from = c(A1FREQ,BETA,SE,P,R2,top_snp,ALLELE0,ALLELE1))

  write_csv(plot_dat_wide,paste0("./outputs/ld_locus_",gene_name,".csv"))
	pal = ggsci::pal_locuszoom("default")(7)

  # refactor
  plot_dat$anc = factor(plot_dat$anc,levels = c("IMSGC-EUR","ADAMS-EUR","ADAMS-SAS","ADAMS-AFR"),ordered=T)

  plot_dat %>%
    group_by(anc) %>%
    slice_min(P) %>%
    dplyr::select(SNP,A1FREQ,BETA,SE,P,anc) %>%
    print()
	plot_dat$r2_bin = Hmisc::cut2(plot_dat$R2,cuts = c(0,0.2,0.4,0.6,0.8,1))
	names(pal) = rev(levels(plot_dat$r2_bin))

	ggplot(plot_dat,aes(BP,-log10(P),fill=r2_bin))+
  	geom_point(data = plot_dat %>% filter(top_snp != "top_snp"), shape=21,color="black",size=3)+
    geom_point(data = plot_dat %>% filter(top_snp == "top_snp"), shape=23,color="black",fill="purple",size=3)+
  	geom_hline(yintercept=-log10(1e-5),linetype="dashed",color="pink")+
  	theme_bw()+
    facet_wrap(~anc,ncol=1)+
		scale_fill_manual(values = (pal))+
		labs(fill = bquote(R^2),y=bquote(-log[10]~P))

}

png("./outputs/ncr3_locus_plot.png",res=900,units="in",width=6,height=8)
make_locus_plot(
chr=6,
pos=41320217,
window=5e5)
dev.off()

png("./outputs/nrxn_locus_plot.png",res=900,units="in",width=6,height=8)
make_locus_plot(
chr=14,
pos=79410873,
gene_name="NRXN3",
window=5e5)
dev.off()

png("./outputs/armc3_locus_plot.png",res=900,units="in",width=6,height=8)
make_locus_plot(
chr=10,
pos=22886937,
gene_name="ARMC3",
window=5e5)
dev.off()


png("./outputs/dysf_locus_plot.png",res=900,units="in",width=6,height=8)
make_locus_plot(chr=2,pos=71449869,window=1e6)
dev.off()

png("./outputs/dnm3_locus_plot.png",res=900,units="in",width=6,height=8)
make_locus_plot(chr=1,pos=200258565,window=5e5)
dev.off()

# global beta beta plots
sig_imsgc_hits = all_res_with_imsgc %>% filter(anc=="IMSGC" & P < 1e-5) %>% dplyr::select(chrpos,SNP,P)
all_res_sig_imsgc_hits = all_res_with_imsgc %>% filter(SNP %in% sig_imsgc_hits$SNP)


## BEN PICK UP HERE
# look at cam severity
cam_cov = read_tsv("/data/home/hmy117/ADAMS/genotypes/Cambridge/sas_covars_with_pcs.tsv")
cam_pheno = read_tsv("/data/home/hmy117/ADAMS/genotypes/Cambridge/sas_pheno.tsv")
cam_all = cam_cov %>%
  left_join(cam_pheno,by=c("FID","IID")) %>%
  filter(study=="CAM") %>%
  filter(!is.na(gARMSS))

chr=6
pos=41320217
A1 = "T"
A2 = "C"
snp_to_test = paste0("chr",chr,":",pos,":",A1,":",A2)
cmd = paste0("~/plink2 --bfile /data/home/hmy117/ADAMS/genotypes/Cambridge/filtered_genotypes_for_sev_gwas_sas_hg38 ",
"--recode AD --snp ",snp_to_test," ",
"--out ./outputs/snp_genos_cam")
system(cmd)

# read in
genos = read_table("./outputs/snp_genos_cam.raw")

# join
flip_alleles=F
dat = genos %>%
  filter(IID %in% cam_all$IID) %>%
  left_join(cam_all,by=c("IID","FID"))

# get colname
snp_colname = colnames(dat)[7]

dat$snp_geno_raw = dat[[snp_colname]]
dat$snp_geno = if(flip_alleles==T){
  2-as.numeric(dat$snp_geno_raw)
  } else {
  as.numeric(dat$snp_geno_raw)
  }
print(dat$snp_geno)
phenotype="gARMSS"  
dat = dat %>% filter(!is.na(snp_geno))
dat$phenocol = dat[[phenotype]]
counts = dat %>%
  dplyr::count(snp_geno)
medians = dat %>%
  group_by(snp_geno) %>%
  summarise(median = median(phenocol,na.rm=T))
label="gARMSS"
ylim=10
    p = ggplot(dat,
      aes(factor(snp_geno),
      phenocol,fill=factor(snp_geno)))+
      geom_boxplot(width=0.1,position=position_dodge(width=0.7),fatten=3)+
      geom_violin(alpha=0.1,position=position_dodge(width=0.7),show.legend=F)+
      geom_dotplot(binaxis="y",stackdir="center",dotsize = 0.5,position=position_dodge(width=0.7))+
      theme_bw()+
      scale_fill_brewer(palette="Set1")+
      geom_text(data = counts,aes(y = ylim*1.1,label=paste0("N=",n)),position=position_dodge(width=0.7))+
      geom_text(data = medians,aes(y = ylim*1.05,label=paste0("Median=",round(median,1))),position=position_dodge(width=0.7))+
      labs(y=label,fill="Genotype")+
      theme(legend.position="none")+
      labs(x=snp_colname)

    png(
      paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/",str_replace_all(snp_to_test,":","_"),"pheno_anc_",
      phenotype,
      ".png"),
      res=900,units="in",width=10,height=4)
    print(p)
    dev.off()


summary(lm(data = cam_dat,
  gARMSS ~ `chr6:41320217:T:C_T`) )


````

#### Impact of EUR variants - ben finish 08-02
### recessive
### fliup alleles
### other snps
### other phenotypes
````R
library(tidyverse)
setwd("~/ADAMS/genotypes/QMUL_Aug_23/")
pheno = read_tsv("./pheno/adams_pheno.tsv")
cov = read_tsv("./pheno/adams_covars_fixed.tsv")
anc = read_tsv("./outputs/ancestry_calls.tsv")

make_pheno_plot = function(phenotype,label,chr,pos,A1,A2,flip_alleles=T,ylim=10){

  snp_to_test = paste0("chr",chr,":",pos,":",A1,":",A2)
  cmd = paste0("~/plink2 --bfile ./outputs/combined_adams_imputed_qc ",
  "--recode AD --snp ",snp_to_test," ",
  "--out ./outputs/snp_genos")
  system(cmd)

  # read in
  genos = read_table("./outputs/snp_genos.raw")

  # join
  dat = genos %>%
    filter(IID %in% anc$IID & IID %in% cov$IID & IID %in% pheno$IID) %>%
  	left_join(pheno,by=c("IID")) %>%
  	left_join(cov,by=c("IID")) %>%
  	left_join(anc,by="IID")

  # get colname
  snp_colname = colnames(dat)[7]

  dat$snp_geno_raw = dat[[snp_colname]]
  dat$snp_geno = if(flip_alleles==T){
    2-as.numeric(dat$snp_geno_raw)
    } else {
    as.numeric(dat$snp_geno_raw)
    }
  print(dat$snp_geno)


	dat$phenocol = dat[[phenotype]]
	counts = dat %>%
	  filter(predicted_ancestry %in% c("CSA","AFR","EUR")) %>%
	  dplyr::count(predicted_ancestry,snp_geno)
	medians = dat %>%
	  filter(predicted_ancestry %in% c("CSA","AFR","EUR")) %>%
		group_by(predicted_ancestry,snp_geno) %>%
	  summarise(median = median(phenocol,na.rm=T))


  	p = ggplot(dat %>% filter(predicted_ancestry %in% c("CSA","AFR","EUR")),
  		aes(factor(snp_geno),
      phenocol,fill=predicted_ancestry))+
  		facet_wrap(~predicted_ancestry)+
  		geom_boxplot(width=0.1,position=position_dodge(width=0.7),fatten=3)+
  		geom_violin(alpha=0.1,position=position_dodge(width=0.7),show.legend=F)+
  		geom_dotplot(binaxis="y",stackdir="center",dotsize = 0.5,position=position_dodge(width=0.7))+
  	  theme_bw()+
  	  scale_fill_brewer(palette="Set1")+
  	  geom_text(data = counts,aes(y = ylim*1.1,label=paste0("N=",n)),position=position_dodge(width=0.7))+
  		geom_text(data = medians,aes(y = ylim*1.05,label=paste0("Median=",round(median,1))),position=position_dodge(width=0.7))+
  	  labs(y=label,fill="Genetic ancestry")+
  		theme(legend.position="none")+
      labs(x=snp_colname)

  	png(
  		paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/",str_replace_all(snp_to_test,":","_"),"pheno_anc_",
  		phenotype,
  		".png"),
  		res=900,units="in",width=10,height=4)
  	print(p)
  	dev.off()
}

make_pheno_plot("gARMSS","Global ARMSS",
chr=9,
pos=26257387,
A1 = "A",
A2 = "T")

make_pheno_plot("gARMSS","Global ARMSS",
chr=2,
pos=44473990,
A1 = "G",
A2 = "A")

make_pheno_plot("age_at_dx","Age at diagnosis",
chr=2,
pos=63933730,
ylim=90,
A1="G",
A2="C")



make_pheno_plot("gARMSS","Global ARMSS",
chr=6,
pos=41320217,
A1 = "T",
A2 = "C")
make_pheno_plot("edss","EDSS")
make_pheno_plot("msis_physical_normalised","MSIS-Physical",ylim=100)

make_pheno_plot("gARMSS","Global ARMSS",
chr=14,
pos=79410873,
A1 = "G",
A2 = "T")

make_pheno_plot("gARMSS","Global ARMSS",
chr=10,
pos=22886937,
A1 = "A",
A2 = "G")



# plot
ggplot(dat %>%
	filter(predicted_ancestry %in% c("CSA","AFR","EUR")),
	aes(factor(rs10191329_A),gARMSS,fill=predicted_ancestry))+
	geom_boxplot()

````

#### HLA & Age at dx
````unix
module load bcftools
bcftools view -i 'ID="HLA_DRB1*15:01"' /data/scratch/hmy117/adams_imputed/hla/chr6.dose.vcf.gz -Oz > /data/scratch/hmy117/adams_imputed/hla/drb1_15_01.vcf.gz
~/plink --vcf /data/scratch/hmy117/adams_imputed/hla/drb1_15_01.vcf.gz \
--double-id \
--recode A include-alt \
--out ~/ADAMS/genotypes/QMUL_Aug_23/hla_drb1_genos
````

````R

library(tidyverse)
setwd("~/ADAMS/genotypes/QMUL_Aug_23/")
pheno = read_tsv("./pheno/adams_pheno.tsv")
cov = read_tsv("./pheno/adams_covars_fixed.tsv")
anc = read_tsv("./outputs/ancestry_calls.tsv")



# read in
genos = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/hla_drb1_genos.raw") %>%
separate(IID,sep="_",into = c("x","y","IID")) %>%
dplyr::select(IID,`HLA_DRB1*15:01_T(/A)`) %>%
mutate(DRB1_1501_alleles = `HLA_DRB1*15:01_T(/A)`)

# join
dat = genos %>%
  mutate(IID = as.numeric(IID)) %>%
  filter(IID %in% anc$IID & IID %in% cov$IID & IID %in% pheno$IID) %>%
	left_join(pheno,by=c("IID")) %>%
	left_join(cov,by=c("IID")) %>%
	left_join(anc,by="IID")


counts = dat %>%
  filter(predicted_ancestry %in% c("CSA","AFR","EUR")) %>%
  dplyr::count(predicted_ancestry,DRB1_1501_alleles)
medians = dat %>%
  filter(predicted_ancestry %in% c("CSA","AFR","EUR")) %>%
	group_by(predicted_ancestry,DRB1_1501_alleles) %>%
  summarise(median = median(age_at_dx,na.rm=T))

ylim = 90

p = ggplot(dat %>% filter(predicted_ancestry %in% c("CSA","AFR","EUR")),
aes(factor(DRB1_1501_alleles),
age_at_dx,fill=predicted_ancestry))+
facet_wrap(~predicted_ancestry)+
geom_boxplot(width=0.1,position=position_dodge(width=0.7),fatten=3)+
geom_violin(alpha=0.1,position=position_dodge(width=0.7),show.legend=F)+
geom_dotplot(binaxis="y",stackdir="center",dotsize = 0.5,position=position_dodge(width=0.7))+
theme_bw()+
scale_fill_brewer(palette="Set1")+
geom_text(data = counts,aes(y = ylim*1.1,label=paste0("N=",n)),position=position_dodge(width=0.7))+
geom_text(data = medians,aes(y = ylim*1.05,label=paste0("Median=",round(median,1))),position=position_dodge(width=0.7))+
labs(y="Age at diagnosis",fill="Genetic ancestry")+
theme(legend.position="none")+
labs(x="HLA-DRB1*15:01 alleles")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/age_at_dx_drb1.png",
res=900,units="in",width=10,height=4)
print(p)
dev.off()

````

## Susceptibility GWAS

### CAM QC - run on CAM HPC
````unix
cd /rds/project/sjs1016/rds-sjs1016-msgen/Genotypes/Montreal/Imputed/
module load plink/2.00-alpha

# snp qc
for i in {1..22};
  do
    plink2 --pfile MS_GT_only_chr$i \
    --maf 0.05 \
    --hwe 1e-5 \
    --geno 0.1 \
    --make-bed \
    --out ~/rds/hpc-work/filtered_chr$i\_plink1 &
  done

# merge chroms
rm merge_filelist
for i in {2..22};
  do
    echo ~/rds/hpc-work/filtered_chr$i\_plink1.bed ~/rds/hpc-work/filtered_chr$i\_plink1.bim ~/rds/hpc-work/filtered_chr$i\_plink1.fam >> merge_filelist
  done

module unload plink
module load plink-1.9-gcc-5.4.0-sm3ojoi

plink --bfile ~/rds/hpc-work/filtered_chr1_plink1 \
--merge-list merge_filelist \
--make-bed \
--out merged_plink1_filtered

# filter merged file
plink --bfile merged_plink1_filtered \
--mind 0.1 \
--make-bed \
--out ./ancestry/outputs/merged_qcd_all_chroms


````

### Liftover to hg19 (for UKB)
````unix

# lift over to hg19
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/

awk '{print "chr"$1,$4-1,$4,$2}' combined_adams_imputed_qc.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions
~/plink --bfile combined_adams_imputed_qc_unrelated \
--update-map hg19_snp_positions \
--make-bed \
--out combined_adams_imputed_qc_temp_hg19

# reset var names to cpra
~/plink2 --bfile combined_adams_imputed_qc_temp_hg19 \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--out combined_adams_imputed_qc_hg19_cpra
````

### Merge with UKB & CAM
#### Find compatible SNPs
````R
library(tidyverse)

# read ukb snps
ukb_snps = list()
for(i in c(1:22)){
	ukb_snps[[i]] = read_table(paste0("/data/Wolfson-PNU-dementia/UKB/imputed_genotypes/plink2_files/chr_",i,".pvar"))
}
ukb_snps = do.call("bind_rows",ukb_snps)

# read adams snps
adams_snps  = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc_hg19_cpra.bim",col_names = F)
colnames(adams_snps) = c("#CHROM","ID","X3","POS","REF","ALT")

# find intersection of chr:pos
adams_snps = adams_snps %>% mutate(chrpos = paste0(`#CHROM`,":",POS))
ukb_snps = ukb_snps %>% mutate(chrpos = paste0(`#CHROM`,":",POS))

combo_snps = inner_join(adams_snps,ukb_snps,by="chrpos")

# clean up columns
combo_snps = combo_snps %>% dplyr::select(-`#CHROM.y`,-X3,-POS.y)

# check compatible alleles
combo_snps = combo_snps %>%
	mutate(compatible = ifelse(
		(REF.x == REF.y & ALT.x == ALT.y) |
		(REF.x == ALT.y & ALT.x == REF.y),
		"compatible","incompatible"))
combo_snps = combo_snps %>%
		filter(compatible == "compatible")


# read Cam SNPs
cam_snps = read_table("/data/scratch/hmy117/cambridge_genos/ancestry_inference/merged_qcd_all_chroms.bim",col_names=F)
cam_snps = cam_snps %>% mutate(chrpos = paste0(X1,":",X4))

# find intersection
combo_snps_with_cam = inner_join(cam_snps,combo_snps,by="chrpos")

# find compatible SNPs
combo_snps_with_cam = combo_snps_with_cam %>%
	mutate(compatible_with_cam = ifelse(
		(X5 == REF.x & X6 == ALT.x) |
		(X5 == ALT.x & X6 == REF.x),
		"compatible","incompatible")) %>%
		filter(compatible_with_cam == "compatible")

# write snp list to keep
adams_snps_to_keep = combo_snps_with_cam %>% dplyr::select(ID.x) %>% distinct()
ukb_snps_to_keep = combo_snps_with_cam %>% dplyr::select(ID.y) %>% distinct()
ukb_snps_names_to_update = combo_snps_with_cam %>% dplyr::select(ID.y,ID.x) %>% distinct(ID.y,.keep_all=T)
cam_snps_to_keep = combo_snps_with_cam %>% dplyr::select(X2) %>% distinct()
cam_snps_names_to_update = combo_snps_with_cam %>% dplyr::select(X2,ID.x) %>% distinct(X2,.keep_all=T)


write_tsv(adams_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_snps_to_keep.tsv",col_names = F)
write_tsv(ukb_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_keep.tsv",col_names = F)
write_tsv(ukb_snps_names_to_update,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_update_names.tsv",col_names = F)
write_tsv(cam_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/cam_snps_to_keep.tsv",col_names = F)
write_tsv(cam_snps_names_to_update,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/cam_snps_to_update_names.tsv",col_names = F)

````

#### Filter to compatible SNPs
##### UKB
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_ukb_files.sh
````
##### CAM
````unix

# set wd
cd /data/scratch/hmy117

# Filter Cam SNPs to ADAMS SNPs
~/plink2 \
--bfile /data/scratch/hmy117/cambridge_genos/ancestry_inference/merged_qcd_all_chroms \
--extract /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/cam_snps_to_keep.tsv \
--make-bed \
--rm-dup exclude-all \
--out /data/scratch/hmy117/cam_filtered_tmp

# update names
~/plink2 \
--bfile  /data/scratch/hmy117/cam_filtered_tmp \
--update-name /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/cam_snps_to_update_names.tsv \
--out /data/scratch/hmy117/cam_filtered_hg19_cpra \
--make-bed
````
##### ADAMS
````unix
~/plink2 \
--bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc_hg19_cpra \
--extract /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_snps_to_keep.tsv \
--out /data/scratch/hmy117/adams_filtered_hg19_cpra \
--make-bed

````

#### Merge Cam files + ADAMS genotypes
````unix
~/plink --bfile /data/scratch/hmy117/cam_filtered_hg19_cpra \
--bmerge /data/scratch/hmy117/adams_filtered_hg19_cpra \
--out /data/scratch/hmy117/adams_cam_hg19_merge \
--make-bed
````

#### Merge Cam-ADAMS with UKB files across chromosomes
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_ukb_files.sh
````


#### QC per chromosome
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merged_ukb_cam_adams_qc_per_chrom.sh
````

#### Clean up
````unix
# clean up
for i in {1..22}; do rm /data/scratch/hmy117/ukb_cam_adams_merged_genotypes_chr$i\.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/ukb_adams_merged_genotypes_chr$i\.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\_hg19_cpra.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\_hg19_cpra_nodups.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/ukb_chr$i\_filtered_tmp.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\.* ; done
````

#### Prepare HGDP-KG for joint PCA on combined dataset
````unix

cd /data/scratch/hmy117/hgdp_1kg_genomes/

# liftover
awk '{print "chr"$1,$4-1,$4,$2}' combined_hgdp_1kg_unrelated.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions
~/plink --bfile combined_hgdp_1kg_unrelated \
--update-map hg19_snp_positions \
--make-bed \
--out combined_hgdp_1kg_unrelated_hg19

# reset var names to cpra
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19 \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--out combined_hgdp_1kg_unrelated_hg19_cpra

# get list of all ADAMS-CAM-UKB vars
rm adams_cam_ukb_snps_to_keep
for i in {1..22};
do
	awk '{print $2}' /data/scratch/hmy117/ukb_cam_adams_merged_genotypes_chr$i\_no_palindromes.bim >> adams_cam_ukb_snps_to_keep
done

# filter to UKB-CAM-ADAMS vars
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra \
--extract adams_cam_ukb_snps_to_keep \
--out combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry \
--make-bed
````


#### Prune each chromosome (for PCA)
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merged_ukb_cam_adams_prune_for_pca.sh
````

#### Merge chromosomes for PCA
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_all_chroms_adams_ukb.sh
````

#### Do PCA
````unix
cd /data/scratch/hmy117/hgdp_1kg_genomes/

# prune
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry \
--extract /data/scratch/hmy117/ukb_cam_adams_merged_genotypes_all_chroms_for_pca.bim \
--indep-pairwise 1000 500 0.05 \
--out combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_prune

~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry \
--out combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_pruned \
--extract combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_prune.prune.in \
--make-bed

# calculate PCs
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_pruned \
--pca allele-wts 50 \
--freq \
--out hgdp_kg_pcs_with_ukb

# project ADAMS-UKB-CAM samples
~/plink2 --bfile /data/scratch/hmy117/ukb_cam_adams_merged_genotypes_all_chroms_for_pca \
--read-freq hgdp_kg_pcs_with_ukb.afreq \
--score hgdp_kg_pcs_with_ukb.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize list-variants \
--score-col-nums 6-55 \
--out ukb_cam_adams_pcs

# project original dataset
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_pruned \
--read-freq hgdp_kg_pcs_with_ukb.afreq \
--extract ukb_cam_adams_pcs.sscore.vars \
--score hgdp_kg_pcs_with_ukb.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize \
--score-col-nums 6-55 \
--out hgdp_kg_pcs_with_ukb_rescored
````

#### Identify ancestry groupings
````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")

# read in data
adams_ukb = read_table("ukb_cam_adams_pcs.sscore",col_types = cols(.default = "c"))

kg_hgdp = read_table("hgdp_kg_pcs_with_ukb_rescored.sscore")
meta = read_tsv("gnomad_meta_v1.tsv") %>%
  dplyr::select(2,hgdp_tgp_meta.Population,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("IID","pop","superpop")
kg_hgdp = kg_hgdp %>%
left_join(meta,by="IID")

# build RF on HGDP data
# filter out those with missing labels
kg_hgdp = kg_hgdp %>%
  filter(!is.na(superpop)) %>%
  dplyr::select(superpop,contains("PC"))
library(caret)

rf_fit = train(superpop ~ .,
                      data=kg_hgdp,
                      method='rf',
                      metric='Accuracy')


# save rf
saveRDS(rf_fit,"rf_fit_with_ukb.rds")

# can reload here
# rf_fit = readRDS("rf_fit_with_ukb.rds")

# predict
adams_ukb$predicted_ancestry = predict(rf_fit,adams_ukb)

# validate these ancestry calls
# read in original adams ancestry calls
ancestry_calls_hgdp_1kg = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_detailed.tsv") %>%
	dplyr::select(IID,contains("ancestry"))
cov = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv") %>% dplyr::select(-FID)
all_cov = ancestry_calls_hgdp_1kg %>% left_join(cov,by="IID")
adams_fam = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc_hg19_cpra.fam",col_names=F)

# read in cam pheno
cam_pheno = read_csv("/data/home/hmy117/ADAMS/genotypes/Cambridge/ClinicalData_V1.csv")

# read ukb pheno
ukb_pheno = read_csv("/data/Wolfson-PNU-dementia/UKB/PRS_April_2023/phenodata/PRS_59138r672130.csv") %>%
	dplyr::select(EID,sex_f31_0_0,age_at_recruitment_f21022_0_0,genetic_ethnic_grouping_f22006_0_0,contains("g35"))

# bring in ethnicity data
ukb_pheno_ethnicity = readRDS("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/78867r672482/78867r672482_FO.rds") %>%
	tibble %>%
	dplyr::select(eid,contains("ethnic"))
id_bridge = read_csv("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/59138_78867/Bridge_eids_59138_78867.csv")

ukb_pheno_ethnicity = ukb_pheno_ethnicity %>%
  dplyr::rename("eid_78867" = eid) %>%
  left_join(id_bridge,by="eid_78867") %>%
  dplyr::rename("EID" = eid_59138)
ukb_pheno_ethnicity = ukb_pheno_ethnicity %>% dplyr::select(EID,contains("ethnic"))
ukb_pheno = ukb_pheno %>% left_join(ukb_pheno_ethnicity,by="EID")

# define cohort
adams_ukb = adams_ukb %>%
	mutate(cohort = case_when(
		IID %in% ukb_pheno$EID ~ "UKB",
		IID %in% cam_pheno$SampleID ~ "CAM",
		IID %in% adams_fam$X2 ~ "ADAMS"
		))

# define case-control status
cam_cases = cam_pheno %>% filter(AffectionStatus == 2)
ukb_cases = ukb_pheno %>% filter(!is.na(source_of_report_of_g35_multiple_sclerosis_f131043_0_0))

adams_ukb = adams_ukb %>%
	mutate(MS_status = case_when(
		cohort == "UKB" & IID %in% ukb_cases$EID ~ "MS",
		cohort == "CAM" & IID %in% cam_cases$SampleID ~ "MS",
		cohort == "ADAMS" ~ "MS",
		.default = as.character("Control")
		))


# counts
counts = adams_ukb %>%
filter(!is.na(cohort)) %>%
	group_by(cohort) %>%
	dplyr::count(predicted_ancestry,MS_status) %>%
	mutate(pct = n/sum(n)*100) %>%
	mutate(n_pct = paste0(n," (",round(pct,0),"%)")) %>%
	ungroup() %>%
	pivot_wider(id_cols = c(cohort,MS_status),
	names_from = predicted_ancestry,
	values_from = c(n_pct))
write_csv(counts,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_calls_tbl.csv")

# plot
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_calls_prop.png",res=900,units="in",width=6,height=4)
ggplot(adams_ukb %>% filter(!is.na(cohort)),
	aes(cohort,fill=predicted_ancestry))+
	geom_bar(position="fill",color="black")+
	theme_bw()+
	scale_fill_brewer(palette="Set1")+
	labs(x="Cohort",y="Proportion",fill="Genetic ancestry")
dev.off()


# join with UKB pheno
adams_ukb = adams_ukb %>%
	dplyr::rename("EID" = IID) %>%
	left_join(ukb_pheno %>% mutate(EID = as.character(EID)),by="EID")

# coerce PCs back to numbers
adams_ukb = adams_ukb %>%
mutate_at(.vars = vars(contains("^PC")),as.numeric)

# plot overall PC space
p1 = ggplot(adams_ukb %>% filter(!is.na(cohort)),aes(PC1_AVG,PC2_AVG,fill=predicted_ancestry))+
  geom_point(shape=21,color="black",size=3)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
	facet_wrap(~cohort) +
  labs(x="PC1",y="PC2",fill="Genetic ancestry")

p2 = ggplot(adams_ukb %>% filter(!is.na(cohort)),aes(PC3_AVG,PC4_AVG,fill=predicted_ancestry))+
  geom_point(shape=21,color="black",size=3)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
	facet_wrap(~cohort) +
  labs(x="PC3",y="PC4",fill="Genetic ancestry")

p3 = ggplot(adams_ukb %>% filter(!is.na(cohort)),aes(PC5_AVG,PC6_AVG,fill=predicted_ancestry))+
  geom_point(shape=21,color="black",size=3)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
	facet_wrap(~cohort) +
  labs(x="PC5",y="PC6",fill="Genetic ancestry")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_adams_pc_plot_with_ref_pc1_6.png",res=900,units="in",width=8,height=12)
cowplot::plot_grid(p1,p2,p3,align="v",ncol=1)
dev.off()

# validate with self-reported ethnicity
adams_ukb_with_ethnicity = adams_ukb %>%
	filter(!is.na(adams_ukb$ethnic_background_f21000_0_0))%>%
	filter(predicted_ancestry %in% c("AFR","EUR","CSA"))
ethnicities_to_plot = adams_ukb_with_ethnicity %>%
	dplyr::count(ethnic_background_f21000_0_0) %>%
	filter(n>200)

p = ggplot(adams_ukb_with_ethnicity %>% filter(ethnic_background_f21000_0_0 %in% ethnicities_to_plot$ethnic_background_f21000_0_0),
aes(ethnic_background_f21000_0_0,fill=predicted_ancestry))+
  geom_bar(position="fill",color="black")+
  theme_bw()+
	coord_flip()+
	scale_fill_brewer(palette="Set1")+
	labs(y="Proportion of individuals",x="Self-reported ethnicity",fill="Inferred genetic ancestry")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_adams_pcs_vs_self_report.png",res=900,units="in",width=8,height=6)
p
dev.off()

# repeat with adams ethnicities
adams_with_ethnicity = adams_ukb %>%
	mutate("IID" = as.numeric(EID)) %>%
	left_join(all_cov,by="IID") %>%
	filter(!is.na(ethnicity_clean))%>%
	filter(predicted_ancestry.x %in% c("AFR","EUR","CSA"))

caret::confusionMatrix(factor(adams_with_ethnicity$predicted_ancestry.x),
factor(adams_with_ethnicity$predicted_ancestry.y))

p = ggplot(adams_with_ethnicity,
aes(ethnicity_clean,fill=predicted_ancestry.x))+
  geom_bar(position="fill",color="black")+
  theme_bw()+
	coord_flip()+
	scale_fill_brewer(palette="Set1")+
	labs(y="Proportion of individuals",x="Self-reported ethnicity",fill="Inferred genetic ancestry")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_pcs_vs_self_report_risk.png",res=900,units="in",width=8,height=4)
p
dev.off()

# repeat with cam ethnicities
cam_with_ethnicity = adams_ukb %>%
	mutate("SampleID" = as.character(EID)) %>%
	left_join(cam_pheno,by="SampleID") %>%
	filter(!is.na(ethnicity))%>%
	filter(predicted_ancestry %in% c("AFR","EUR","CSA"))

p = ggplot(cam_with_ethnicity %>%
	filter(ethnicity != "X"),
aes(ethnicity,fill=predicted_ancestry))+
  geom_bar(position="fill",color="black")+
  theme_bw()+
	coord_flip()+
	scale_fill_brewer(palette="Set1")+
	labs(y="Proportion of individuals",x="Self-reported ethnicity",fill="Inferred genetic ancestry")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/cam_pcs_vs_self_report_risk.png",res=900,units="in",width=8,height=4)
p
dev.off()

# write all ancestry calls
adams_ukb = adams_ukb %>%
	filter(!is.na(cohort))

anc_calls = adams_ukb %>%
	dplyr::rename("IID" = EID) %>%
  dplyr::select(`#FID`,IID,predicted_ancestry)
write_tsv(anc_calls,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_all_ukb_adams.tsv")

# write pheno - covar
pheno_cov = adams_ukb %>%
	dplyr::rename("IID" = EID) %>%
	dplyr::select(-contains("PC")) %>%
	dplyr::select(IID,MS_status,sex_f31_0_0,age_at_recruitment_f21022_0_0,cohort) %>%
	dplyr::rename("age" = age_at_recruitment_f21022_0_0, "sex" = sex_f31_0_0) %>%
	mutate(sex = ifelse(sex == "Male",1,2))

# process cam pheno
cam_pheno_filtered = cam_pheno %>%
	dplyr::select(SampleID,sex,YOB,YO_edss) %>%
	mutate(yob = as.numeric(YOB),yo_edss = as.numeric(YO_edss)) %>%
	mutate(age = ifelse(is.na(yo_edss),2010 - yob,yo_edss - yob)) %>%
	dplyr::select(SampleID,age,sex) %>%
	dplyr::rename("IID" = SampleID) %>%
	na.omit()

# process adams pheno
adams_cov = cov %>%
	mutate(IID = as.character(IID)) %>%
 	dplyr::select(IID,ageatedss,Sex) %>%
	mutate(sex = ifelse(Sex == "Male",1,2)) %>%
	dplyr::rename("age" = ageatedss) %>%
	dplyr::select(IID,age,sex) %>%
	na.omit()

cov_adams_cam = bind_rows(cam_pheno_filtered,adams_cov)

non_ukb = pheno_cov %>%
		filter(cohort != "UKB") %>%
		dplyr::select(IID,cohort,MS_status) %>%
		left_join(cov_adams_cam,by="IID")
ukb = pheno_cov %>%
		filter(cohort == "UKB")
combo_dat = bind_rows(non_ukb,ukb)

# remove NAs
combo_dat = combo_dat %>% na.omit()
write_tsv(combo_dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_cam_adams_pheno_cov.tsv")

# demographics
combo_dat = combo_dat %>%
	left_join(anc_calls,by="IID")

# make basic demographics table
library(compareGroups)
combo_dat = combo_dat %>%
	mutate(gender = ifelse(
		sex == 1,"Male","Female"
		)) %>%
	mutate("Genetic_ancestry" = ifelse(
		predicted_ancestry %in% c("MID","OCE","AMR"),"Other",predicted_ancestry
		))
tbl = compareGroups(data = combo_dat,
	MS_status ~ cohort + age + gender + Genetic_ancestry,
	method = c(3,2,3,3)
	) %>%
	createTable()

export2csv(tbl,"./outputs/risk_gwas_demographics.csv")

````

#### Check plausible AF differences
````R
library(tidyverse)

or = 3
beta = log(or)
baseline_risk = 0.002
baseline_log_odds = log(baseline_risk / (1-baseline_risk) )
het_log_odds = baseline_log_odds + beta
hom_log_odds = baseline_log_odds + 2*beta   

# get disease probs per genotype
convert_log_odds_to_prob = function(x){
	exp(x) * (1-baseline_risk)
}
het_prob = convert_log_odds_to_prob(het_log_odds)
hom_prob = convert_log_odds_to_prob(hom_log_odds)

# simulate genotypes with different mafs
mafs = seq(0.01,0.5,by=0.01)
delta_afs = list()
for(i in c(1:length(mafs))){

	maf = mafs[i]

	geno_df = data.frame(genos = rbinom(n=1e7,size=2,prob = maf))

	# for each geno, simulate disease status
	common_hom = geno_df %>%
		filter(genos == 0 )
	common_hom$dis_status = rbinom(n = nrow(common_hom),size=1,prob = baseline_risk)
	het = geno_df %>%
		filter(genos == 1 )
	het$dis_status = rbinom(n = nrow(het),size=1,prob = het_prob)

	rare_hom = geno_df %>%
		filter(genos == 2 )
	rare_hom$dis_status = rbinom(n = nrow(rare_hom),size=1,prob = hom_prob)

	all_dat = bind_rows(common_hom,het,rare_hom)

	# find mafs
	all_dat = all_dat %>%
		group_by(dis_status) %>%
		dplyr::count(genos) %>%
		mutate(total_alleles = sum(n)*2) %>%
		mutate(total_rare_alleles = sum(genos*n)) %>%
		mutate(maf = total_rare_alleles / total_alleles) %>%
		distinct(dis_status,maf) %>%
		pivot_wider(names_from=dis_status,values_from=maf)

	delta_af = abs( all_dat$`0` - all_dat$`1`)
	delta_afs[[i]] = delta_af

}

df = data.frame(mafs,delta_af = unlist(delta_afs))

p = ggplot(df,aes(mafs,delta_af))+
	theme_bw()+
	geom_line()+
	geom_point(size=3,color="black",shape=21,fill="white")+
	labs(x="MAF",y="Delta allele frequency")+
	scale_y_continuous(limits = c(0,0.3),breaks=c(0,0.1,0.2,0.3))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/maf_vs_delta_af.png",res=900,units="in",width=4,height=4)
p
dev.off()

````


#### Identify ancestry groupings (part 2)

##### Split into ancestry groups
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/split_ancestry_groups.sh
````

##### Prune, combine across chroms per ancestry, and run PCA per ancestry
###### Filter out PCA outliers within each ancestry
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/run_pca_per_ancestry_risk.sh
````

##### Standardise PCs in covariate files
````R
library(tidyverse)

# function to z score
z_score = function(x){
	y = ( x - mean(x) ) / sd(x)
	y
}

# standardize age and PCs 1 - 4
for(ancestry in c("sas","afr","eur")){
	in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_",ancestry,"_covars_with_pcs.tsv")
	pcs = read_tsv(in_file) %>%
		mutate_at(.vars = c("age","PC1","PC2","PC3","PC4"),z_score)
	write_tsv(pcs,in_file)
}
````

#### GWAS prep
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_susceptibility_prep.sh
````

#### GWAS
#### Run GWAS of ADAMS EUR vs UKB EUR
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
awk '{if($3==2) print $1,$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_eur_pheno.tsv > /data/scratch/hmy117/ms_cases_only

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_eur \
--keep /data/scratch/hmy117/ms_cases_only \
--remove /data/scratch/hmy117/cam_filtered_hg19_cpra.fam \
--make-pheno /data/scratch/hmy117/adams_filtered_hg19_cpra.fam \* \
--make-bed \
--allow-no-sex \
--out /data/scratch/hmy117/ms_cases_only_just_ukb_adams

~/plink --bfile /data/scratch/hmy117/ms_cases_only_just_ukb_adams \
--assoc \
--extract /data/scratch/hmy117/snps_to_keep_for_gwaseur.tsv \
--allow-no-sex \
--out /data/scratch/hmy117/ms_cases_ukb_vs_adams

awk '{if($9<0.05 && !($1==6 && $3>25e6 && $3<35e6)) print $2}' /data/scratch/hmy117/ms_cases_ukb_vs_adams.assoc > /data/scratch/hmy117/ms_cases_ukb_vs_adams_snps_to_exclude

````

#### Run GWAS

````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_susceptibility.sh
````

#### Cojo
````unix

for ancestry_out in sas afr;
do
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
# make ma file
awk 'NR==1{print "SNP","A1","A2","freq","b","se","p","N"};NR>1{print $3,$5,$4,$6,$13,$14,$9}' ./outputs/case_control_gwas_$ancestry_out\_MS_status.regenie > ./outputs/case_control_gwas_cojo_$ancestry_out  

~/gcta --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry_out \
--extract /data/scratch/hmy117/snps_to_keep_for_gwas$ancestry_out\.tsv \
--cojo-file ./outputs/case_control_gwas_cojo_$ancestry_out \
--cojo-slct \
--cojo-p 1e-5 \
--chr 6 \
--extract-region-bp 6 25000000 10000 \
--out ./outputs/case_control_gwas_cojo_cond_$ancestry_out
done

````
##### Explore in R
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")


make_conditioned_mhc_plot = function(ancestry){
  gwas = read_table(paste0("./outputs/case_control_gwas_",ancestry,"_MS_status.regenie"))
  dat = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_cojo_cond_",ancestry,".cma.cojo"))
  dat = dat %>% filter(bp > 25e6 & bp < 35e6)
  gwas = gwas %>% filter(CHROM==6 & GENPOS > 25e6 & GENPOS < 35e6)
  cond_snps = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_cojo_cond_",ancestry,".jma.cojo"))
  gwas = gwas %>% mutate(cond_snp = ifelse(ID %in% cond_snps$SNP,"yes","no"))

  p1 = ggplot(gwas,aes(GENPOS,LOG10P))+
  geom_point(data = gwas %>% filter(LOG10P>5),fill="orange",shape=21)+
  geom_point(data = gwas %>% filter(LOG10P<5),fill="grey",alpha=0.2,shape=21)+  
  theme_bw()+
  ggrepel::geom_text_repel(data = gwas %>% filter(cond_snp=="yes"),mapping=aes(label=ID),min.segment.length=0,direction="x",nudge_y=1)+
  geom_hline(yintercept=-log10(1e-5),linetype="dashed",color="pink")+
  ggtitle("Non-conditioned")+
  labs(x="Position on chromosome 6 (hg19)",y=bquote(-log[10]~P))
  p2 = ggplot(dat,aes(bp,-log10(pC)))+
  geom_point(data = dat %>% filter(pC<1e-5),fill="orange",shape=21)+
  geom_point(data = dat %>% filter(pC>1e-5),fill="grey",alpha=0.2,shape=21)+  
  theme_bw()+
  geom_hline(yintercept=-log10(1e-5),linetype="dashed",color="pink")+
  ggtitle("Conditioned")+
  labs(x="Position on chromosome 6 (hg19)",y=bquote(-log[10]~P))

  outfile=paste0("./outputs/cond_mhc_plot_",ancestry,".png")
  png(outfile,res=900,height=8,width=6,units="in")
  print(cowplot::plot_grid(p1,p2,align="v",ncol=1))
  dev.off()
}

make_conditioned_mhc_plot("afr")
make_conditioned_mhc_plot("sas")

````


#### HLA imputation
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hibag_imputation.sh
````
##### Process HLA results
````R

library(tidyverse)


overall_res_table = list()
overall_allele_freqs = list()
make_hla_forest = function(ancestry,allele,model="additive"){
  message(ancestry)
  message(allele)

  in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_allele_",allele,"_",ancestry,".tsv")

  allele_calls = read_tsv(in_file,col_types="cccdd")
  colnames(allele_calls) = c("IID","allele1","allele2","prob","match")
  allele_counts = allele_calls %>%
    dplyr::select(c(1:3)) %>%
    pivot_longer(cols = c(2,3)) %>%
    dplyr::select(-name) %>%
    group_by(IID) %>%
    dplyr::count(value)

    # join with phenot & covars

  pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_",ancestry,"_pheno.tsv"),col_types = cols(.default="c"))
  cov = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_",ancestry,"_covars_with_pcs.tsv"),col_types = "ccdddddddddddd")

  allele_counts = allele_counts %>%
    ungroup() %>%
    left_join(pheno,by="IID") %>%
    left_join(cov,by="IID")

  # get overall AFs
  all_allele_freqs = allele_counts %>%
    group_by(value) %>%
    summarise(ac = sum(n)) %>%
    ungroup() %>%
    mutate(total = sum(ac)) %>%
    mutate(af = ac/total) %>%
    arrange(af)
  all_allele_freqs$gene = allele
  all_allele_freqs$anc = ancestry  
  overall_allele_freqs[[length(overall_allele_freqs)+1]] <<- all_allele_freqs

  common_alleles = all_allele_freqs %>% filter(af > 0.05)


  # regress vs MS status for each allele
  alleles = unique(common_alleles$value)

  all_res = list()
  for(i in c(1:length(alleles))){
    message(i)
    this_allele = alleles[i]

    # filter
    this_allele_counts = allele_counts %>% filter(value == this_allele)  

    # get people with 0 counts
    this_allele_counts_zeroes = allele_counts %>% filter(!IID %in% this_allele_counts$IID)  %>%
    distinct(IID,.keep_all=T) %>%
    mutate(n = 0)

    all_dat = this_allele_counts %>%
      bind_rows(this_allele_counts_zeroes)

    all_dat$MS_status = relevel(factor(all_dat$MS_status),ref="1")


    # dominant & recessive coding
    all_dat = all_dat %>%
        mutate(n_rec = ifelse(n <=1 ,0, 1)) %>%
        mutate(n_dom = ifelse(n >=1 ,1, 0))

    # make all models and compare
      glm(data=all_dat,factor(MS_status) ~ age + sex + PC1 + PC2 + PC3 + PC4 + n + n_rec + n_dom, family=binomial(link="logit"))

    hla_model = if(model == "additive"){
      glm(data=all_dat,factor(MS_status) ~ age + sex + PC1 + PC2 + PC3 + PC4 + n, family=binomial(link="logit"))
    } else if(model=="recessive"){
      glm(data=all_dat, factor(MS_status) ~ age + sex + PC1 + PC2 + PC3 + PC4 + n_rec, family=binomial(link="logit"))
    } else if(model=="dominant"){
      glm(data=all_dat, factor(MS_status) ~ age + sex + PC1 + PC2 + PC3 + PC4 + n_dom, family=binomial(link="logit"))
    }


    if(nrow(summary(hla_model)$coefficients)==8){
      all_res[[i]] = summary(hla_model)$coefficients[8,]
    } else {
      na_vec = c(NA,NA,NA,NA)    
      names(na_vec) = colnames(summary(hla_model)$coefficients)
      all_res[[i]] = na_vec
    }
  }
  all_res = do.call("bind_rows",all_res)
  colnames(all_res) = c("beta","se","z","p")
  all_res$hla_allele = alleles
  all_res$gene = allele
  all_res$anc = ancestry
  overall_res_table[[length(overall_res_table)+1]] <<- all_res

  all_res = all_res %>% arrange(desc(beta))

  all_res$hla_allele = factor(all_res$hla_allele,ordered=T,levels = all_res$hla_allele)
  outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_allele_",allele,"_",ancestry,".png")
  p = ggplot(all_res,aes(beta,hla_allele))+
    geom_errorbarh(mapping = aes(xmin = beta - 1.96*se,xmax = beta +1.96*se,y=hla_allele),height=0.1)+
    geom_point(size=3,shape=21,fill="orange")+
    geom_vline(xintercept=0,linetype="dashed",color="pink")+
    labs(y=paste0(allele," allele"),x="Beta (log Odds Ratio)")+
    theme_bw()+
    ggtitle(paste0(toupper(ancestry),"\n",allele))


  png(outfile,res=900,units="in",width=4,height=4)
  print(p)
  dev.off()
}  

# make parameter table
alleles = c("A","B","C","DQB1","DQA1","DPA1","DPB1","DRB1","DRB3","DRB4","DRB5")
ancestries = c("sas","afr")
param_tbl = expand.grid(alleles,ancestries)

## ADDITIVE
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i])
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",hla_allele))

bonf = 0.05/11
overall_res_table$direction = ifelse(overall_res_table$beta>0,"up","down")
p =  ggplot(overall_res_table,aes(full_allele,-log10(p),fill=toupper(anc),shape=direction))+
    geom_point(data = overall_res_table %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(bonf),linetype="dashed",color="pink")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="Ancestry",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table %>% filter(p<0.05),mapping = aes(label = full_allele))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all.png",res=900,units="in",width=12,height=6)
p
dev.off()

# add freqs
overall_allele_freqs = do.call("bind_rows",overall_allele_freqs)
overall_allele_freqs = overall_allele_freqs %>%
 mutate(full_allele = paste0(gene,"*",value)) %>%
 dplyr::select(anc,full_allele,af)

overall_res_table_additive = overall_res_table %>%
  mutate(OR = exp(beta), lower_ci = exp(beta - 1.96*se),upper_ci = exp(beta+1.96*se)) %>%
  dplyr::select(gene,full_allele,anc,OR,lower_ci,upper_ci,p)  %>%
  left_join(overall_allele_freqs,by=c("full_allele","anc")) %>%
  pivot_wider(id_cols = c(1:2),names_from=anc,values_from=c(OR,lower_ci,upper_ci,p,af)) %>%
  dplyr::select(gene,full_allele,contains("sas"),contains("afr"))

## RECESSIVE

overall_res_table = list()
overall_allele_freqs = list()
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i],model="recessive")
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",hla_allele))

bonf = 0.05/11
overall_res_table$direction = ifelse(overall_res_table$beta>0,"up","down")
p =  ggplot(overall_res_table,aes(full_allele,-log10(p),fill=toupper(anc),shape=direction))+
    geom_point(data = overall_res_table %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(bonf),linetype="dashed",color="pink")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="Ancestry",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table %>% filter(p<0.05),mapping = aes(label = full_allele))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all_recessive.png",res=900,units="in",width=12,height=6)
p
dev.off()

# save
overall_res_table_recessive = overall_res_table

## DOMINANT

overall_res_table = list()
overall_allele_freqs = list()
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i],model="dominant")
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",hla_allele))

bonf = 0.05/11
overall_res_table$direction = ifelse(overall_res_table$beta>0,"up","down")
p =  ggplot(overall_res_table,aes(full_allele,-log10(p),fill=toupper(anc),shape=direction))+
    geom_point(data = overall_res_table %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(bonf),linetype="dashed",color="pink")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="Ancestry",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table %>% filter(p<0.05),mapping = aes(label = full_allele))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all_dominant.png",res=900,units="in",width=12,height=6)
p
dev.off()

overall_res_table_dominant = overall_res_table

# combine
overall_res_table_combined = overall_res_table_additive %>%
left_join(
 overall_res_table_dominant %>% dplyr::select(full_allele,anc,p) %>% dplyr::rename("p_dom" = p) %>% pivot_wider(id_cols = full_allele,names_from=anc,names_prefix="p_dom_",values_from=p_dom),
 by="full_allele"
 ) %>%
 left_join(
overall_res_table_recessive %>% dplyr::select(full_allele,anc,p) %>% dplyr::rename("p_rec" = p) %>% pivot_wider(names_prefix="p_rec_",id_cols = full_allele,names_from=anc,values_from=p_rec),
by="full_allele"
)

overall_res_table_combined = overall_res_table_combined %>%
  rowwise() %>%
  mutate(min_p_afr = min(p_afr,p_dom_afr,p_rec_afr)) %>%
  ungroup() %>%
  mutate(afr_best_model = case_when(
    p_afr == min_p_afr ~ "additive",
    p_dom_afr == min_p_afr ~ "dominant",
    p_rec_afr == min_p_afr ~ "recessive"  
    )) %>%
  rowwise() %>%
  mutate(min_p_sas = min(p_sas,p_dom_sas,p_rec_sas)) %>%
  ungroup() %>%
  mutate(sas_best_model = case_when(
    p_sas == min_p_sas ~ "additive",
    p_dom_sas == min_p_sas ~ "dominant",
    p_rec_sas == min_p_sas ~ "recessive"  
    )) %>%
    dplyr::select(gene,full_allele,contains("sas"),contains("afr"))

write_csv(overall_res_table_combined,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all.csv")

# repeat with 2 digit resolution

overall_res_table = list()
overall_allele_freqs = list()
make_hla_forest_twodigit = function(ancestry,allele){
  message(ancestry)
  message(allele)

  in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_allele_",allele,"_",ancestry,".tsv")
  # truncate to two-digit
  allele_calls = read_tsv(in_file,col_types="cccdd")
  allele_calls$allele1 = substr(allele_calls$allele1,0,2)
  allele_calls$allele2 = substr(allele_calls$allele2,0,2)



  colnames(allele_calls) = c("IID","allele1","allele2","prob","match")
  allele_counts = allele_calls %>%
    dplyr::select(c(1:3)) %>%
    pivot_longer(cols = c(2,3)) %>%
    dplyr::select(-name) %>%
    group_by(IID) %>%
    dplyr::count(value)

    # join with phenot & covars

  pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_",ancestry,"_pheno.tsv"),col_types = cols(.default="c"))
  cov = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_",ancestry,"_covars_with_pcs.tsv"),col_types = "ccdddddddddddd")

  allele_counts = allele_counts %>%
    ungroup() %>%
    left_join(pheno,by="IID") %>%
    left_join(cov,by="IID")

  # get overall AFs
  all_allele_freqs = allele_counts %>%
    group_by(value) %>%
    summarise(ac = sum(n)) %>%
    ungroup() %>%
    mutate(total = sum(ac)) %>%
    mutate(af = ac/total) %>%
    arrange(af)
  all_allele_freqs$gene = allele
  all_allele_freqs$anc = ancestry  
  overall_allele_freqs[[length(overall_allele_freqs)+1]] <<- all_allele_freqs

  common_alleles = all_allele_freqs %>% filter(af > 0.05)


  # regress vs MS status for each allele
  alleles = unique(common_alleles$value)

  all_res = list()
  for(i in c(1:length(alleles))){

  this_allele = alleles[i]

  # filter
  this_allele_counts = allele_counts %>% filter(value == this_allele)  

  # get people with 0 counts
  this_allele_counts_zeroes = allele_counts %>% filter(!IID %in% this_allele_counts$IID)  %>%
  distinct(IID,.keep_all=T) %>%
  mutate(n = 0)

  all_dat = this_allele_counts %>%
    bind_rows(this_allele_counts_zeroes)

  all_dat$MS_status = relevel(factor(all_dat$MS_status),ref="1")

  model = glm(data=all_dat,factor(MS_status) ~ age + sex + PC1 + PC2 + PC3 + PC4 + n, family=binomial(link="logit"))
  all_res[[i]] = summary(model)$coefficients[8,]

  }
  all_res = do.call("bind_rows",all_res)
  colnames(all_res) = c("beta","se","z","p")
  all_res$hla_allele = alleles
  all_res$gene = allele
  all_res$anc = ancestry
  overall_res_table[[length(overall_res_table)+1]] <<- all_res

  all_res = all_res %>% arrange(desc(beta))

  all_res$hla_allele = factor(all_res$hla_allele,ordered=T,levels = all_res$hla_allele)
  outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_allele_twodigit",allele,"_",ancestry,".png")
  p = ggplot(all_res,aes(beta,hla_allele))+
    geom_errorbarh(mapping = aes(xmin = beta - 1.96*se,xmax = beta +1.96*se,y=hla_allele),height=0.1)+
    geom_point(size=3,shape=21,fill="orange")+
    geom_vline(xintercept=0,linetype="dashed",color="pink")+
    labs(y=paste0(allele," allele"),x="Beta (log Odds Ratio)")+
    theme_bw()+
    ggtitle(paste0(toupper(ancestry),"\n",allele))


  png(outfile,res=900,units="in",width=4,height=4)
  print(p)
  dev.off()
}  

# make parameter table
alleles = c("A","B","C","DQB1","DQA1","DPA1","DPB1","DRB1","DRB3","DRB4","DRB5")
ancestries = c("sas","afr")
param_tbl = expand.grid(alleles,ancestries)

# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest_twodigit(param_tbl$Var2[i],param_tbl$Var1[i])
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",hla_allele))

bonf = 0.05/11
overall_res_table$direction = ifelse(overall_res_table$beta>0,"up","down")
p =  ggplot(overall_res_table,aes(full_allele,-log10(p),fill=toupper(anc),shape=direction))+
    geom_point(data = overall_res_table %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(bonf),linetype="dashed",color="pink")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="Ancestry",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table %>% filter(p<0.05),mapping = aes(label = full_allele))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_twodigit_all.png",res=900,units="in",width=12,height=6)
p
dev.off()



````

#### HLA-TAPAS
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_imputation.sh
````


#### HLA associations
````unix

cd /data/home/hmy117/HLA-TAPAS
for ancestry in sas afr;
do
module load R/3.6.1
module load python/3.8

# make plink file with phased vcf
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus

awk 'NR>1{print 0,$2,$3-1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_omnibus.tsv

awk 'NR==1{print "FID","IID","MS_status"};NR>1{print 0,$2,$3}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv

awk 'NR==1{print "FID",$2,$3,$4,$5,$6,$7,$8};NR>1{print 0,$2,$3,$4,$5,$6,$7,$8}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_with_pcs.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv

awk 'NR==1{print "FID",$2,$3,$4,$5,$6,$7,$8};NR>1{print 0,$2,$3,$4,$5,$6,$7,$8}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_with_pcs.tsv > \
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_logistic.tsv


python3 -m HLAassoc OMNIBUS_LOGISTIC \
--vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--bim /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus.bim \
--fam  /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus.fam \
--covars /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_omnibus_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_omnibus.tsv \
--aa-only \
--maf-threshold 0.05

# pick up here - only done in AFR
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05

done

# condition on rs3130041
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond1_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition rs3130041

echo rs3130041 > condition_list
echo rs2516652 >> condition_list
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond2_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition-list condition_list


# condition on rs3130242
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond1_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition rs3130242

# condition on rs3130242 & rs2534803
echo rs3130242 > condition_list
echo rs2534803 >> condition_list
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond2_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition-list condition_list

# condition on rs3130242 & rs2534803 & rs2534808
echo rs2534808 >> condition_list
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond3_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition-list condition_list


echo rs111455094 >> condition_list
~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_cond4_$ancestry \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_pheno_for_logistic.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_$ancestry\_covars_for_omnibus.tsv \
--glm  hide-covar \
--maf 0.05 \
--condition-list condition_list




````

##### Analyze omnibus
````unix

library(tidyverse)
ancestry="sas"

logistic_infile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_logistic_",ancestry,".MS_status.glm.logistic.hybrid")
dat = read_table(logistic_infile)

# define variation type
dat = dat %>%
mutate(type = case_when(
  grepl("^HLA",ID) ~ "HLA allele",
  grepl("^AA",ID) ~ "Amino acid change",
  .default = as.character("SNP")
  ))

# just show HLA

ggplot(dat %>% filter(type=="HLA allele"),aes(POS,-log10(P),label=ID))+
geom_point()+
theme_bw()+
labs(x="Position on chr6 (hg19)",y=bquote(-log[10]~P))+
geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
geom_hline(yintercept=-log10(1e-5),linetype="dashed",color="pink")+
ggrepel::geom_text_repel(data = dat %>% filter(type=="HLA allele" & P < 0.05) ,min.segment.length=0)+
ggtitle(toupper(ancestry))


ggplot(dat,aes(POS,-log10(P),label=ID))+
geom_point()+
theme_bw()+
facet_wrap(~type,nrow=3)+
labs(x="Position on chr6 (hg19)",y=bquote(-log[10]~P))+
geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
geom_hline(yintercept=-log10(1e-5),linetype="dashed",color="pink")+
ggrepel::geom_text_repel(data = dat %>% slice_min(P),min.segment.length=0)+
ggtitle(toupper(ancestry))


cond_list = read_table("condition_list",col_names = F)





in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_omnibus_",ancestry,".txt")
dat = read_table(in_file) %>%
mutate(fdr = p.adjust(PVALUE,method="fdr"))

ggplot(dat %>% filter(grepl("^AA_",AA_ID)),aes(AA_POS,-log10(PVALUE),label=AA_POS))+
geom_point()+
facet_wrap(~GENE)+
theme_bw()+
labs(x="Amino acid position",y=bquote(-log[10]~P))+
geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
geom_hline(yintercept=-log10(0.05/8),linetype="dashed",color="pink")+
ggrepel::geom_text_repel(data = dat %>% filter(grepl("^AA_",AA_ID)) %>% filter(PVALUE < 0.05))+
ggtitle(toupper(ancestry))


ggplot(dat,aes(POS,-log10(PVALUE)))+
geom_point()+
theme_bw()+
labs(x="Position on chr6 (hg19")+
geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")

````



##### Define haplotypes
````unix
for ancestry in sas afr;
do
# get HLA alleles
cd /data/scratch/hmy117
grep HLA /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus.bim > hla_alleles
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus \
--extract hla_alleles \
--make-bed \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla

# recode as additive
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla \
--recode A include-alt \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla_ad
done




~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla \
--r2 in-phase dprime \
--ld-window-kb 999999 \
--ld-window-r2 0

~/plink2 --vcf /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--export haps

````


#### Compare HIBAG vs SNP2HLA
````R

library(tidyverse)


all_res = list()
check_hla_allele_agreement = function(ancestry,allele){
    message(ancestry)
    message(allele)

    in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_allele_",allele,"_",ancestry,".tsv")
    hibag_calls = read_tsv(in_file,col_types="cccdd")
    colnames(hibag_calls) = c("IID","allele1","allele2","prob","match")
    hibag_calls$allele1 = substr(hibag_calls$allele1,0,2)
    hibag_calls$allele2 = substr(hibag_calls$allele2,0,2)

    in_file2 = paste0("/data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_",ancestry,"_for_hla_imp_omnibus_just_hla_ad.raw")
    snp2hla_calls = read_table(in_file2)

  # process hibag
    hibag_calls = hibag_calls %>%
      dplyr::select(c(1:3)) %>%
      pivot_longer(cols = c(2,3)) %>%
      dplyr::select(-name) %>%
      group_by(IID) %>%
      dplyr::count(value) %>%
      ungroup() %>%
      mutate(name = paste0("HLA_",allele,"*",value))

      snp2hla_calls = snp2hla_calls %>%
        pivot_longer(cols = contains("HLA")) %>%
        dplyr::select(IID,name,value) %>%
        filter(grepl(allele,name)) %>%
        mutate(name = str_remove_all(name,"_A\\(\\/T\\)$"))

  # get two field res
  combo = snp2hla_calls %>%
        filter(name %in% hibag_calls$name) %>%
        mutate(value = 2-value) %>%
        dplyr::rename("snp2hla_dose" = value) %>%
        mutate(IID = as.character(IID)) %>%
        left_join(hibag_calls,by=c("IID","name")) %>%
        mutate(hibag_dose = ifelse(is.na(n),0,n)) %>%
        dplyr::select(1,2,3,6)

  combo = combo %>%
      mutate(agreement = snp2hla_dose == hibag_dose) %>%
      group_by(name) %>%
      dplyr::count(agreement) %>%
      mutate(prop = n/sum(n)) %>%
      filter(agreement == T) %>%
      mutate(gene = allele,anc=ancestry)
  all_res[[length(all_res)+1]] <<- combo
}

# make parameter table
alleles = c("A","B","C","DQB1","DQA1","DPA1","DPB1","DRB1","DRB3","DRB4","DRB5")
ancestries = c("sas","afr")
param_tbl = expand.grid(alleles,ancestries)

# loop through and check agreement
for(i in c(1:nrow(param_tbl))){
  check_hla_allele_agreement(param_tbl$Var2[i],param_tbl$Var1[i])
}


# plot agreement
all_res = do.call("bind_rows",all_res)
all_res$name =  str_remove_all(all_res$name,"HLA_")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_snp2hla_vs_hibag.png",res=900,width=14,height=4,units="in")
ggplot(all_res,aes(name,prop,fill=gene))+
  geom_col(color="black")+
  facet_wrap(~toupper(anc),nrow=2)+
  scale_fill_brewer(palette="Paired")+
  labs(fill="Gene",x="HLA allele",y="Proportion agreement\n(SNP2HLA vs HIBAG)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))
dev.off()

````


#### VEP
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

module load R/4.2.2


# liftover
# for ancestry in sas eur afr
for ancestry in sas afr
do
	Rscript ./scripts/liftover_hg19_to_hg38.R ./outputs/case_control_gwas_$ancestry\_MS_status.regenie
	awk 'NR>1{if($15>3) print $1,$18,$18,$3"/"$4,"+",$2}' ./outputs/case_control_gwas_$ancestry\_MS_status.regenie_hg38 > ./outputs/susceptibility_vep_input_$ancestry

# annotate
module load ensembl-vep
~/ensembl-vep/vep -i ./outputs/susceptibility_vep_input_$ancestry -o ./outputs/snp_annotations_nearest_susceptibility_$ancestry \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence,Existing_variation"
done

````


#### Combine with IMSGC risk GWAS
##### Liftover to hg38
````unix

# lift over to hg38
cd /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS
awk 'NR>1{print "chr"$1,$2-1,$2,$3}' discovery_metav3.0.meta > hg19_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg19_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
hg38_bedfile \
unmapped
````

##### Liftover to hg38
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS")
imsgc_risk = read_table("discovery_metav3.0.meta")
hg38_positions = read_table("hg38_bedfile",col_names=F)

# prepare for merge
hg38_positions = hg38_positions %>%
			dplyr::select(X3,X4) %>%
			dplyr::rename("SNP" = X4,"BP_hg38" = X3)

# find duplicate positions
dups = hg38_positions %>%
	dplyr::count(SNP) %>%
	filter(n>1)

# filter to biallelics
hg38_positions = hg38_positions %>%
	filter(!SNP %in% dups$SNP)

# filter
imsgc_risk = imsgc_risk %>%
	filter(SNP %in% hg38_positions$SNP)

# merge
imsgc_risk = imsgc_risk %>%
	left_join(hg38_positions,by="SNP")

# add beta
imsgc_risk = imsgc_risk %>%
  mutate(BETA = log(OR))

# add se
imsgc_risk = imsgc_risk %>%
  mutate(
    Z = qnorm(1 - (P / 2)  ) ,
    SE = BETA / Z
    )

# replace hg19 positions
imsgc_risk = imsgc_risk %>%
  dplyr::select(-BP,-SNP) %>%
  dplyr::rename("BP" = BP_hg38) %>%
  mutate(SNP = paste0(CHR,":",BP))

# write to file
write_tsv(imsgc_risk,"imsgc_ms_risk_discovery_hg38.tsv")
````

#### Explore GWAS results
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
# read in GWAS
read_regenie_gwas = function(x){
  read_table(x, col_types = cols_only(
		  CHROM = col_double(),
		  GENPOS = col_double(),
		  ID = col_character(),
		  ALLELE0 = col_character(),
		  ALLELE1 = col_character(),
		  A1FREQ = col_double(),
      A1FREQ_CASES = col_double(),
      A1FREQ_CONTROLS = col_double(),
		  N = col_double(),
      N_CONTROLS = col_double(),
      N_CASES = col_double(),
		  TEST = col_character(),
		  BETA = col_double(),
		  SE = col_double(),
		  CHISQ = col_double(),
		  LOG10P = col_double())) %>%
  filter(TEST == "ADD" & !is.na(LOG10P)) %>%
  dplyr::rename("CHR" = CHROM, "BP" = GENPOS, "SNP" = ID) %>%
  dplyr::select(-TEST) %>%
	mutate(P = 10^-LOG10P)
}
read_plink_gwas = function(x){
  read_table(x, col_types = cols_only(
  `#CHROM` = col_double(),
  POS = col_double(),
  ID = col_character(),
  REF = col_character(),
  ALT = col_character(),
  A1 = col_character(),
  `FIRTH?` = col_character(),
  TEST = col_character(),
  OBS_CT = col_double(),
  OR = col_double(),
  `LOG(OR)_SE` = col_double(),
  Z_STAT = col_double(),
  P = col_double(),
  ERRCODE = col_character()
)) %>%
  filter(TEST == "ADD" & !is.na(P)) %>%
  dplyr::rename("CHR" = `#CHROM`, "BP" = POS, "SNP" = ID) %>%
  dplyr::select(-TEST)
}

# define ancestries
# ancestries=c("afr","sas","eur")
ancestries=c("afr","sas")
# ancestries=c("afr")

# read in plink gwas
all_res = purrr::map(ancestries,function(x){
	gwas_res = read_plink_gwas(
		paste0("./outputs/plink_risk_gwas_",x,".MS_status.glm.logistic.hybrid")) %>%
	mutate(anc = x) %>%
	filter(!is.na(ERRCODE))
	gwas_res
})

# combine
all_res_plink = do.call("bind_rows",all_res)

# get gc
all_res_plink$chisq = qchisq(all_res_plink$P,df=1,lower.tail=F)
lambdas = all_res_plink %>%
	group_by(anc) %>%
	summarise(median_chsq = median(chisq)) %>%
	mutate(lambda = median_chsq / qchisq(0.5,df=1,lower.tail=F))
all_res_plink = all_res_plink %>%
	left_join(lambdas,by="anc") %>%
	mutate(chsq_adj = chisq / lambda) %>%
	mutate(P_adj = pchisq(chsq_adj,df=1,lower.tail=F))

qqman::manhattan(all_res_plink %>% filter(anc=="afr"),main="AFR")

qqman::manhattan(all_res_plink %>% filter(anc=="sas") %>% mutate(P = P+1e-100),main="SAS")

# read in gwas
all_res = purrr::map(ancestries,function(x){
	gwas_res = read_regenie_gwas(
		paste0("./outputs/case_control_gwas_",x,"_MS_status.regenie_hg38")) %>%
	mutate(anc = x)
	gwas_res
})


# combine
all_res = do.call("bind_rows",all_res)


qqman::manhattan(all_res %>% filter(anc=="sas"),main="SAS")
qqman::manhattan(all_res %>% filter(anc=="afr"),main="AFR")

# read in gwas (hg19)
all_res_hg19 = purrr::map(ancestries,function(x){
	gwas_res = read_regenie_gwas(
		paste0("./outputs/case_control_gwas_",x,"_MS_status.regenie")) %>%
	mutate(anc = x)
	gwas_res
})


# combine
all_res_hg19 = do.call("bind_rows",all_res_hg19)

# comparison of plink vs regenie

all_res_combo = all_res_plink %>%
  left_join(all_res_hg19,by=c("SNP","anc"))

# define MHC SNPs
all_res_combo = all_res_combo %>%
  mutate(in_mhc = ifelse(CHR.x==6 & BP.x > 25e6 & BP.y < 35e6,"MHC","Non-MHC"))
png("./outputs/plink_vs_regenie_risk.png",res=900,units="in",width=8,height=4)
ggplot(all_res_combo,aes(-log10(P.x),-log10(P.y),fill=in_mhc))+
  geom_abline(slope=1,intercept=0,color="pink")+
  geom_hline(yintercept=5,color="pink")+
  geom_vline(xintercept=5,color="pink")+
  geom_point(shape=21,color="black")+
  facet_wrap(~toupper(anc))+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(fill="In MHC?",x=bquote(-log[10]~P~Fixed~effects(PLINK)),y=bquote(-log[10]~P~Mixed~effects(REGENIE)) )
dev.off()


# manhattans
manhattan_plots = list()
make_manhattans = function(ancestry){
	gwas_res = all_res %>% filter(anc == ancestry)

	# read annotations

	anno = read_table(paste0("./outputs/snp_annotations_nearest_susceptibility_",ancestry),skip=31,col_types = cols(.default="c"))
  anno = anno %>%
		dplyr::select(1,NEAREST)
	colnames(anno) = c("SNP","Gene")
	anno = anno %>% distinct()

  # get coords
	coords = gwas_res %>%
		group_by(CHR) %>%
		summarise(min_bp = min(BP),max_bp=max(BP),median_bp = median(BP))

	new_coords = coords %>%
		mutate(CHR = CHR+1) %>%
		mutate(cumbp = cumsum(max_bp))

	midpoints = coords$median_bp+c(0,new_coords$cumbp[1:21])
	gwas_res = gwas_res %>%
			left_join(new_coords,by="CHR") %>%
			mutate(cum_bp = ifelse(is.na(cumbp),BP,BP+cumbp))

	# colors
	gwas_res = gwas_res %>%
		mutate(significance = case_when(
			P < 1e-5 ~ "sig",
			CHR %% 2 == 0 ~ "even",
			CHR %% 2 != 0 ~ "odd"
			))

	col_pal = c("lavenderblush1","lavenderblush2","orange")
	names(col_pal) = c("even","odd","sig")

	plot_dat = gwas_res %>% filter(P < 0.5)
	p=ggplot(plot_dat,aes(cum_bp,-log10(P),col=significance))+
		geom_hline(yintercept = -log10(1e-5),color="blue",linetype="dashed",alpha=0.5)+
		geom_point()+
		scale_color_manual(values = col_pal)+
		scale_x_continuous(breaks = midpoints,labels = c(1:22))+
		ggrepel::geom_text_repel(plot_dat %>%
			filter(P < 1e-4 & SNP %in% anno$SNP) %>%
			left_join(anno,by="SNP") %>%
			distinct(Gene,.keep_all=T),
		mapping = aes(label = Gene),min.segment.length = 0,color="black",nudge_y=1,direction = "x")+
		labs(x="Genomic position (hg38)",y=bquote(-log[10]~P))+
		ggtitle(paste0("GWAS of MS risk in ",toupper(ancestry)))+
		theme_bw()+
		theme(legend.position="none",panel.grid=element_blank())


	png(paste0("./outputs/manhattan_plot_susceptibility",ancestry,".png"),res=900,units="in",width=10,height=4)
	print(p)
	dev.off()
	manhattan_plots[[length(manhattan_plots)+1]] <<- p
  # check allele counts
  gwas_res = gwas_res %>%
  mutate(AC_case = A1FREQ_CASES * (2 * N_CASES ) ) %>%
  mutate(AC_cont = A1FREQ_CONTROLS * (2 * N_CONTROLS ) )

  # filter to MAC 20 in cases and controls
  gwas_res = gwas_res %>% filter(AC_case >= 20 & AC_cont >= 20)

  # write annotations
  gwas_res %>%
    left_join(anno,by="SNP") %>%
    filter(P<5e-8) %>%
    dplyr::select(SNP,ALLELE0,ALLELE1,A1FREQ,BETA,P,Gene) %>%
    group_by(Gene) %>%
    slice_min(P) %>%
    arrange(P)

}

sapply(ancestries,make_manhattans)
# print all together
library(gridExtra)
png("./outputs/all_manhattan_plots_susceptibility.png",res=900,units="in",width=10,height=8)
cowplot::plot_grid(plotlist = manhattan_plots,align="v",ncol=1)
dev.off()


# read all annotations

all_anno = list()
for(ancestry in c("sas","afr")){
	anno = read_table(paste0("./outputs/snp_annotations_nearest_susceptibility_",ancestry),skip=31,col_types = cols(.default="c"))
  anno = anno %>%
		dplyr::select(1,NEAREST)
	colnames(anno) = c("SNP","Gene")
	anno = anno %>% distinct() %>%
  mutate(anc = ancestry)
  all_anno[[length(all_anno)+1]] = anno
}
all_anno = do.call("bind_rows",all_anno)
sig_hits = all_res %>%
  filter(P < 1e-5) %>%
  left_join(all_anno,by=c("SNP","anc"))
# save sig snps
write_csv(sig_hits,"./outputs/sig_snps_risk.csv")



# inflation
qqplots = list()
make_qq_plot = function(ancestry){
	gwas_res = all_res %>% filter(anc == ancestry)
	lambda = median(qchisq(gwas_res$P,df=1,lower.tail=F)) /  qchisq(0.5,lower.tail=F,df=1)
	message(ancestry)
	message(lambda)
	pval_dat = data.frame(observed = -log10(gwas_res$P)) %>%
		arrange(observed)
	pval_dat = pval_dat %>%
		mutate(expected = -log10(seq(1,1/nrow(pval_dat),by=-1/nrow(pval_dat))))

	#	downsample but ensure top snps are included
	pval_dat$index = seq(1,nrow(pval_dat),by=1)
	tophits = tail(pval_dat,round(nrow(pval_dat)*0.01,0))
	sampled_pval_dat = sample_frac(pval_dat,0.05) %>% bind_rows(tophits) %>% distinct()

	# plot
	p = ggplot(sampled_pval_dat,aes(expected,observed))+
		geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+
		geom_point()+
		labs(x="Expected -log10(P)",y="Observed -log10(P)")+
		ggtitle(paste0("Ancestry: ",toupper(ancestry),"\nLambda =  ",round(lambda,2)))+
		theme_bw()

	qqplots[[length(qqplots)+1]] <<- p
	png(paste0("./outputs/qq_plot_susceptibility_",ancestry,".png"),res=900,units="in",width=4,height=4)
	print(p)
	dev.off()

}
sapply(ancestries,make_qq_plot)

# print all together
library(gridExtra)
png("./outputs/all_qq_plots.png",res=900,units="in",width=4,height=8)
cowplot::plot_grid(plotlist = qqplots,align="v",ncol=1)
dev.off()


# print all together
png("./outputs/all_manhattans_plots.png",res=900,units="in",width=10,height=10)
cowplot::plot_grid(plotlist = manhattans,align="v",ncol=1)
dev.off()

# compare with IMSGC GWAS
# read IMSGC risk GWAS (hg38 )
imsgc_risk_hg38 = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv",col_types="dccdddddddc")

# filter & join
# join on chr:pos
all_res_with_imsgc = imsgc_risk_hg38 %>%
		dplyr::rename("ALLELE1" = A1,"ALLELE0" = A2) %>%
		mutate(chrpos = paste0(CHR,":",BP)) %>%
		mutate(anc = "IMSGC") %>%
		bind_rows(all_res %>% mutate(chrpos = paste0(CHR,":",BP)))

# get sig snps in IMSGC
all_res_with_imsgc %>% filter(anc=="IMSGC" & P < 5e-7) %>% dplyr::select(chrpos,SNP,P)

# forest
make_forest = function(snp){
  plot_dat = all_res_with_imsgc %>% filter(chrpos == snp)

  # Add study label
  plot_dat = plot_dat %>%
  	mutate(anc = ifelse(anc == "IMSGC","IMSGC-EUR",paste0("ADAMS-",toupper(anc))))

  # flip beta if necessary
  plot_dat = plot_dat %>%
  	mutate(BETA = ifelse(ALLELE1 == plot_dat$ALLELE1[1],BETA,BETA*-1))

  # refactor
  plot_dat$anc = factor(plot_dat$anc,levels = c("IMSGC-EUR","ADAMS-EUR","ADAMS-SAS","ADAMS-AFR"),ordered=T)
  pictured_allele = paste0(snp,"-",plot_dat$ALLELE1[1])
  ggplot(plot_dat,aes(BETA,anc,fill=anc))+
  	geom_errorbarh(mapping = aes(xmin = BETA - 1.96*SE,xmax = BETA + 1.96*SE,y=anc),height=0.3)+
  	geom_point(shape=21,color="black",size=3)+
  	geom_vline(xintercept=0,linetype="dashed")+
  	theme_bw()+
  	labs(x=paste0("Effect of ",pictured_allele," on gARMSS"),y="Study & ancestry")+
  	theme(legend.position="none")
}

# PICK UP HERE 05-1
snp = "2:71449869"
png("./outputs/dysf_snp.png",res=900,units="in")
make_forest("2:71449869")
dev.off()

png("./outputs/edem_snp.png",res=900,units="in",width=6,height=6)
make_forest("3:5290006")
dev.off()

# region plot

make_locus_plot = function(chr,pos,window=10e6){
  plot_dat = all_res_with_imsgc %>% filter(CHR == chr & BP > pos - window/2 & BP < pos + window/2)

  # Add study label
  plot_dat = plot_dat %>%
  	mutate(anc = ifelse(anc == "IMSGC","IMSGC-EUR",paste0("ADAMS-",toupper(anc))))

  # refactor
  plot_dat$anc = factor(plot_dat$anc,levels = c("IMSGC-EUR","ADAMS-EUR","ADAMS-SAS","ADAMS-AFR"),ordered=T)

  p1 = ggplot(plot_dat %>% filter(anc=="IMSGC-EUR"),aes(BP,-log10(P),fill=anc))+
  	geom_point(shape=21,color="black",size=3)+
  	geom_hline(yintercept=-log10(1e-5),linetype="dashed")+
  	theme_bw()+
  	theme(legend.position="none")
	p2 = ggplot(plot_dat %>% filter(anc=="ADAMS-SAS"),aes(BP,-log10(P),fill=anc))+
  	geom_point(shape=21,color="black",size=3)+
  	geom_hline(yintercept=-log10(1e-5),linetype="dashed")+
  	theme_bw()+
  	theme(legend.position="none")
	p3 = ggplot(plot_dat %>% filter(anc=="ADAMS-AFR"),aes(BP,-log10(P),fill=anc))+
  	geom_point(shape=21,color="black",size=3)+
  	geom_hline(yintercept=-log10(1e-5),linetype="dashed")+
  	theme_bw()+
  	theme(legend.position="none")
	cowplot::plot_grid(p1,p2,p3,align="v",ncol=1)

}
png("./outputs/mhc_locus.png",res=900,units="in",width=6,height=6)
make_locus_plot(chr=6,pos=30e6)
dev.off()

# beta beta plot
ancestry = "sas"

plots = list()
make_beta_beta_plot = function(ancestry){
  all_res$chrpos = paste0(all_res$CHR,":",all_res$BP)

  filtered_imsgc = imsgc_risk_hg38 %>%
  	filter(SNP %in% all_res$chrpos)

  filtered_imsgc_this_anc = filtered_imsgc %>%
  	dplyr::rename("chrpos" = SNP) %>%
  	left_join(all_res %>%
  		filter(anc == ancestry),by="chrpos") %>%
  		filter(P.x < 5e-8) %>%
  		mutate(BETA.y = ifelse(A1 == ALLELE1,BETA.y,BETA.y*-1))

  spearman_cor = cor.test(filtered_imsgc_this_anc$BETA.x,filtered_imsgc_this_anc$BETA.y,method="spearman")
  print(spearman_cor)
  cor = cor.test(filtered_imsgc_this_anc$BETA.x,filtered_imsgc_this_anc$BETA.y)
  print(cor)
  lab = paste0(toupper(ancestry),"\nr=",round(cor$estimate,2))
  p = ggplot(filtered_imsgc_this_anc,aes(BETA.x,BETA.y))+
    geom_point()+
    labs(x="Beta (IMSGC)",y="Beta (ADAMS)")+
    geom_abline(slope=1,intercept=0,linetype="dashed",color="pink")+
    geom_vline(xintercept=0,linetype="dashed",color="pink")+
    geom_hline(yintercept=0,linetype="dashed",color="pink")+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    ggtitle(lab)
  plots[[length(plots)+1]] <<- p
}

make_beta_beta_plot("sas")
make_beta_beta_plot("afr")

png("./outputs/beta_beta_plots.png",res=900,units="in",width=8,height=4)
cowplot::plot_grid(plotlist = plots,align="h",ncol=2)
dev.off()

````
