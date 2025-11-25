
# Genetic analysis of Multiple Sclerosis in diverse ancestries
### A genetic _A_ssociation Study of individuals of _D_iverse _A_ncestral backgrounds with _M_ultiple _S_clerosis (ADAMS)<a href =https://app.mantal.co.uk/adams><p><center> <img src="https://app.mantal.co.uk/files/db2440daa093e6b2d688b8a02b82cff2.png" width=20% height=40%></a></p></center>


This repo contains code used to analyze data from the ADAMS project - a genotype-phenotype cohort of people with Multiple Sclerosis focussed on those from diverse ancestral backgrounds. The baseline cohort description is [here](https://bmjopen.bmj.com/content/13/5/e071656.full).

# Code
Prior to imputation genotypes were called using Illumina Genome Studio v2.0. Raw PLINK binary files (.map and .ped) were exported from Genome Studio with calls coerced to the forward strand.

## Imputation
### Initial genotype quality control
First the .ped and .map files are converted to PLINK1 binary format with some light QC.
- Conversion to PLINK1 binary
- Basic QC pre-imputation: autosomes, missingess for genotypes and people <10%, MAC > 1, HWE deviation

````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23

~/plink --file /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/raw/ADAMS_19_06_24 \
--make-bed \
--chr 1-22 \
--geno 0.1 \
--hwe 1e-10 \
--mac 1 \
--mind 0.1 \
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

cut -f1,2 ./outputs/king.con > ./outputs/ids_to_exclude.tsv

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

awk '{print $4,$1":"$3}' hg19_bedfile > hg19_snps
awk '{print $4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions and IDs
~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid_inpheno_nodups \
--update-map hg19_snp_positions \
--make-bed \
--out adams_temp_hg19_for_bigsnpr

~/plink2 --bfile adams_temp_hg19_for_bigsnpr \
--set-all-var-ids chr@:#:\$r\:\$a \
--make-bed \
--rm-dup exclude-all \
--out adams_hg19_for_bigsnpr_cpra

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
for i in {1..22};
do
	~/plink --bfile /data/scratch/hmy117/hgdp_1kg_genomes/adams_hg19_for_bigsnpr_cpra \
	--chr $i \
	--recode vcf \
	--out ./imputation_raw_files/chr$i &
done

### Sort with bcf tools
module load bcftools
for i in {1..22};
do
  bcftools sort ./imputation_raw_files/chr$i\.vcf \
  -Oz -o ./imputation_raw_files/sorted_chr$i\.vcf.gz &
done
````

### First pass imputation
- Attempt imputation at MIS
- This will fail due to strand flips
- Download excluded snps and use strand flips to flip

````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/imputation_raw_files/
# download snps_excluded.txt

awk 'FS="\t"{if($2=="Strand flip") print $1}' snps-excluded.txt > snps_to_flip
awk 'FS="\t"{if($2!="Strand flip") print $1}' snps-excluded.txt | uniq > snps_to_exclude

# set var IDs to chr:pos:ref:alt
for i in {1..22};
do
	~/plink2 --vcf chr$i\.vcf \
	--set-all-var-ids @:#:\$r\:\$a \
	--make-bed \
	--rm-dup exclude-all \
	--out updated_ids_chr$i &
done

# recode VCFs
for i in {1..22};
do
	~/plink --bfile updated_ids_chr$i \
	--exclude snps_to_exclude \
	--out flipped_chr$i \
	--recode vcf &
done

# flip strands
for i in {1..22};
do
	~/plink --vcf flipped_chr$i\.vcf --double-id \
	--flip snps_to_flip \
	--out dbl_flipped_chr$i \
	--recode vcf &
done

	# sort with bcf tools
module load bcftools
for i in {1..22};
do
	bcftools sort dbl_flipped_chr$i\.vcf \
	 -Oz -o sorted_chr$i\.vcf.gz &
done
````

###########################
### Download imputed data
###########################

### 2nd pass imputation
- Performed imputation again
- Download imputed data
- Rsq filter 0.3
- MIS

````unix
cd /data/scratch/hmy117/adams_imputed

wget https://imputationserver.sph.umich.edu/share/results/8b684d5743ea0d75b84192845fb11b760b1dba65391eb606893f8d637bb1c815/chr_1.zip &
wget https://imputationserver.sph.umich.edu/share/results/a0591a8396007a67d4abff86b1133638529e40f5a4d6dd64c590edb1f3c980e8/chr_10.zip &
wget https://imputationserver.sph.umich.edu/share/results/3adfd245fc04c8b0c1d1cfde04b5040f97172486ac5642f2d47f456562cb94d6/chr_11.zip &
wget https://imputationserver.sph.umich.edu/share/results/82abc9f8b5cd763e72b6f61e19fdc98256c5003d7e6b812e241f8e345d7fed40/chr_12.zip &
wget https://imputationserver.sph.umich.edu/share/results/b8b8a11b3b88bbc0d057d27f3700cbd80566a0cf6838e545c6e3604d8cd183ca/chr_13.zip &
wget https://imputationserver.sph.umich.edu/share/results/bbaa5525b6403a4ce89bf807b7332c4499a5327cb5ecdee4f332400dcf2e8474/chr_14.zip &
wget https://imputationserver.sph.umich.edu/share/results/cfa7a86d5448ac6dea63db4a463fe31ca88c502f0b4574c4110eb1f8dd3a0407/chr_15.zip &
wget https://imputationserver.sph.umich.edu/share/results/4d2390886c5974e7defb979bd7b40970e995ae2b45e78f02e247cd3a98a1f544/chr_16.zip &
wget https://imputationserver.sph.umich.edu/share/results/c3217964bb276a0cb06237e2f767757deebf97e42c09842dee9af35b58857279/chr_17.zip &
wget https://imputationserver.sph.umich.edu/share/results/275717953a09f44661fd1b922cb713eafe40c98cc5cc0ba890e19de027744fd0/chr_18.zip &
wget https://imputationserver.sph.umich.edu/share/results/4547efd5c45ce766486f7088136cf88d33d46be6480350c3f2ec5d0a43f52a1e/chr_19.zip &
wget https://imputationserver.sph.umich.edu/share/results/baf4bf599a9a2e36c0e2abc644d1f283faa9aa01799509e74b0c052b1777e0ca/chr_2.zip &
wget https://imputationserver.sph.umich.edu/share/results/bdc76ad02e442faa44b8d158675cbf7981247039fa438080dec7478a6442be1c/chr_20.zip &
wget https://imputationserver.sph.umich.edu/share/results/aa271ad77a984858586b54c4d9fa5bc6ab6c5959d09993127e66a0e641a37d87/chr_21.zip &
wget https://imputationserver.sph.umich.edu/share/results/cd373d0e7426f5ed5afec9213ed957066965fa6f9c909dbd8a7c5c5bb798dbe4/chr_22.zip &
wget https://imputationserver.sph.umich.edu/share/results/ad160bc24ca6c3ca7806ae385b221c8de76f5b13df9841d951ad894222165399/chr_3.zip &
wget https://imputationserver.sph.umich.edu/share/results/836831970bfab56eec74f0798bc3466553935bec59ed9c3d8e0ef485ddcd3fea/chr_4.zip &
wget https://imputationserver.sph.umich.edu/share/results/cacd2b4a688f0ce48666fe91752cbe6dc825a9f72d4ba5125a02da706368b032/chr_5.zip &
wget https://imputationserver.sph.umich.edu/share/results/664d67ba617bc6218a37112593fff2e5f6031b0d5f1bcec6617edb2897d709f6/chr_6.zip &
wget https://imputationserver.sph.umich.edu/share/results/c899694b6950fe036942b309da5a0b4fb66ec2857b2d1985e353aa8309efcf29/chr_7.zip &
wget https://imputationserver.sph.umich.edu/share/results/19068676962bf362598800b539cbb1fbdc167155c736afc1ae2dc59bf845d0ed/chr_8.zip &
wget https://imputationserver.sph.umich.edu/share/results/d2953be197b4836d17c3cf25422767223386951ca735260eae84faa4ffe5b7ec/chr_9.zip &
wget https://imputationserver.sph.umich.edu/share/results/93e474d7c04c40c6a23e983655eba7212a45bbb3eb6db4fdb041ddef1e8449b2/results.md5

# unzip
for i in {1..22};
do
  unzip -o -P "FqOPnzaNsYa>70" chr_$i\.zip &
done

````

## Ancestry inference


### Download reference data
````unix
# running 21-06
# run with qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/download_hgdp_1kg.sh
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

### Download UKB SNP list

cd /data/scratch/hmy117/adams_imputed/
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_mfi.tgz
tar -xvf ukb_imp_mfi.tgz
#### Inpsect in R & get SNPs with high INFO & MAF
````R
Rscript "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/info_snp_filtering.R"
````


### SNP QC
- SNP QC on individual VCFs
- Filter to INFO > 0.3 & MAF > 0.01 <10% and HWE P>1e-10
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


## Susceptibility GWAS


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
adams_snps  = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc.bim",col_names = F)
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


# write snp list to keep
adams_snps_to_keep = combo_snps %>% dplyr::select(ID.x) %>% distinct()
ukb_snps_to_keep = combo_snps %>% dplyr::select(ID.y) %>% distinct()
ukb_snps_names_to_update = combo_snps %>% dplyr::select(ID.y,ID.x) %>% distinct(ID.y,.keep_all=T)


write_tsv(adams_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_snps_to_keep.tsv",col_names = F)
write_tsv(ukb_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_keep.tsv",col_names = F)
write_tsv(ukb_snps_names_to_update,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_update_names.tsv",col_names = F)

````

#### Filter to compatible SNPs
##### UKB
````unix
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_ukb_genos_step1.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/update_ukb_snp_names.sh"
````
##### ADAMS
````unix
~/plink2 \
--bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc \
--extract /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_snps_to_keep.tsv \
--out /data/scratch/hmy117/adams_filtered_hg19_cpra \
--make-bed

````


#### Merge Cam-ADAMS with UKB files across chromosomes
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_ukb_files.sh
````


#### QC per chromosome
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merged_ukb_adams_qc_per_chrom.sh
````

#### Clean up
````unix
# clean up
for i in {1..22}; do rm /data/scratch/hmy117/ukb_adams_merged_genotypes_chr$i\.* ; done
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
	awk '{print $2}' /data/scratch/hmy117/ukb_adams_merged_genotypes_chr$i\_no_palindromes.bim >> adams_ukb_snps_to_keep
done

# filter to UKB-ADAMS vars
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra \
--extract adams_ukb_snps_to_keep \
--out combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry \
--make-bed
````


#### Prune each chromosome (for PCA)
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merged_ukb_adams_prune_for_pca.sh
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
--extract /data/scratch/hmy117/ukb_adams_merged_genotypes_all_chroms_for_pca.bim \
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

# project ADAMS-UKB samples
~/plink2 --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_all_chroms_for_pca \
--read-freq hgdp_kg_pcs_with_ukb.afreq \
--score hgdp_kg_pcs_with_ukb.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize list-variants \
--score-col-nums 6-55 \
--out ukb_adams_pcs

# project original dataset
~/plink2 --bfile combined_hgdp_1kg_unrelated_hg19_cpra_for_ancestry_pruned \
--read-freq hgdp_kg_pcs_with_ukb.afreq \
--extract ukb_adams_pcs.sscore.vars \
--score hgdp_kg_pcs_with_ukb.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize \
--score-col-nums 6-55 \
--out hgdp_kg_pcs_with_ukb_rescored
````

#### Identify ancestry groupings
````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")

# read in data
adams_ukb = read_table("ukb_adams_pcs.sscore")

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

predicted_ancestry_confidence = predict(rf_fit,adams_ukb,type="prob")
predicted_ancestry_confidence$IID = adams_ukb$IID

predicted_ancestry_confidence = tibble(predicted_ancestry_confidence) %>%
	mutate(max_prob = pmax(AFR,AMR,CSA,EAS,EUR,OCE,MID))

ambiguous_calls =  predicted_ancestry_confidence %>% filter(max_prob < 0.80)
adams_ukb = adams_ukb %>%
		mutate(ambiguous = ifelse(IID %in% ambiguous_calls$IID,"yes","no"))

#ggplot(adams_ukb,aes(PC1_AVG,PC2_AVG,fill=ambiguous))+geom_point(shape=21)


cov = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv") %>% dplyr::select(-FID)
adams_fam = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc.fam",col_names=F)

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
		IID %in% adams_fam$X2 ~ "ADAMS"
		))

# define case-control status
ukb_cases = ukb_pheno %>% filter(!is.na(source_of_report_of_g35_multiple_sclerosis_f131043_0_0))

adams_ukb = adams_ukb %>%
	mutate(MS_status = case_when(
		cohort == "UKB" & IID %in% ukb_cases$EID ~ "MS",
		cohort == "ADAMS" ~ "MS",
		.default = as.character("Control")
		))

# ancestry plots
p = ggplot()+
	geom_point(data = adams_ukb %>% filter(ambiguous=="no" & cohort == "UKB"),mapping = aes(PC1_AVG,PC2_AVG,fill=predicted_ancestry),shape=22,alpha=0.1)+
	geom_point(data = adams_ukb %>% filter(ambiguous=="no" & cohort == "ADAMS"),mapping = aes(PC1_AVG,PC2_AVG,fill=predicted_ancestry),shape=21)+
	theme_bw()+
	labs(x="PC1",y="PC2",fill="Ancestry")+
	scale_fill_brewer(palette="Set2")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_exclusions_all_anc.png",res=900,units="in",width=4,height=4)
p
dev.off()


p = ggplot()+
	geom_point(data = adams_ukb %>% filter(!is.na(cohort) & predicted_ancestry %in% c("CSA","AFR","EUR") & ambiguous=="yes"),mapping = aes(PC1_AVG,PC2_AVG),fill="grey",alpha=0.1,shape=21)+
	geom_point(data = adams_ukb %>% filter(!is.na(cohort) & predicted_ancestry %in% c("CSA","AFR","EUR") & ambiguous=="no"),mapping = aes(PC1_AVG,PC2_AVG,fill=cohort),shape=21)+
	facet_wrap(~predicted_ancestry)+
	theme_bw()+
	labs(x="PC1",y="PC2",fill="Cohort")+
	scale_fill_brewer(palette="Set1")
# plot
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_exclusions.png",res=900,units="in",width=8,height=3)
p
dev.off()


ggplot(adams_ukb,aes(PC1_AVG,PC2_AVG,fill=cohort))+
	geom_point(shape=21)

# ambiguous counts
counts_inc_exclusions = adams_ukb %>%
	group_by(cohort) %>%
	dplyr::count(ambiguous,MS_status,predicted_ancestry) %>%
	ungroup() %>%
	pivot_wider(id_cols = c(cohort,MS_status,predicted_ancestry),names_from = ambiguous, values_from=n)
write_csv(counts_inc_exclusions,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_exclusions_tbl.csv")

# remove ambiguous calls
adams_ukb = adams_ukb %>%
	filter(!is.na(cohort) & ambiguous=="no")

# counts
counts = adams_ukb %>%
filter(!is.na(cohort) & ambiguous=="no") %>%
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
	mutate("EID" = as.character(IID)) %>%
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

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_adams_pc_plot_with_ref_pc1_4.png",res=900,units="in",width=8,height=8)
cowplot::plot_grid(p1,p2,align="v",ncol=1)
dev.off()

# validate with self-reported ethnicity
adams_ukb_with_ethnicity = adams_ukb %>%
	filter(!is.na(ethnic_background_f21000_0_0))%>%
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
	left_join(cov,by="IID") %>%
	filter(!is.na(ethnicity_clean))%>%
	filter(predicted_ancestry %in% c("AFR","EUR","CSA"))

table(factor(adams_with_ethnicity$ethnicity_clean),
factor(adams_with_ethnicity$predicted_ancestry))

p = ggplot(adams_with_ethnicity,
aes(ethnicity_clean,fill=predicted_ancestry))+
  geom_bar(position="fill",color="black")+
  theme_bw()+
	coord_flip()+
	scale_fill_brewer(palette="Set1")+
	labs(y="Proportion of individuals",x="Self-reported ethnicity",fill="Inferred genetic ancestry")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_pcs_vs_self_report_risk.png",res=900,units="in",width=8,height=4)
p
dev.off()

# write all ancestry calls
adams_ukb = adams_ukb %>%
	filter(!is.na(cohort))

anc_calls = adams_ukb %>%
  dplyr::select(-IID) %>%
	dplyr::rename("IID" = EID) %>%
  dplyr::select(`#FID`,IID,predicted_ancestry)
write_tsv(anc_calls,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_all_ukb_adams.tsv")

# write pheno - covar
pheno_cov = adams_ukb %>%
  dplyr::select(-IID) %>%
	dplyr::rename("IID" = EID) %>%
	dplyr::select(-contains("PC")) %>%
	dplyr::select(IID,MS_status,sex_f31_0_0,age_at_recruitment_f21022_0_0,cohort) %>%
	dplyr::rename("age" = age_at_recruitment_f21022_0_0, "sex" = sex_f31_0_0) %>%
	mutate(sex = ifelse(sex == "Male",1,2))

# process adams pheno
adams_cov = cov %>%
	mutate(IID = as.character(IID)) %>%
 	dplyr::select(IID,ageatedss,Sex) %>%
	mutate(sex = ifelse(Sex == "Male",1,2)) %>%
	dplyr::rename("age" = ageatedss) %>%
	dplyr::select(IID,age,sex)

non_ukb = pheno_cov %>%
		filter(cohort != "UKB") %>%
		dplyr::select(IID,cohort,MS_status) %>%
		left_join(adams_cov,by="IID")
ukb = pheno_cov %>%
		filter(cohort == "UKB")
combo_dat = bind_rows(non_ukb,ukb)

# remove NAs
combo_dat = combo_dat
write_tsv(combo_dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv")

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
		predicted_ancestry %in% c("MID","OCE","AMR"),"Other",as.character(predicted_ancestry)
		))
tbl = compareGroups(data = combo_dat,
	MS_status ~ cohort + age + gender + Genetic_ancestry,
	method = c(3,2,3,3)
	) %>%
	createTable()

export2csv(tbl,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_gwas_demographics.csv")

````


#### Identify ancestry groupings (part 2)

##### Split into ancestry groups
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/split_ancestry_groups.sh
````
##### Prune, combine across chroms per ancestry, and run PCA per ancestry
###### Filter out PCA outliers within each ancestry
###### RUNNING 23-08

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
		mutate_at(.vars = c("age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),z_score)
	write_tsv(pcs,in_file)
}
````

#### GWAS within ancestry - using whole ancestral cluster - running 2208 1145
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_prep.sh # all within ancestry

# pick up here - running 2208 1630

qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_step1.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_step2.sh"
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
--make-pheno /data/scratch/hmy117/adams_filtered_hg19_cpra.fam \* \
--make-bed \
--allow-no-sex \
--out /data/scratch/hmy117/ms_cases_only_just_ukb_adams

~/plink --bfile /data/scratch/hmy117/ms_cases_only_just_ukb_adams \
--assoc \
--extract /data/scratch/hmy117/snps_to_keep_for_gwaseur.tsv \
--allow-no-sex \
--out /data/scratch/hmy117/ms_cases_ukb_vs_adams

awk '{if($9<5e-8) print}' /data/scratch/hmy117/ms_cases_ukb_vs_adams.assoc > /data/scratch/hmy117/ms_cases_ukb_vs_adams_snps_to_exclude

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
