
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

~/plink --file /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/raw/ADAMS_20_12_23 \
--make-bed \
--chr 1-22 \
--geno 0.1 \
--hwe 1e-10 \
--mac 1 \
--mind 0.1 \
--out ./outputs/ADAMS_geno

````

### Export VCFs

````unix
# flip strand
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
for i in {1..22};
	do
		~/plink --bfile ./outputs/ADAMS_geno \
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

### 2nd pass imputation
- Performed imputation again
- Download imputed data
- Rsq filter 0.3
- TOPMED-r3

````unix
cd /data/scratch/hmy117/adams_imputed
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/15793a0b7deeecd4cc11d7e26952d3f2ecd4519c44b1ca8b567c654e4141574c/chr_1.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/46a7e2523d01a63c7d1dfe50bb37438ac354b1f9815b2df6eda5a3e65f20daad/chr_10.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/fb52c451339eaa2a9b87baaf948cb58dcb2539f149e6625e28f8175e659cb94c/chr_11.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/8840f6502f5e31d2b56d26e2858d794d494745791027089e1e40b31c0c000b67/chr_12.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/8cb1c840ea01cef94aa2fdd6670ecd867964b9526199e4ca6f8270aefc31d47b/chr_13.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/761dc0d69ce4331e2ae0fd38033cbbf43d323ab7f2ea591e725e9a8c34c7ce4a/chr_14.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/f70faa6e000d0aba9a4b508795203985c4af6589831731af25385d2a842d091b/chr_15.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/81bc637da9b201b66de316ae3218dc0bde549fa9978ee24a07618d85bb10719f/chr_16.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/ca7d62b9c0240129297ed74f47dc799495f7ce05747ca2ba52497debca896f34/chr_17.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/18735e2d68adfbff7a7aa55f68618164661d93cf3ee1acc4a3955c632e9be579/chr_18.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/f0ff47444188b5940a5d8e43724d174d54b63821287f28eeeff277b1cf33ca7c/chr_19.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/6b0deeb56a0d439446839a808f2219a6dc6725a775ab291b7a027f1585da8e62/chr_2.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/51abf847119619d6cb7e9d73e05cf06d4411281d3dff3ec0086b8a1d84658105/chr_20.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/fc9b88f649c3ee8695077a68bb3648d9e691842c7cc335adf0bd711f5b065514/chr_21.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/7eb66498b931557796d0349b558d74550f6cae45ac58cb3878ebb140853d779c/chr_22.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/39aa125eafd168e400a06d78f091cfa7f5d7660dc1479e198c8e17bf106efb2e/chr_3.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/b5fb37dfdcb5909adb703d2c0643e6e57652f7e0f3cadc938c6a39de4c2ec9ae/chr_4.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/d51197fe4c386e23c685c9836d990494b64473b8b5f725a0af918184281ed182/chr_5.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/b4a89b06ac38a9dcb7602b57acbbc08064f37745d348e855129dd6ee1f6c613d/chr_6.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/2afc4d699385c493db1a39170cbcd036f5355add2e3defa44e70fd87da7e08e3/chr_7.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/25dd9a461b219d2a6e8447f0e903a20f84b2b61c4348cf211b759b40f65edc36/chr_8.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/a9d1c7759fee50b0ec41f0aea07f8644da6f11aaab238fd5e85c1e3cb0c141e7/chr_9.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/f2486e6eb5fe8c94c5de14c6d6b9595327cba5e5cbbc1e159f6b6425995e87c5/results.md5

# unzip
for i in {1..22};
  do
    unzip -P wSqJCdOGLh87xf chr_$i\.zip &
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
awk '{print "chr"$1,$4-1,$4,$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$1":"$3":"$4}' hg19_bedfile > hg19_snps
awk '{print $1":"$3":"$4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions and IDs
~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno \
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
df = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno.bim",col_names=F) %>%
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
#qsub ./scripts/filter_hgdp_to_adams_vars.sh
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
--out combined_hgdp_1kg
````

### Do same QC for adams
````unix
~/plink --bfile ~/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno \
--maf 0.05 \
--geno 0.01 \
--mind 0.1 \
--hwe 1e-10 \
--extract combined_hgdp_1kg.bim \
--out ADAMS_qc_genos_for_pcs \
--make-bed
````

### Missingness filter
````unix
~/plink --bfile combined_hgdp_1kg \
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

table(adams$predicted_ancestry)
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
  slice_max(pct,n=3) %>%
  print(n=30)
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
# rm info_stats_all_snps
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
write_tsv(snps_qcd,"snps_to_keep_maf1e-2_info_0.3.tsv",col_names=F)

snps_qcd_regenie = snps_qcd_regenie %>%
	mutate(snp_name = paste0(CHR,":",BP,":",REF,":",ALT)) %>%
	dplyr::select(snp_name)
write_tsv(snps_qcd_regenie,"snps_to_keep_step1_regenie.tsv",col_names=F)

# cut into MAF bins
snps$maf_bin = Hmisc::cut2(snps$MAF,cuts = c(0,0.01,0.05,0.1,0.5))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_maf.png",res = 300, units="in",width=4, height=4)
ggplot(snps,aes(maf_bin,INFO,fill=maf_bin))+
	geom_boxplot()+
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

````

### SNP QC
- SNP QC on individual VCFs
- Filter to INFO > 0.7 & MAF > 0.01 & missingness <10% and HWE P>1e-5
- Conversion back to plink
````unix
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_snp_qc_imputed.sh
````

### Update fam files
- Modify IDs in R
````R
library(tidyverse)
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

# exclude indels
snps = snps %>%
	filter(nchar(X5)==1 & nchar(X6)==1)

# exclude duplicate positions (i.e. non-biallelics)
dups = snps %>% dplyr::count(X1,X4)
dups =  dups %>% filter(n>1)

snps = snps %>% left_join(dups,by=c("X1","X4"))
snps = snps %>% filter(is.na(n)) %>% dplyr::select(-n)
write_tsv(snps,"non_duplicated_no_indel_snps.tsv",col_names=F)
````

### QC
````unix
cd /data/scratch/hmy117/adams_imputed/

~/plink --bfile combined_adams_imputed \
--extract non_duplicated_no_indel_snps.tsv \
--make-bed \
--maf 0.01 \
--hwe 1e-5 \
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
anc = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls.tsv")
het = het %>% left_join(anc,by="IID")


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

#### Validate duplicates
````R
library(tidyverse)
kinship  = read_table("king.kin0",col_types=cols(.default="c"))
dups = kinship %>%
	filter(InfType == "Dup/MZ" )
data.frame(ID = c(dups$ID1,dups$ID2))  %>%
distinct(ID)

covars = read_tsv("~/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv",col_types=cols(.default="c"))
matches = dups %>%
	dplyr::select(ID1,ID2,Kinship) %>%
	left_join(covars %>%
			mutate(ID1 = IID) %>% dplyr::select(-IID,-FID),
			by = "ID1") %>%
	left_join(covars %>%
			mutate(ID2 = IID) %>% dplyr::select(-IID,-FID),
			by = "ID2")
matches

# write file for exclusion
exclusion = dups %>% dplyr::select(FID2,ID2) %>% dplyr::rename("FID"=FID2,"IID"=ID2)
write_tsv(exclusion,"dups_to_exclude")
````
#### Exclusion
````unix
~/plink --bfile combined_adams_imputed_qc \
--remove dups_to_exclude \
--make-bed \
--out combined_adams_imputed_qc_unrelated

# copy to home folder
cp combined_adams_imputed_qc_unrelated* ~/ADAMS/genotypes/QMUL_Aug_23/outputs/
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
for ancestry in sas eur afr;
	do
		~/plink2 \
		--bfile ./outputs/combined_adams_imputed_qc_unrelated \
		--keep ./outputs/$ancestry\_iids \
		--maf 0.05 \
		--geno 0.01 \
		--indep-pairwise 1000 100 0.01 \
		--out ./outputs/$ancestry\_pruned

		# filter
		~/plink2 \
		--bfile ./outputs/combined_adams_imputed_qc_unrelated \
		--keep ./outputs/$ancestry\_iids \
		--extract ./outputs/$ancestry\_pruned.prune.in \
		--pca 2 \
		--out ./outputs/$ancestry\_pcs
	done
````

#### Make covar & pheno files
````R
library(tidyverse)
library(RNOmni)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

ancestry_calls_hgdp_1kg = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_detailed.tsv") %>% dplyr::select(IID,contains("ancestry"))
cov = read_tsv("./pheno/adams_covars.tsv") %>% dplyr::select(-FID)
all_cov = ancestry_calls_hgdp_1kg %>% left_join(cov,by="IID")
pheno = read_tsv("./pheno/adams_pheno.tsv") %>% dplyr::select(-FID)

# loop through each ancestry to find outliers & process phenotype data

# find outliers
find_outliers_process_pheno = function(ancestry,sd_num=3){

	pcs = read_table(paste0("./outputs/",ancestry,"_pcs.eigenvec")) %>% dplyr::select(-1)

	mean_pc1 = mean(pcs$PC1)
	mean_pc2 = mean(pcs$PC2)
	sd_pc1 = sd(pcs$PC1)
	sd_pc2 = sd(pcs$PC2)
	upper_pc1 = mean_pc1 + sd_num*sd_pc1
	upper_pc2 = mean_pc2 + sd_num*sd_pc2
	lower_pc1 = mean_pc1 - sd_num*sd_pc1
	lower_pc2 = mean_pc2 - sd_num*sd_pc2

	pcs = pcs %>%
		mutate(outlier = ifelse(PC1 > upper_pc1 | PC1 < lower_pc1 | PC2 > upper_pc2 | PC2 < lower_pc2,
			"outlier",
			"keep"))
	message("PCA outliers:")
	print(pcs %>% dplyr::count(outlier))		
	# save outliers to file
	outliers = pcs %>% filter(outlier=="outlier") %>% dplyr::select(IID) %>% mutate(FID = IID)
	write_tsv(outliers,paste0("./outputs/pca_outliers_",ancestry,".tsv"))

	# combine with main covar file
	pcs = pcs  %>% left_join(all_cov,by="IID")
	p = ggplot(pcs %>% filter(outlier=="keep"),aes(PC1,PC2,fill = bigsnpr_ancestry))+
		geom_point(size=3,data = pcs %>% filter(outlier!="keep"),alpha=0.5,shape=21)+
		geom_point(size=3,color="black",shape=21)+
		theme_bw()+
		scale_fill_brewer(palette="Paired")+
		labs(fill="Ancestry group")+
		ggtitle(toupper(ancestry))

	png(paste0("./outputs/pca_outliers_",ancestry,".png"),res=900,units="in",width=6,height=4)
	print(p)
	dev.off()

	# write covar file
	pcs = pcs %>% filter(!is.na(PC1)) %>% dplyr::select(-contains("ancestry")) %>% mutate(FID = IID) %>% dplyr::select(FID,IID,contains("age"),contains("sex"),contains("PC")) %>% na.omit()
	write_tsv(pcs,paste0("./pheno/",ancestry,"_covars_with_pcs.tsv"))

	# make pheno file
	pheno = pheno %>%
		filter(IID %in% pcs$IID)

	# join with pcs
	pheno = pheno %>%
		dplyr::select(IID,gARMSS,edss) %>%
		na.omit() %>%
		mutate(FID = IID) %>%
		dplyr::select(FID,IID,gARMSS,edss)

	write_tsv(pheno,paste0("./pheno/",ancestry,"_pheno.tsv"))

}

find_outliers_process_pheno(ancestry="sas")
find_outliers_process_pheno(ancestry="afr")
find_outliers_process_pheno(ancestry="eur")

````

#### Run GWAS
````unix
# regenie GWAS
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# remove PCA outliers and filter to MAF 0.05 within each cluster
for ancestry in sas eur afr;
do
	# filter in PLINK
	~/plink2 --bfile ./outputs/combined_adams_imputed_qc_unrelated \
	--keep ./outputs/$ancestry\_iids \
	--remove ./outputs/pca_outliers_$ancestry\.tsv \
	--maf 0.05 \
	--make-bed \
	--out ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38
done

for ancestry in sas eur afr;
do
# step 1 REGENIE
	~/regenie \
--step 1 \
--extract /data/scratch/hmy117/adams_imputed/snps_to_keep_step1_regenie.tsv \
--bed ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile ./pheno/$ancestry\_covars_with_pcs.tsv \
--phenoFile ./pheno/$ancestry\_pheno.tsv \
--bsize 1000 \
--qt --lowmem --apply-rint \
--lowmem-prefix tmp_rg \
--out step1_$ancestry\_fit \
--catCovarList Sex

# step 2 REGENIE
	~/regenie \
--step 2 \
--bed ./outputs/filtered_genotypes_for_sev_gwas_$ancestry\_hg38 \
--covarFile ./pheno/$ancestry\_covars_with_pcs.tsv \
--phenoFile ./pheno/$ancestry\_pheno.tsv \
--bsize 1000 \
--qt --lowmem --apply-rint \
--lowmem-prefix tmp_rg \
--pred step1_$ancestry\_fit_pred.list \
--catCovarList Sex \
--out ./outputs/sev_gwas_$ancestry

done
````


#### VEP
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

# print hits at P < 1e-3
rm ./outputs/vep_input
for ancestry in sas eur afr
do
	awk 'NR>1{if($12>3) print $1,$2,$2,$4"/"$5,"+",$3}' ./outputs/sev_gwas_$ancestry\_gARMSS_INT.regenie >> ./outputs/vep_input
done

# annotate all hits with P <1e-3
module load ensembl-vep
~/ensembl-vep/vep -i ./outputs/vep_input -o ./outputs/snp_annotations \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--check_existing \
--most_severe \
--tab --fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,Consequence,Existing_variation"

~/ensembl-vep/vep -i ./outputs/vep_input -o ./outputs/snp_annotations_nearest \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence,Existing_variation"


````

#### Combine with IMSGC severity GWAS
##### Liftover to hg19 (for UKB)
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
##### Liftover to hg19 (for UKB)
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

# write to file
write_tsv(imsgc_sev %>% dplyr::select(-POS) %>% dplyr::rename("POS" = BP_hg38),"imsgc_mssev_discovery_hg38.tsv")
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

# read in gwas
all_res = purrr::map(ancestries,function(x){
	gwas_res = read_regenie_gwas(
		paste0("./outputs/sev_gwas_",x,"_gARMSS.regenie")) %>%
	mutate(anc = x)
	gwas_res
})

# combine
all_res = do.call("bind_rows",all_res)

# read in edss gwas
all_res_edss = purrr::map(ancestries,function(x){
	gwas_res = read_regenie_gwas(
		paste0("./outputs/sev_gwas_",x,"_edss.regenie")) %>%
	mutate(anc = x)
	gwas_res
})

# combine
all_res_edss = do.call("bind_rows",all_res_edss)

# compare armss gwas vs edss
all_res_edss_and_armss = all_res %>% left_join(all_res_edss, by = c("CHR","BP","SNP","ALLELE0","ALLELE1","anc"))
sampled_dat = all_res_edss_and_armss %>% filter(P.x < 0.05 | P.y < 0.05)

# split by ancestry
eur_adams = all_res %>% filter(anc == "eur")
afr_adams = all_res %>% filter(anc == "afr")
sas_adams = all_res %>% filter(anc == "sas")

# add annotations
anno = read_table("./outputs/snp_annotations_nearest",skip=31,col_types = cols(.default="c")) %>%
	dplyr::select(1,NEAREST)
colnames(anno) = c("SNP","Gene")
anno = anno %>% distinct()

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
	png(paste0("./outputs/qq_plot_",ancestry,".png"),res=900,units="in",width=4,height=4)
	print(p)
	dev.off()

}
sapply(ancestries,make_qq_plot)

# print all together
library(gridExtra)
png("./outputs/all_qq_plots.png",res=900,units="in",width=4,height=8)
cowplot::plot_grid(plotlist = qqplots,align="v",ncol=1)
dev.off()


# manhattans
manhattans =list()
make_manhattans = function(ancestry){
	gwas_res = all_res %>% filter(anc == ancestry)

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
			group_by(CHR) %>%
			slice_min(P,with_ties=F) %>%
			distinct(Gene,.keep_all=T),
		mapping = aes(label = Gene),min.segment.length = 0,color="black",nudge_y=1,direction = "x")+
		labs(x="Genomic position (hg38)",y="-log10(P)")+
		ggtitle(paste0("GWAS of gARMSS in ",toupper(ancestry)))+
		theme_bw()+
		theme(legend.position="none")

	manhattans[[length(manhattans)+1]] <<- p

	png(paste0("./outputs/manhattan_plot_garmss",ancestry,".png"),res=900,units="in",width=10,height=4)
	print(p)
	dev.off()		
}

sapply(ancestries,make_manhattans)

# print all together
png("./outputs/all_manhattans_plots.png",res=900,units="in",width=10,height=10)
cowplot::plot_grid(plotlist = manhattans,align="v",ncol=1)
dev.off()

# compare with IMSGC GWAS
imsgc_hg38 = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_mssev_discovery_hg38.tsv")

# join on chr:pos
all_res_with_imsgc = imsgc_hg38 %>%
		dplyr::rename("BP"=POS,"ALLELE1" = A1,"ALLELE0" = A2,"A1FREQ" = AF1) %>%
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


````

## Susceptibility GWAS


### Liftover to hg19 (for UKB)
````unix

# lift over to hg19
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/

awk '{print "chr"$1,$4-1,$4,$2}' combined_adams_imputed_qc_unrelated.bim > hg38_bedfile

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

### Merge with UKB
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
		"compatible","incompatible")) %>%
		filter(compatible == "compatible")

# write snp list to keep
ukb_snps_to_keep = combo_snps %>% dplyr::select(ID.y)
ukb_snps_names_to_update = combo_snps %>% dplyr::select(ID.y,ID.x)

write_tsv(ukb_snps_to_keep,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_keep.tsv",col_names = F)
write_tsv(ukb_snps_names_to_update,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_snps_to_update_names.tsv",col_names = F)

````

#### Filter to compatible SNPs

````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_ukb_files.sh
````

#### Merge UKB files and ADAMS genotypes across chromosomes
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_ukb_files.sh
````

#### Merge all UKB-ADAMS chroms together
````unix

# clean up
for i in {1..22}; do rm /data/scratch/hmy117/ukb_chr$i\_cpra.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\_hg19_cpra.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\_hg19_cpra_nodups.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/ukb_chr$i\_filtered_tmp.* ; done
for i in {1..22}; do rm /data/scratch/hmy117/filtered_ukb_chr$i\.* ; done

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_all_chroms_adams_ukb.sh
````

#### Joint PCA on combined dataset
````unix

~/plink --bfile /data/scratch/hmy117/ukb_adams_merged_genotypes_all_chroms \
--maf 0.05 \
--hwe 1e-5 \
--geno 0.1 \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca \
--make-bed

~/plink --bfile /data/scratch/hmy117/adams_ukb_snps_for_pca \
--indep-pairwise 1000 500 0.001 \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned


# restrict to pruned snps
~/plink --bfile /data/scratch/hmy117/adams_ukb_snps_for_pca \
--extract /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned.prune.in \
--make-bed \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned

# calculate PCs
~/plink2 --bfile /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned \
--pca approx \
--out /data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_pca
````

#### Identify ancestry groupings
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

# read in ukb-adams pcs
ukb_adams_pcs = read_tsv("/data/scratch/hmy117/adams_ukb_snps_for_pca_pruned_pca.eigenvec")

# read in adams ancestry calls
ancestry_calls_hgdp_1kg = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_detailed.tsv") %>% dplyr::select(IID,contains("ancestry"))
cov = read_tsv("./pheno/adams_covars.tsv") %>% dplyr::select(-FID)
all_cov = ancestry_calls_hgdp_1kg %>% left_join(cov,by="IID")

# merge
ukb_adams_pcs = ukb_adams_pcs %>%
	mutate(study = ifelse(IID %in% all_cov$IID,"ADAMS","UKB")) %>%
	left_join(all_cov,by="IID")

# read ukb pheno
# get phenotype data for UKB
ukb_pheno = read_csv("/data/Wolfson-PNU-dementia/UKB/PRS_April_2023/phenodata/PRS_59138r672130.csv") %>%
	dplyr::select(EID,sex_f31_0_0,age_at_recruitment_f21022_0_0,genetic_ethnic_grouping_f22006_0_0,contains("g35"))

# bring in ethnicity data
pheno = readRDS("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/78867r672482/78867r672482_FO.rds")
pheno = pheno %>% tibble %>% dplyr::select(eid,contains("ethnic"))
id_bridge = read_csv("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/59138_78867/Bridge_eids_59138_78867.csv")
pheno = pheno %>%
  dplyr::rename("eid_78867" = eid) %>%
  left_join(id_bridge,by="eid_78867") %>%
  dplyr::rename("EID" = eid_59138)
pheno = pheno %>% dplyr::select(EID,contains("ethnic"))
ukb_pheno = ukb_pheno %>% left_join(pheno,by="EID")


# loop through each ancestry to find outliers & process phenotype data

# find outliers
find_outliers_process_pheno = function(ancestry,sd_num=2,out_ancestry){

	# first remove case outliers (derived earlier from within-ancestry PCA)
	infile = paste0("./outputs/pca_outliers_",out_ancestry,".tsv")
	outliers_this_anc = read_tsv(infile)

	# restrict to ADAMS cases who are non-outliers for this ancestry
	pcs = ukb_adams_pcs %>%
		filter(study=="ADAMS" & predicted_ancestry == ancestry & !(IID %in% outliers_this_anc$IID))

	# get pc limits
	mean_pc1 = mean(pcs$PC1)
	mean_pc2 = mean(pcs$PC2)
	sd_pc1 = sd(pcs$PC1)
	sd_pc2 = sd(pcs$PC2)
	upper_pc1 = mean_pc1 + sd_num*sd_pc1
	upper_pc2 = mean_pc2 + sd_num*sd_pc2
	lower_pc1 = mean_pc1 - sd_num*sd_pc1
	lower_pc2 = mean_pc2 - sd_num*sd_pc2

	ukb_adams_pcs = ukb_adams_pcs %>%
		mutate(outlier = ifelse(
			(study == "ADAMS" & (IID %in% outliers_this_anc$IID | predicted_ancestry != ancestry)) |
			(study == "UKB" & (PC1 > upper_pc1 | PC1 < lower_pc1 | PC2 > upper_pc2 | PC2 < lower_pc2)),
			"outlier",
			"keep"))
	message("PCA outliers:")
	print(ukb_adams_pcs %>% dplyr::count(study,outlier))		

	# save keepers to file
	non_outliers = ukb_adams_pcs %>% filter(outlier!="outlier") %>% dplyr::select(IID) %>% mutate(FID = IID)
	write_tsv(non_outliers,paste0("./outputs/ukb_pca_non_outliers_",out_ancestry,".tsv"))


	# plot pcs with and without outliers
	p = ggplot(ukb_adams_pcs,aes(PC1,PC2,fill = outlier))+
		geom_point(size=3,color="black",shape=21)+
		theme_bw()+
		scale_fill_brewer(palette="Set1")+
		labs(fill="Ancestry group")+
		ggtitle(toupper(ancestry))
	png(paste0("./outputs/ukb_adams_pca_outliers_",ancestry,".png"),res=900,units="in",width=6,height=4)
	print(p)
	dev.off()

	# plot just the people to keep
	p1 = ggplot(ukb_adams_pcs %>% filter(outlier =="keep"),aes(PC1,PC2,fill = bigsnpr_ancestry))+
		geom_point(size=3,color="black",shape=21)+
		theme_bw()+
		facet_wrap(~study)+
		scale_fill_brewer(palette="Set1")+
		labs(fill="Ancestry group")+
		ggtitle(toupper(ancestry))


	# combine with main covar file
	ukb_adams_pcs = ukb_adams_pcs %>%
		filter(outlier=="keep") %>%
		dplyr::select(-bigsnpr_ancestry,-predicted_ancestry,-outlier)

	# filter to this ancestry
	ukb_pheno_this_ancestry = ukb_pheno %>%
		filter(EID %in% ukb_adams_pcs$IID) %>%
		dplyr::rename("IID" = EID,"Sex" = sex_f31_0_0,"age_at_recruitment" = age_at_recruitment_f21022_0_0)


	# split covar / pheno into UKB and non-UKB
	ukb_pheno_this_ancestry_tmp = ukb_adams_pcs %>%
		filter(study =="UKB") %>%
		dplyr::select(-Sex,-age_at_recruitment) %>%
		left_join(ukb_pheno_this_ancestry,by="IID")

	ukb_pheno_this_ancestry_tmp = ukb_pheno_this_ancestry_tmp	%>%
		mutate(MS_status = ifelse(!is.na(ukb_pheno_this_ancestry_tmp$source_of_report_of_g35_multiple_sclerosis_f131043_0_0),2,1))

	adams_pheno_tmp = ukb_adams_pcs %>%
		filter(study == "ADAMS") %>%
		mutate(MS_status = 2)

	# recombine
	all_pheno_cov = bind_rows(adams_pheno_tmp,ukb_pheno_this_ancestry_tmp)


	# plot just the people to keep
	ethnicity_counts = all_pheno_cov %>% dplyr::count(ethnic_background_f21000_0_0) %>% filter(n>100)

	p2 = ggplot(all_pheno_cov %>% filter(
		(study=="UKB" & ethnic_background_f21000_0_0 %in% ethnicity_counts$ethnic_background_f21000_0_0) |
		study == "ADAMS"),
		aes(PC1,PC2,fill = ethnic_background_f21000_0_0))+
		geom_point(size=3,color="black",shape=21)+
		theme_bw()+
		facet_wrap(~study)+
		labs(fill="Self-reported ethnicity (UKB)")+
		ggtitle(toupper(ancestry))
	png(paste0("./outputs/ukb_adams_pca_outliers_with_ethnicity_ancestry_",ancestry,".png"),res=900,units="in",width=6,height=8)
	print(cowplot::plot_grid(p1,p2,align="v",ncol=1))
	dev.off()

	# write covar file
	all_cov = all_pheno_cov %>%
		filter(!is.na(PC1)) %>%
		dplyr::select(-contains("genetic"),-contains("g35")) %>%
		mutate(FID = IID) %>%
		dplyr::select(FID,IID,contains("age"),contains("sex"),contains("PC")) %>%
		na.omit()


	write_tsv(all_cov,paste0("./pheno/susceptibility_",out_ancestry,"_covars_with_pcs.tsv"))

	# make pheno file
	all_pheno = all_pheno_cov %>%
		mutate(FID = IID) %>%
		dplyr::select(FID,IID,MS_status) %>%
		na.omit()

	write_tsv(all_pheno,paste0("./pheno/susceptibility_",out_ancestry,"_pheno.tsv"))

	all_pheno_cov %>% dplyr::count(MS_status,study) %>% print()


}

find_outliers_process_pheno(ancestry="CSA",out_ancestry="sas")
find_outliers_process_pheno(ancestry="AFR",out_ancestry="afr")
find_outliers_process_pheno(ancestry="EUR",out_ancestry="eur")

````

#### GWAS
##### Split
````unix
view /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_diff_miss_snps.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_diff_miss_snps.sh

````

#### Run GWAS
````unix

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas_susceptibility.sh

````

# get rid of palindromes 
#########################################
# PICK UP HERE
# FIX tmp file REGENIE
# THINK ABOUT
# - BETTER PCA INFERENCE
# - other causes inflation in case-control
# - downstream analysis for severity
# - recalc PCs within clusters?
# - redo VEP with hg19 (liftover)
#########################################


#### VEP
````unix
# prepare snps for vep
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

module load R/4.2.2


# liftover
for ancestry in sas eur afr
do
	Rscript ./scripts/liftover_hg19_to_hg38.R ./outputs/susceptibility_gwas_$ancestry\_MS_status.regenie
	head ./outputs/susceptibility_gwas_$ancestry\_MS_status.regenie_hg38
	awk 'NR>1{if($11>8) print $1,$13,$13,$3"/"$4,"+",$2}' ./outputs/susceptibility_gwas_$ancestry\_MS_status.regenie_hg38 > ./outputs/susceptibility_vep_input_$ancestry

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

# read in gwas
all_res = purrr::map(ancestries,function(x){
	gwas_res = read_regenie_gwas(
		paste0("./outputs/susceptibility_gwas_",x,"_MS_status.regenie_hg38")) %>%
	mutate(anc = x)
	gwas_res
})

# combine
all_res = do.call("bind_rows",all_res)

# manhattans
make_manhattans = function(ancestry){
	gwas_res = all_res %>% filter(anc == ancestry)

	# read annotations
	anno = read_table(paste0("./outputs/snp_annotations_nearest_susceptibility_",ancestry),skip=31,col_types = cols(.default="c")) %>%
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
			P < 5e-8 ~ "sig",
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
			filter(P < 5e-8 & SNP %in% anno$SNP) %>%
			left_join(anno,by="SNP") %>%
			group_by(CHR) %>%
			slice_min(P,with_ties=F) %>%
			distinct(Gene,.keep_all=T),
		mapping = aes(label = Gene),min.segment.length = 0,color="black",nudge_y=1,direction = "x")+
		labs(x="Genomic position (hg38)",y="-log10(P)")+
		ggtitle(paste0("GWAS of MS risk in ",toupper(ancestry)))+
		theme_bw()+
		theme(legend.position="none")


	png(paste0("./outputs/manhattan_plot_susceptibility",ancestry,".png"),res=900,units="in",width=10,height=4)
	print(p)
	dev.off()		
}

sapply(ancestries,make_manhattans)

# print all together
png("./outputs/all_manhattans_plots.png",res=900,units="in",width=10,height=10)
cowplot::plot_grid(plotlist = manhattans,align="v",ncol=1)
dev.off()

# compare with IMSGC GWAS
imsgc_hg38 = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_mssev_discovery_hg38.tsv")

# join on chr:pos
all_res_with_imsgc = imsgc_hg38 %>%
		dplyr::rename("BP"=POS,"ALLELE1" = A1,"ALLELE0" = A2,"A1FREQ" = AF1) %>%
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
make_forest("1:200258565")

````
