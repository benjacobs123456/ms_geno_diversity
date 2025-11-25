
# Go to directory
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
module load plink/1.9-170906
````
#########################
# IMPUTATION PREP
#########################

# Step 1
- Conversion to PLINK1 binary
- Basic QC pre-imputation
- remove controls
````unix
qsub ./scripts/convert_to_plink_binary.sh
````
# Step 1.1
- Ancestry inference
````unix
qsub ./scripts/download_hgdp_1kg.sh
```

# filter to ADAMS variants
```unix
module load R/3.6.1
Rscript ./scripts/get_adams_vars.R
qsub ./scripts/filter_hgdp_to_adams_vars.sh

# download ancestry calls
~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/data_intersection/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv \
/data/scratch/hmy117/hgdp_1kg_genomes/

# merge hgdp-1kg
cd /data/scratch/hmy117/hgdp_1kg_genomes/
for i in {2..22};
  do
    echo filtered_1kg_hgdp_1kg_chr$i >> merge_filelist
  done

module load plink/1.9-170906
cd /data/scratch/hmy117/hgdp_1kg_genomes/

# try merge
plink --bfile filtered_1kg_hgdp_1kg_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg

# remove duplicates & do QC
for i in {1..22};
  do
    plink --bfile filtered_1kg_hgdp_1kg_chr$i \
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
plink --bfile hgdp_1kg_nodups_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg

# do same QC for adams
plink --bfile ~/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno \
--maf 0.05 \
--geno 0.01 \
--mind 0.1 \
--hwe 1e-10 \
--extract combined_hgdp_1kg.bim \
--out ADAMS_qc_genos_for_pcs \
--make-bed

# missingness filter
plink --bfile combined_hgdp_1kg \
--mind 0.1 \
--out combined_hgdp_1kg_nonmissing \
--make-bed

# get intersection of nonpalindromic, compatible snps
````R
library(tidyverse)
adams_snps = read_table2("ADAMS_qc_genos_for_pcs.bim",col_names=F)
kg_snps = read_table2("combined_hgdp_1kg_nonmissing.bim",col_names=F)

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
```

````unix
# filter to compatible snps
# remove related indivs
~/king -b combined_hgdp_1kg_nonmissing.bed --degree 3
awk 'NR>1{print $1,$2}' king.kin0 > indivs_to_remove
plink --bfile combined_hgdp_1kg_nonmissing \
--extract snps_to_keep.tsv \
--remove indivs_to_remove \
--out combined_hgdp_1kg_unrelated \
--make-bed

# prune 1kg
plink --bfile combined_hgdp_1kg_unrelated \
--indep-pairwise 1000 100 0.01 \
--out pruned_snps_hgdp_1kg

# filter
plink --bfile combined_hgdp_1kg_unrelated \
--out pruned_hgdp_1kg \
--extract pruned_snps_hgdp_1kg.prune.in \
--make-bed


# calculate PCs
~/plink2 --bfile pruned_hgdp_1kg \
--pca allele-wts 50 \
--freq \
--out hgdp_kg_pcs

# project original dataset
~/plink2 --bfile pruned_hgdp_1kg \
--read-freq hgdp_kg_pcs.afreq \
--score hgdp_kg_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize \
--score-col-nums 6-45 \
--out hgdp_kg_pcs_rescored

# project ADAMS samples
~/plink2 --bfile ADAMS_qc_genos_for_pcs \
--read-freq hgdp_kg_pcs.afreq \
--score hgdp_kg_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize \
--score-col-nums 6-45 \
--out adams_pcs
````

# Download metadata
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv

# Prepare for ADMIXTURE
cd /data/scratch/hmy117/hgdp_1kg_genomes
module load plink/1.9-170906
plink --bfile pruned_hgdp_1kg --bmerge ADAMS_qc_genos_for_pcs --out combined_adams_kg_for_admixture --make-bed
plink --bfile combined_adams_kg_for_admixture --geno 0.01 --out combined_adams_kg_for_admixture_qc --make-bed

````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")
adams = read_table2("adams_pcs.sscore")
kg_hgdp = read_table2("hgdp_kg_pcs_rescored.sscore")
meta = read_tsv("gnomad_meta_v1.tsv") %>%
  dplyr::select(2,hgdp_tgp_meta.Population,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("IID","pop","superpop")
kg_hgdp = kg_hgdp %>%
left_join(meta,by="IID")

# prep for admixture
fam_file = read_table2("combined_adams_kg_for_admixture_qc.fam",col_names = F)
colnames(fam_file)[2] = "IID"
fam_file = fam_file %>%
  left_join(kg_hgdp %>% dplyr::select(IID,superpop),by="IID") %>%
  dplyr::select(IID,superpop) %>%
  replace(is.na(.),"-") %>%
  dplyr::select(superpop)
write_tsv(fam_file,"combined_adams_kg_for_admixture_qc.pop",col_names=F)  

# raw plots
ggplot(kg_hgdp,aes(PC1_AVG,PC2_AVG,col=superpop))+geom_point()
ggplot(kg_hgdp,aes(PC1_AVG,PC2_AVG,col=pop))+geom_point()


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
table(adams$predicted_ancestry)
p=ggplot(adams,aes(PC1_AVG,PC2_AVG,col=predicted_ancestry))+
  geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Inferred genetic ancestry")

p2=ggplot(adams,aes(PC1_AVG,PC2_AVG,col=predicted_ancestry))+
  geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Inferred genetic ancestry")+
  ggtitle("ADAMS")+
  theme(legend.position = c(0.9,0.9))



p3 = ggplot(kg_hgdp,aes(PC1_AVG,PC2_AVG,col=superpop))+
  geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Genetic ancestry")+
  ggtitle("1kg + HGDP reference")+
  theme(legend.position = c(0.9,0.9))
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pc_plot_with_ref.png",res=900,units="in",width=8,height=12)
gridExtra::grid.arrange(p2,p3,nrow=2)
dev.off()


png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pc_plot.png",res=900,units="in",width=6,height=4)
p
dev.off()


# read in pheno file
pheno = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_pheno.tsv", col_types = cols(.default = "c"))

# merge
adams = adams %>% left_join(pheno,by="IID")

# code NAs as other
adams = adams %>%
  mutate(Ethnicity = ifelse(is.na(Ethnicity),"Other ethnic group",Ethnicity))

# plot vs self-reported ethnicity
p=ggplot(adams,aes(PC1_AVG,PC2_AVG,col=Ethnicity))+
  geom_point()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  labs(x="PC1",y="PC2",col="Self-reported ethnicity")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pc_plot_vs_ethnicity.png",res=900,units="in",width=6,height=4)
p
dev.off()


# write list of sas
sas = adams %>%
  filter(predicted_ancestry=="CSA") %>%
  dplyr::select(1,2)
write_tsv(sas,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/sas_iids.tsv")

# write all ancestry calls
anc_calls = adams %>%
  dplyr::select(1,2,predicted_ancestry)
write_tsv(anc_calls,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls.tsv")

````


````
library(tidyverse)

# modify pheno file
fam_file = read_table2("./outputs/ADAMS_geno.fam", col_types = cols(.default = "c"),col_names=F)
armss = read_tsv("./pheno/adams_armss.tsv",col_names = F, col_types = cols(.default = "c"))
covars = read_tsv("./pheno/adams_covars.tsv",col_names = F, col_types = cols(.default = "c"))

# filtering
armss = armss %>%
  filter(X1 %in% fam_file$X2) %>%
  dplyr::rename("IID" = X1) %>%
  left_join(fam_file %>% dplyr::rename("IID" = X2), by = "IID") %>%
  dplyr::select(X1,IID,X2) %>%
  distinct(IID,.keep_all=T) %>%
  mutate(armss_z = RNOmni::rankNorm(X2)) %>%
  dplyr::select(-X2) %>%
  mutate(X1 = IID)


covars = covars %>%
  filter(X1 %in% fam_file$X2) %>%
  dplyr::rename("IID" = X1) %>%
  left_join(fam_file %>% dplyr::rename("IID" = X2,"FID"=X1), by = "IID") %>%
  dplyr::select(FID,IID,X2,X5) %>%
  distinct(IID,.keep_all=T) %>%
  dplyr::rename("age" = X2,"sex" = X5) %>%
  mutate(FID=IID)


# rank norm
write_tsv(armss,"./pheno/adams_armss_modified.tsv")
write_tsv(covars,"./pheno/adams_covars_modified.tsv")

````


# Run ADMIXTURE
~/dist/admixture_linux-1.3.0/admixture combined_adams_kg_for_admixture_qc.bed 7 \
--supervised \
-j10

# Inspect ADMIXTURE results
````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes/")
q_inference = read_table2("combined_adams_kg_for_admixture_qc.7.Q",col_names=F)
pop_file = read_table2("combined_adams_kg_for_admixture_qc.pop",col_names=F)
colnames(pop_file)[1] = "pop"
fam_file = read_table2("combined_adams_kg_for_admixture_qc.fam",col_names = F) %>% dplyr::select(2) %>%
dplyr::rename("IID" = X2)

df = bind_cols(pop_file,q_inference,fam_file)
codex = df %>% pivot_longer(cols = contains("X")) %>%
  group_by(pop,name) %>%
  summarise(prop = mean(value)) %>%
  slice_max(prop) %>%
  filter(pop != "-") %>%
  dplyr::select(-prop)

df = df %>%
  pivot_longer(cols = contains("X")) %>%
  filter(pop == "-") %>%
  left_join(codex,by="name")


iid_order = df %>% filter(pop.y=="CSA") %>% arrange(value)
df$IID = factor(df$IID,levels = unique(iid_order$IID),ordered=T)

# find dominant pop

p=ggplot(df,aes(IID,value,fill=pop.y))+
  geom_col(color="black",size=0.05)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank() )+
  labs(x="Sample",y="Proportion of each ancestry",fill="Superpopulation")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/admixture_plot.png",res=900,units="in",width=10,height=4)
p
dev.off()

````

# Step 2
- Export as VCFs

````unix
# flip strand
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
module load plink/1.9-170906
for i in {1..22};
	do
		plink --bfile ./outputs/ADAMS_geno \
		--chr $i \
		--recode vcf \
		--output-chr chrMT \
		--out ./imputation_raw_files/chr$i &
	done
````

# Step 3
- Sort with bcf tools
````unix
module load bcftools
for i in {1..22};
	do
		bcftools sort ./imputation_raw_files/chr$i\.vcf \
		 -Oz -o ./imputation_raw_files/sorted_chr$i\.vcf.gz &
	done
````

# Step 4
- Attempt imputation at TOPMED server (TOPMED-R2)
- This will fail due to strand flips
- Download excluded snps and use strand flips to flip
# wget https://imputation.biodatacatalyst.nhlbi.nih.gov/results/job-20230901-122606-755/statisticDir/snps-excluded.txt

````unix
awk 'FS="\t"{if($2=="Strand flip") print $1}' ./imputation_raw_files/snps-excluded.txt > ./imputation_raw_files/snps_to_flip
awk 'FS="\t"{if($2!="Strand flip") print $1}' ./imputation_raw_files/snps-excluded.txt | uniq > ./imputation_raw_files/snps_to_exclude

# set var IDs to chr:pos:ref:alt
for i in {1..22};
	do
		module unload plink
		module load plink/2.0-20220603
		plink2 --vcf ./imputation_raw_files/chr$i\.vcf \
		--set-all-var-ids chr@:#:\$r\:\$a \
		--make-bed \
		--out ./imputation_raw_files/updated_ids_chr$i &
	done


for i in {1..22};
	do
		module unload plink
		module load plink/1.9-170906
		plink --bfile ./imputation_raw_files/updated_ids_chr$i \
		--exclude ./imputation_raw_files/snps_to_exclude \
		--out ./imputation_raw_files/flipped_chr$i \
		--output-chr chrMT \
		--recode vcf &
	done

# flip
	for i in {1..22};
		do
			module unload plink
			module load plink/1.9-170906
			plink --vcf ./imputation_raw_files/flipped_chr$i\.vcf --double-id \
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

# Step 5
- Performed imputation again
- Download imputed data
````unix
mkdir /data/scratch/hmy117/adams_imputed
cd /data/scratch/hmy117/adams_imputed
qsub ~/ADAMS/genotypes/QMUL_Aug_23/download_imputed_data.sh


# unzip
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/unzip.sh
````

#########################
# PROCESS IMPUTED DATA
#########################

# Step 1  
- Check imputation quality
- Examine snps in R
- Plot info in MAF bins
````unix
Rscript ~/ADAMS/genotypes/QMUL_Aug_23/scripts/check_imputation_quality.R
````

# Step 2
- Filter to INFO > 0.5 & MAF > 0.01

````unix
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_imputed_data.sh
````

# Step 3
- SNP QC on individual VCFs
- Conversion back to plink
````unix
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_snp_qc_imputed.sh
````

# Step 4
- Modify IDs in R
````R
library(tidyverse)
for(i in c(1:22)){
df = read_table2(paste0("ADAMS_imputed_qc_chr",i,".fam"),col_names=F) %>%
	mutate(X1 = str_remove(X1,"^0_")) %>%
	tidyr::separate(X1,sep = "_",into = c("part1","oragene")) %>%
	dplyr::select(X2,oragene) %>%
	mutate(oldfid = X2, oldiid = X2, newfid = oragene, newiid = oragene) %>%
	dplyr::select(oldfid,oldiid,newfid,newiid)
write_tsv(df,paste0("chr",i,"_newids.tsv"),col_names = F)
}
````
- Merge across chromosomes in PLINK
- Rename sample IDs
````unix
# rename IDs
for i in {22..1}; do plink --bfile ADAMS_imputed_qc_chr$i --update-ids chr$i\_newids.tsv --out chr$i\_combined_adams_imputed_newids --make-bed; done

# Merge
cd /data/scratch/hmy117/adams_imputed/
rm filelist_for_merge
for i in {2..22}; do echo chr$i\_combined_adams_imputed_newids >> filelist_for_merge; done
module load plink/1.9-170906
plink --bfile chr1_combined_adams_imputed_newids \
--merge-list filelist_for_merge \
--out combined_adams_imputed \
--make-bed
````
# QC
cd /data/scratch/hmy117/adams_imputed/

plink --bfile combined_adams_imputed \
--make-bed \
--mind 0.1 \
--maf 0.01 \
--hwe 1e-15 \
--geno 0.1 \
--out combined_adams_imputed_qc \
--snps-only \
--chr 1-22

# Step 5 - liftover to hg19
- Copy files
````unix

# lift over to hg19
awk '{print "chr"$1,$4-1,$4,$2}' combined_adams_imputed_qc.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$1":"$3":"$4}' hg19_bedfile > hg19_snps
awk '{print $1":"$3":"$4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions and IDs
plink --bfile combined_adams_imputed_qc \
--update-name hg19_snps \
--make-bed \
--out combined_adams_imputed_qc_temp_hg19

plink --bfile combined_adams_imputed_qc_temp_hg19 \
--update-map hg19_snp_positions \
--make-bed \
--out combined_adams_imputed_qc_hg19

~/plink2 --bfile combined_adams_imputed_qc_hg19 \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--out combined_adams_imputed_qc_hg19_chrpos

# rm dups
~/plink2 --bfile combined_adams_imputed_qc_hg19_chrpos \
--rm-dup \
--make-bed

~/plink2 --bfile combined_adams_imputed_qc_hg19_chrpos \
--exclude plink2.rmdup.mismatch \
--make-bed \
--out combined_adams_imputed_qc_hg19_chrpos_nodups




cp combined_adams_imputed_qc_hg19_chrpos_nodups* /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs
````

#########################
# GWAS
#########################

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/


# regression on EDSS - SAS
awk 'NR>1{if($3=="CSA") print $2,$2}' ./outputs/ancestry_calls.tsv | grep -v Control > ./outputs/sas_iids

# get PCs
~/plink2 \
--bfile ./outputs/combined_adams_imputed_qc_hg19_chrpos_nodups \
--keep ./outputs/sas_iids \
--maf 0.05 \
--geno 0.01 \
--indep-pairwise 1000 100 0.01 \
--out sas_pruned

# filter
~/plink2 \
--bfile ./outputs/combined_adams_imputed_qc_hg19_chrpos_nodups \
--keep ./outputs/sas_iids \
--extract sas_pruned.prune.in \
--pca \
--out ./outputs/sas_pcs

# join with covars
```R
library(tidyverse)
pcs = read_table("sas_pcs.eigenvec") %>% dplyr::select(-1)
cov = read_tsv("./pheno/adams_covars.tsv")

cov = cov %>% left_join(pcs,by="IID") %>% filter(!is.na(PC1))
write_tsv(cov,"./pheno/sas_covars_with_pcs.tsv")

```

# run GWAS - SAS
~/plink2 --bfile ./outputs/combined_adams_imputed_qc_hg19_chrpos \
--keep ./outputs/sas_iids \
--pheno ./pheno/adams_pheno.tsv \
--glm hide-covar \
--covar ./pheno/sas_covars_with_pcs.tsv \
--maf 0.01

```R
library(tidyverse)

# read in GWAS
read_gwas = function(x){
  read_tsv(x, col_types = cols_only(
  `#CHROM` = col_double(),
  POS = col_double(),
  ID = col_character(),
  REF = col_character(),
  ALT = col_character(),
  A1 = col_character(),
  TEST = col_character(),
  OBS_CT = col_double(),
  BETA = col_double(),
  SE = col_double(),
  P = col_double()
)) %>%
  filter(TEST == "ADD" & !is.na(P)) %>%
  dplyr::rename("CHR" = `#CHROM`, "BP" = POS, "SNP" = ID) %>%
  dplyr::select(-TEST)
}
sas = read_gwas("plink2.gARMSS.glm.linear")%>% mutate(ancestry = "SAS")

qqman::manhattan(eur %>% filter(P<0.005))
```




#########################
# PCs
#########################


# Step 1
- Download 1kg genotypes
- Recode to plink
````unix
qsub ~/ADAMS/genotypes/QMUL_Jun_23/scripts/download_1kg.sh
````

# filter 1kg snps

module load plink/2.0-20220603
for i in {1..22}
	do
    plink2 --bfile /data/scratch/hmy117/1kg_chr$i --set-all-var-ids @:# --make-bed --out /data/scratch/hmy117/updated_snp_ids_1kg_chr$i &
	done

# Extract ADAMS SNPs from 1kg data
awk '{print $2}' ~/ADAMS/genotypes/QMUL_Jun_23/outputs/ADAMS_imputed_qc_hg19_chrpos.bim > adams_snps
for i in {1..22}
	do
    plink2 --bfile /data/scratch/hmy117/updated_snp_ids_1kg_chr$i \
    --extract adams_snps \
    --out /data/scratch/hmy117/filtered_1kg_chr$i \
    --make-bed &
  done


# merge
rm merge_filelist
for i in {2..22}; do echo /data/scratch/hmy117/filtered_1kg_chr$i >> merge_filelist; done

module unload plink
module load plink/1.9-170906
plink --bfile /data/scratch/hmy117/filtered_1kg_chr1 \
--merge-list merge_filelist \
--out /data/scratch/hmy117/1kg_merged

# repeat excluding triallelics
for i in {1..22}
	do
    plink --bfile /data/scratch/hmy117/filtered_1kg_chr$i \
    --exclude 1kg_merged.missnp \
    --out /data/scratch/hmy117/unique_filtered_1kg_chr$i \
    --make-bed &
  done

# merge 1kg
rm merge_filelist
for i in {2..22}; do echo /data/scratch/hmy117/unique_filtered_1kg_chr$i >> merge_filelist; done

plink --bfile /data/scratch/hmy117/unique_filtered_1kg_chr1 \
--merge-list merge_filelist \
--out /data/scratch/hmy117/1kg_merged

# aggressive QC with ADAMS SNPs
plink --bfile ~/ADAMS/genotypes/QMUL_Jun_23/outputs/ADAMS_imputed_qc_hg19_chrpos \
--maf 0.05 \
--geno 0.99 \
--indep-pairwise 1000 100 0.05 \
--make-bed

# filter
plink --bfile ~/ADAMS/genotypes/QMUL_Jun_23/outputs/ADAMS_imputed_qc_hg19_chrpos \
--extract plink.prune.in \
--out adams_pruned \
--make-bed

# now merge with ADAMS
cd /data/scratch/hmy117/
plink --bfile adams_pruned \
--bmerge /data/scratch/hmy117/1kg_merged \
--out adams_1kg_merged

# try again without dups
plink --bfile adams_pruned \
--exclude adams_1kg_merged.missnp \
--out adams_for_merging \
--make-bed

plink --bfile adams_for_merging \
--bmerge /data/scratch/hmy117/1kg_merged \
--out adams_1kg_merged

# more QC
plink --bfile adams_1kg_merged \
--geno 0.01 \
--make-bed \
--out merged_genos_pca

# pca
plink --bfile merged_genos_pca \
--pca

````R
library(tidyverse)
setwd("/data/scratch/hmy117/")

pops = read_tsv("/data/home/hmy117/ADAMS/deprecated/rfmix.1kg_reference_superpops.tsv",col_names=F)
pcs = read_table2("plink.eigenvec",col_names=F)

pcs = pcs %>% dplyr::rename("IID" = X2) %>% left_join(pops %>% dplyr::rename("IID"=X1,"pop" = X2),by="IID")

pcs = pcs %>% mutate(cohort = ifelse(is.na(pop),"ADAMS","1000 genomes"))

p = ggplot(pcs,aes(X3,X4,col=pop))+geom_point()+facet_wrap(~cohort)+theme_minimal()+labs(x="PC1",y="PC2")
p1 = ggplot(pcs,aes(X5,X6,col=pop))+geom_point()+facet_wrap(~cohort)+theme_minimal()+labs(x="PC3",y="PC4")
gridExtra::grid.arrange(p,p1)


library(caret)

# code as SAS vs non-sas

kg_pcs = pcs %>% filter(cohort=="1000 genomes") %>%
	dplyr::select(-1,-2,-cohort) %>%
	dplyr::select(1:50,pop)

rf_fit = train(factor(pop) ~ ., data = kg_pcs,method = "ranger")

adams_pcs = pcs %>% filter(cohort=="ADAMS") %>% dplyr::select(-1,-2,-cohort)
predictions = predict(rf_fit, newdata = adams_pcs)
adams_pcs$pop = predictions

setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Jun_23/outputs/")
p1 = ggplot(kg_pcs,aes(X3,X4,col=pop))+geom_point()+theme_minimal()+labs(x="PC1",y="PC2")+
ggtitle("1000 Genomes")
p2 = ggplot(adams_pcs,aes(X3,X4,col=pop))+geom_point()+theme_minimal()+labs(x="PC1",y="PC2")+
ggtitle("ADAMS")
p3 = ggplot(adams_pcs,aes(X5,X6,col=pop))+geom_point()+theme_minimal()+labs(x="PC3",y="PC4")+
ggtitle("ADAMS")

png("pcs.png",res=600,units="in",width=8,height=4)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

write_tsv(adams_pcs,"adams_pops.tsv")
write_tsv(adams_pcs %>% filter(pop=="SAS"),"adams_sas.tsv")
write_tsv(adams_pcs %>% filter(pop=="AFR"),"adams_afr.tsv")
````


# combine with non-eur UKB
# find UKB outliers
````R
# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
setwd("/data/Wolfson-UKBB-Dobson/UKB_alcohol")

# read data
df = read_tsv("/data/Wolfson-UKBB-Dobson/UKB_alcohol/ukb_pheno_21_01_22.tsv")

# read in source of report data
source_of_report_data = read_tsv("/data/Wolfson-UKBB-Dobson/ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  `Date G35 first reported (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

# combine
df = df %>% dplyr::select(-MS_status) %>% left_join(source_of_report_data,by="EID")

# exclude participants who have withdrawn
withdrawn = read_tsv("../helper_progs_and_key/excluded_indivs",col_names=FALSE)
df = df %>% filter(!EID %in% withdrawn$X1)

# restrict to genetic cohort

# Remove non-European participants
eur = df %>% filter(`Genetic ethnic grouping.0.0` =="Caucasian")
non_eur = df %>% filter(!EID %in% eur$EID)
non_eur = non_eur %>% dplyr::select(EID,MS_status,`Ethnic background.0.0`)

# get self-reported SAS
table(non_eur$`Ethnic background.0.0`)
sas_ethnicities = c("Any other Asian background",
"Bangladeshi",
"Indian",
"Pakistani")

sas = non_eur %>%
filter(`Ethnic background.0.0` %in% sas_ethnicities)
write_tsv(sas,"~/ADAMS/ukb/ukb_sas_pheno.tsv")
````


# filter to SAS genotypes
````unix
cd ~/ADAMS/ukb
awk 'NR>1{print $1,$1}' ukb_sas_pheno.tsv > samples_to_keep
qsub process_ukb_sas_genotypes.sh
````

# merge
module unload plink
module load plink/1.9-170906

rm merge_filelist
for i in {2..22};
  do
    echo /data/scratch/hmy117/filtered_ukb_sas_ms_snps_chrpos_chr$i>> merge_filelist
  done

plink --bfile /data/scratch/hmy117/filtered_ukb_sas_ms_snps_chrpos_chr1 \
--merge-list merge_filelist \
--out merged_ukb \
--make-bed


module unload plink
module load plink/2.0-20220603
plink2 --bfile merged_ukb \
--set-all-var-ids @:# \
--make-bed \
--out merged_ukb_hg19_chrpos


# filter to ADAMS SNPs
plink2 --bfile merged_ukb_hg19_chrpos \
--extract ~/ADAMS/genotypes/QMUL_Jan_23/outputs/ADAMS_imputed_qc_hg19_chrpos.bim \
--out filtered_ukb_sas_adams_snps \
--make-bed

# filter to UKB SNPs
plink2 --bfile ~/ADAMS/genotypes/QMUL_Jan_23/outputs/ADAMS_imputed_qc_hg19_chrpos \
--extract filtered_ukb_sas_adams_snps.bim \
--out adams_for_merge \
--make-bed

# now merge ukb with ms
module unload plink
module load plink/1.9-170906
plink --bfile adams_for_merge \
--bmerge filtered_ukb_sas_adams_snps \
--out ms_ukb_merged \
--make-bed

# retry
plink --bfile adams_for_merge \
--exclude ms_ukb_merged-merge.missnp \
--out /data/scratch/hmy117/adams_for_merge_no_triallelics \
--make-bed

plink --bfile filtered_ukb_sas_adams_snps \
--exclude ms_ukb_merged-merge.missnp \
--out /data/scratch/hmy117/ukb_for_merge_no_triallelics \
--make-bed

plink --bfile /data/scratch/hmy117/adams_for_merge_no_triallelics \
--bmerge /data/scratch/hmy117/ukb_for_merge_no_triallelics \
--out /data/scratch/hmy117/ms_ukb_merged \
--make-bed


# filter missing genos
cut -f1,2 /data/scratch/hmy117/adams_for_merge_no_triallelics.fam > cases
plink --bfile /data/scratch/hmy117/ms_ukb_merged \
--geno 0.1 \
--make-bed \
--make-pheno cases \* \
--allow-no-sex \
--out /data/scratch/hmy117/ms_ukb_merged_qc

# differential missingness
plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc \
--test-missing \
--allow-no-sex

# find bad SNPs
 awk '{if($5<1e-5) print $2}' plink.missing > diff_miss_snps

# QC
plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc \
--maf 0.01 \
--exclude diff_miss_snps \
--geno 0.1 \
--mind 0.1 \
--hwe 1e-5 \
--make-bed \
--allow-no-sex \
--out /data/scratch/hmy117/ms_ukb_merged_qc_round2


# kinship
./king -b /data/scratch/hmy117/ms_ukb_merged_qc_round2.bed --related --degree 3

# merge with 1kg files
awk '{print $1,$1}' rfmix.1kg_reference_superpops.tsv > kg_samples

# more QC
plink --bfile /data/home/hmy117/ADAMS/merged_genos_pca \
--keep kg_samples \
--make-bed \
--out pca_genos_kg

# rename snps
module unload plink
module load plink
plink2 --bfile pca_genos_kg \
--set-all-var-ids @:# \
--out pca_genos_kg_chrpos \
--make-bed

# merge
plink2 --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2 \
--extract pca_genos_kg_chrpos.bim \
--out /data/scratch/hmy117/ms_ukb_for_pca \
--make-bed

module unload plink
module load plink/1.9-170906
plink --bfile /data/scratch/hmy117/ms_ukb_for_pca \
--bmerge pca_genos_kg_chrpos \
--out /data/scratch/hmy117/ms_ukb_kg_merged_pca \
--make-bed \
--allow-no-sex

# qc
plink --bfile /data/scratch/hmy117/ms_ukb_kg_merged_pca \
--geno 0.01 \
--make-bed \
--out /data/scratch/hmy117/ms_ukb_kg_merged_pca_qc

# pca
./flashpca_x86-64 --bfile /data/scratch/hmy117/ms_ukb_kg_merged_pca_qc

````R
library(tidyverse)

pops = read_tsv("/data/home/hmy117/ADAMS/rfmix.1kg_reference_superpops.tsv",col_names=F)
pheno = read_table2("/data/scratch/hmy117/ms_ukb_merged_qc.fam")
pcs = read_table2("eigenvectors.txt")
pcs = pcs %>% left_join(pops %>% dplyr::rename("IID"=X1,"pop" = X2),by="IID")
pcs = pcs %>%
	mutate(cohort = ifelse(is.na(pop),"ADAMS/UKB","1000 genomes"))


# define limits of SAS clusters
sas_bounds = pcs %>%
	filter(pop=="SAS") %>%
	summarise(pc1_mu = mean(U1,na.rm=T),
	pc2_mu = mean(U2,na.rm=T),
	pc1_sd = sd(U1,na.rm=T),
	pc2_sd = sd(U2,na.rm=T)) %>%
	mutate(lower_pc1 = pc1_mu - 3*pc1_sd,
	lower_pc2 = pc2_mu - 3*pc2_sd,
	upper_pc1 = pc1_mu + 3*pc1_sd,
	upper_pc2 = pc2_mu + 3*pc2_sd)

sas_cc_samples_to_keep = pcs %>%
	filter(cohort == "ADAMS/UKB") %>%
	filter(U1 > sas_bounds$lower_pc1[1] & U1 < sas_bounds$upper_pc1[1]) %>%
	filter(U2 > sas_bounds$lower_pc2[1] & U2 < sas_bounds$upper_pc2[1])

p = ggplot(pcs %>% filter(is.na(pop) | pop == "SAS") %>%
mutate(SAS = ifelse(IID %in% sas_cc_samples_to_keep$IID,"SAS","Non-SAS")),
aes(U1,U2,col=SAS))+
	geom_point()+
	facet_wrap(~cohort)+
	theme_minimal()+labs(x="PC1",y="PC2")
write_tsv(sas_cc_samples_to_keep %>% dplyr::select(1,2),
"ukb_adams_cc_sas_samples_to_keep.tsv")
````

# filter to these samples
plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2 \
--keep ukb_adams_cc_sas_samples_to_keep.tsv \
--out /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--make-bed


# snps for pca
plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--indep-pairwise 1000 100 0.01 \
--out /data/scratch/hmy117/pca_snps \
--allow-no-sex

plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--extract /data/scratch/hmy117/pca_snps.prune.in \
--geno 0.01 \
--maf 0.05 \
--hwe 1e-10 \
--out /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only_for_pca \
--allow-no-sex \
--make-bed

# pca
./flashpca_x86-64 --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only_for_pca


#make covar
cut -f1,2,3,4,5,6,7,8,9,10,11,12 eigenvectors.txt > covars.txt

# gwas
module load plink
plink2 --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--glm hide-covar \
--covar covars.txt \
--covar-variance-standardize

plink2 --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--freq

plink --bfile /data/scratch/hmy117/ms_ukb_merged_qc_round2_sas_only \
--assoc \
--allow-no-sex



````R
library(tidyverse)
library(qqman)

df = read_table2("plink2.PHENO1.glm.logistic.hybrid",col_types="ddccccccdddddc")

# remove palindromes
df = df %>%
	filter(!(REF == "C" & ALT == "G")) %>%
	filter(!(REF == "G" & ALT == "C")) %>%
	filter(!(REF == "A" & ALT == "T")) %>%
	filter(!(REF == "T" & ALT == "A"))

# see if dodgy snps have high diff miss
missing = read_table2("plink.missing",col_types = "dcdddc")
miss_snps = missing %>% filter(P<0.05)

# filter these out
df = df %>% filter(!ID %in% miss_snps$SNP)

# freq
freqs = read_table2("plink2.afreq",col_types = "dcccdd")

df = df %>% left_join(freqs %>% dplyr::select(ID,ALT_FREQS),by="ID")

# remove MAF < 5%
df = df %>% filter(ALT_FREQS > 0.1 & ALT_FREQS < 0.9)

# manhattan
manhattan(df %>%
	dplyr::rename("CHR" = `#CHROM`,"BP"=POS,"SNP"=ID) %>%
	filter(P<0.1),
	col = c("blue4", "orange3"), suggestiveline = F,annotatePval = 5e-8)

# qq
qqman::qq(df$P)

# lambda
dchisq(median(df$P),df=1)/dchisq(0.5,df=1)










# freq of drb1*15
plink --bfile adams_filtered_pca --recode AD include-alt --snp rs3135388 --out hla_drb_1_15_snp_genos

````R
library(tidyverse)

genos = read_table2("hla_drb_1_15_snp_genos.raw",col_names=T)
pops = read_tsv("adams_pops.tsv")

genos$pop = pops$pop
colnames(genos)[7] = "geno"
genos = genos %>% mutate(drb_status = ifelse(geno > 0,"DRB_pos","DRB_neg"))
ggplot(genos,aes(pop,fill= factor(drb_status)))+geom_bar()
