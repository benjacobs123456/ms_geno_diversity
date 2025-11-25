
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

~/plink --file /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/raw/adams_23_05_25 \
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

#### Missingness
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid \
--missing \
--out miss_check

````
#### Heterozygosity
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23
~/plink --bfile ./outputs/ADAMS_geno_fid_iid \
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

fam_file = read_table("./outputs/ADAMS_geno_fid_iid.fam",col_names=F,col_types = "dddddd")
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

fam_file %>% filter(!IID %in% covars$IID) %>% nrow()

as.character(test$IID)

````

### Export VCFs for imputation
````unix
cd /data/scratch/hmy117/

## Liftover to hg19

awk '{print "chr"$1,$4-1,$4,$2}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid.bim > hg38_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg38_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz \
hg19_bedfile \
unmapped

awk '{print $4,$1":"$3}' hg19_bedfile > hg19_snps
awk '{print $4,$3}' hg19_bedfile > hg19_snp_positions

# Update SNP positions and IDs
~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid \
--update-map hg19_snp_positions \
--make-bed \
--out adams_tmp_hg19

~/plink2 --bfile adams_tmp_hg19 \
--set-all-var-ids chr@:#:\$r\:\$a \
--make-bed \
--snps-only just-acgt \
--rm-dup exclude-all \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_hg19_cpra


# check for strand flips prior to imputation

# download refs
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

~/plink --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_hg19_cpra \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/raw_adams_freqs \
--freq

module load perl
perl HRC-1000G-check-bim.pl \
-b /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_hg19_cpra.bim \
-f /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/raw_adams_freqs.frq \
-r HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h \
--threshold 0.4 \
--plink ~/plink

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/
bash Run-plink.sh

### Sort with bcf tools
module load bcftools
for i in {1..22};
do
  bcftools sort /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_hg19_cpra-updated-chr$i\.vcf \
  -Oz -o ../imputation_raw_files/sorted_chr$i\.vcf.gz
done

````
### 1st pass imputation
- Imputation via MIS
- Rsq filter 0.001
- hg19
- Eagle
- HRC 1.1 (hg19)
- Imputation via minimac4

###########################
### Download imputed data
###########################

````unix
cd /data/scratch/hmy117/adams_imputed

# copy wget commands from MIS
wget https://imputationserver.sph.umich.edu/share/results/aa9c371f8d1e34b92a0a4a0fe49006e67b9267ad0fa3a8e2a1785405dfbb1d5d/chr_1.zip &
wget https://imputationserver.sph.umich.edu/share/results/2196df871a9b80cb206cec5cff12d264d57f019c578a5aaed47d2299ec6429ca/chr_10.zip &
wget https://imputationserver.sph.umich.edu/share/results/c8c8dd8e670f175940db02856222977deb883b80816ddd084f909e439c67400c/chr_11.zip &
wget https://imputationserver.sph.umich.edu/share/results/ba834fb685e56487e15956a6894c712d77d423bc961952e66e4f97eb9ed5c818/chr_12.zip &
wget https://imputationserver.sph.umich.edu/share/results/0868c4c66f44f397b6408b008bb92311f0898884d74281b73894f8a4910b26f7/chr_13.zip &
wget https://imputationserver.sph.umich.edu/share/results/5a2e63691e05e357194c482657a056415f3fa8915cc20cbac789a0b4981a9105/chr_14.zip &
wget https://imputationserver.sph.umich.edu/share/results/23ce997e443ee6cd1fa677710d44949b6246cf924456887ca3864bcd8305a7b4/chr_15.zip &
wget https://imputationserver.sph.umich.edu/share/results/0fdca94229f16420ddd41ddcbf17ce647ac86b93b00dbdef0c884eee0551d436/chr_16.zip &
wget https://imputationserver.sph.umich.edu/share/results/b2e3555a4a649551459f834584221a8326fbe2176aa71c07cebb6b66947e9272/chr_17.zip &
wget https://imputationserver.sph.umich.edu/share/results/3a35a9f085cc43dbf97834a693d98d2f087290efc1bfa9b3a72ca0faaafdcaf1/chr_18.zip &
wget https://imputationserver.sph.umich.edu/share/results/c947cf175737d311f9898512c926fb26f90ba9afd8b1b309b503bd3cc88d4dd9/chr_19.zip &
wget https://imputationserver.sph.umich.edu/share/results/dee3b256d930d0aafecc527dd4e9bcc0423406448cd78aab659a0dd720f64f3b/chr_2.zip &
wget https://imputationserver.sph.umich.edu/share/results/94f33650a057576994369902e09ab1ed92376d834c101c1b4d916d8256457fed/chr_20.zip &
wget https://imputationserver.sph.umich.edu/share/results/6ae964267300a88cd969972b8c68252c908d444bb17c3ef0d182ab4f150a3449/chr_21.zip &
wget https://imputationserver.sph.umich.edu/share/results/14190d8107f8eaad3b54c2c32be52b13ef1cc88dbafd0a4c982a7251bab36809/chr_22.zip &
wget https://imputationserver.sph.umich.edu/share/results/dc6d1041b1e98249a4dcc993a7305e547772a433ca0c288c261ab9540c5ec3d5/chr_3.zip &
wget https://imputationserver.sph.umich.edu/share/results/f0d322a729fcf6bda5207351e303d9c12c9d1b90d69fc88483608b4acc4ab4ef/chr_4.zip &
wget https://imputationserver.sph.umich.edu/share/results/b3f28ea53c7cc0f7bde54e34f143f9f05ceddf88abf0c658948c374c9d229608/chr_5.zip &
wget https://imputationserver.sph.umich.edu/share/results/3bf351ed5681d44090ef44dea0a988fcd23e1668ffc99330a9559d03a0cee299/chr_6.zip &
wget https://imputationserver.sph.umich.edu/share/results/02475ec5e87c0f971efae0a0b6de8a88cb43cb432c534c9ebbbf9abab376c090/chr_7.zip &
wget https://imputationserver.sph.umich.edu/share/results/bec8e13664be095fe75467480c64ba9ebfde2753a3f78e9f6e3ac24731df0026/chr_8.zip &
wget https://imputationserver.sph.umich.edu/share/results/0bba7ec97060887e4509d5d6b51017197d8f8d2072e32883f3ac999679b6411c/chr_9.zip &
wget https://imputationserver.sph.umich.edu/share/results/da4444f15dd4c5129d07c19763852139aef355b2c03c84ff92bb3e88f9521e8c/qc_report.txt &
wget https://imputationserver.sph.umich.edu/share/results/0213e66134064b3e91924eb87e2a8508a875fc810030f5250d1f4536db6d1210/quality-control.html &
wget https://imputationserver.sph.umich.edu/share/results/26fd137b2ad3a2a40ca6b520131b4626fdc7014682d21d04a782456ee454f93e/statistics/chunks-excluded.txt &

# unzip
for i in {1..22};
do
  unzip -P {a8FbUSTk2xOf! -o chr_$i\.zip &
done

````



## Ancestry inference


### Download reference data
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/download_hgdp_1kg.sh
````

### Up to here 26-05

### Download ancestry calls
````unix
~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/data_intersection/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv \
/data/scratch/hmy117/hgdp_1kg_genomes/

### Download metadata
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv
````

## QC of imputed data

### Explore imputation quality
- Check imputation quality
- Examine snps in R
- Plot info in MAF bins

### Download UKB SNP list

````unix
cd /data/scratch/hmy117/adams_imputed/
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_mfi.tgz
tar -xvf ukb_imp_mfi.tgz
````

## Assess quality of imputation
````unix
# get stats
cd /data/scratch/hmy117/adams_imputed
rm adams_all_snps_info
for i in {1..22};
  do
    bcftools query -f "%CHROM %POS %REF %ALT %MAF %R2 %ER2"  "/data/scratch/hmy117/adams_imputed/chr$i.info.gz" >> adams_all_snps_info
  done

# repeat for UKB
rm ukb_all_snps_info
for i in {1..22};
  do
    awk '{if($8>0.9 && $6>=0.01) print}' /data/scratch/hmy117/adams_imputed/ukb_mfi_chr$i\_v3.txt >> ukb_all_snps_info
  done

````

#### Inpsect in R & get SNPs with high INFO & MAF
````R
library(tidyverse)

# read ukb snps
ukb_snps = read_table("/data/scratch/hmy117/adams_imputed/ukb_all_snps_info",col_names=F)
colnames(ukb_snps) = c("cpra","rsid","pos","ref","alt","maf","alt_allele","info")

# read adams snps
# examine imputation quality for adams snps
adams_snps = read_table("/data/scratch/hmy117/adams_imputed/adams_all_snps_info",col_names=F)
colnames(adams_snps) = c("CHR","POS","REF","ALT","MAF","R2","ER2")

# plot
adams_snps$MAF_bin = Hmisc::cut2(adams_snps$MAF,cuts = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5))
levels(adams_snps$MAF_bin) = c("<0.01","0.01-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5")
p = ggplot(adams_snps,aes(MAF_bin,R2,fill=MAF_bin))+
  geom_violin()+
  geom_boxplot(width=0.1,alpha=0.5)+
  theme_bw()+
  labs(x="MAF bin",y=bquote(Imputation~R^2))+
  theme(legend.position="none")
p2 = ggplot(adams_snps %>% filter(ER2 != ".") ,aes(MAF_bin,as.numeric(ER2),fill=MAF_bin))+
  geom_violin()+
  geom_boxplot(width=0.1,alpha=0.5)+
  theme_bw()+
  labs(x="MAF bin",y=bquote(Empirical~R^2))+
  theme(legend.position="none")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_first_pass_stats.png",res=900,units="in",width=8,height=6)
cowplot::plot_grid(p,p2,ncol=1,align="v")
dev.off()

#  adams filtering (INFO > 0.9, MAF > 0.01)
adams_snps = adams_snps %>% filter(MAF > 0.01)
adams_snps %>% nrow()
adams_snps = adams_snps %>% filter(R2 > 0.90)
adams_snps %>% nrow()

# get chrpos of adams snps
adams_snps = adams_snps %>% mutate(chrpos = paste0(CHR,":",POS))
ukb_chrpos = str_split_fixed(ukb_snps$cpra,"_",2)
ukb_chrpos = str_split_fixed(ukb_chrpos[,1],":",2)
ukb_chrpos = data.frame(ukb_chrpos)

# add chrpos back in
ukb_snps$chrpos = paste0(ukb_chrpos$X1,":",ukb_chrpos$X2)

chrpos_intersection = ukb_snps %>% inner_join(adams_snps,by=c("chrpos"))

# check ref/alt match
chrpos_intersection = chrpos_intersection %>%
	filter(ref == REF & alt == ALT)

colnames(chrpos_intersection)[c(6:8)] = paste0("UKB_",colnames(chrpos_intersection)[c(6:8)])
colnames(chrpos_intersection)[c(14:16)] = paste0("ADAMS_",colnames(chrpos_intersection)[c(14:16)])

# save full snp lists of intersecting chrpos snps
feather::write_feather(chrpos_intersection,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_adams_snps_maf_info.arrow")
# chrpos_intersection = feather::read_feather("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_adams_snps_maf_info.arrow")

# save ukb snps to keep
write_tsv(chrpos_intersection %>% dplyr::select(rsid),"/data/scratch/hmy117/ukb_snps_to_keep.tsv",col_names=F)
write_tsv(chrpos_intersection %>% dplyr::select(chrpos),"/data/scratch/hmy117/adams_snps_to_keep.tsv",col_names=F)
write_tsv(chrpos_intersection %>%
  mutate(cpra = paste0(CHR,":",POS,":",REF,":",ALT)) %>%
  distinct(rsid,.keep_all=T) %>%
  dplyr::select(rsid,cpra) ,"/data/scratch/hmy117/ukb_snps_to_update.tsv",col_names=F)


# maf maf plot
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/maf_maf_plot_ukb_adams.png",res=900,units="in",width=6,height=6)
 ggplot(chrpos_intersection,aes(UKB_maf,ADAMS_MAF,fill=0.5*(UKB_info+ADAMS_R2))) +
  geom_point(alpha=0.8,shape=21)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  labs(x="UKB MAF",y="ADAMS MAF",fill="Mean INFO")
dev.off()
````


### SNP QC
- SNP QC on individual VCFs
- Conversion back to plink, hard call with threshold 0.1
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_snp_qc_imputed.sh

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
````unix
cd /data/scratch/hmy117/adams_imputed/

~/plink --bfile combined_adams_imputed \
--snps-only just-acgt \
--make-bed \
--mac 1 \
--geno 0.1 \
--out combined_adams_imputed_qc \
--chr 1-22
````

### Individual QC

#### Heterozygosity
````unix
cd /data/scratch/hmy117/adams_imputed/
~/plink --bfile combined_adams_imputed_qc \
--indep-pairwise 1000 500 0.1 \
--maf 0.05 \
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

#### Decide which duplicates to keep and which to discard
````R
library(tidyverse)
dups = read_table("king.con",col_types = cols(.default = "d"))
cov = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv")

dups = dups %>%
  dplyr::rename("IID" = ID1) %>% left_join(cov,by="IID") %>%
  dplyr::rename("ID1" = IID) %>%
  dplyr::rename("IID" = ID2) %>% left_join(cov,by="IID") %>%  
  dplyr::rename("ID2" = IID) %>%
  mutate(person_to_discard = ifelse(is.na(ageatedss.x),ID1,ID2)) %>%
  mutate(concordant_gender = ifelse(Sex.x == Sex.y,"Concordant","Discordant"))

# plot concordance
ggplot(dups,aes(ageatedss.x,ageatedss.y,color=concordant_gender))+geom_point()

# save to file
dups = dups %>%
  mutate(FID = person_to_discard,IID = person_to_discard) %>%
  dplyr::select(FID,IID)
write_tsv(dups,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/dups_to_discard.tsv")

````

### Remove duplicates & copy back to home
````unix

~/plink --bfile combined_adams_imputed_qc \
--remove /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/dups_to_discard.tsv \
--out ~/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed \
--make-bed

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

## Susceptibility GWAS


### Merge with UKB

#### Filter to compatible SNPs
##### UKB
````unix

# convert from bgen to plink
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_ukb_files.sh"

````

#### Merge ADAMS with UKB files across chromosomes
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_ukb_files.sh
````

### Filter HGDP genomes to ADAMS-UKB SNP list
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_hgdp_to_adams_vars.sh
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

# prune
for i in {1..22};
  do
    ~/plink --bfile hgdp_1kg_nodups_chr$i \
    --make-bed \
    --out hgdp_1kg_nodups_pruned_chr$i \
    --indep-pairwise 1000 100 0.1

    ~/plink --bfile hgdp_1kg_nodups_chr$i \
    --extract hgdp_1kg_nodups_pruned_chr$i\.prune.in \
    --make-bed \
    --out hgdp_1kg_nodups_pruned_chr$i
done


# try merge again
rm merge_filelist
for i in {2..22};
  do
    echo hgdp_1kg_nodups_pruned_chr$i >> merge_filelist
  done
~/plink --bfile hgdp_1kg_nodups_pruned_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg_filtered

~/king -b combined_hgdp_1kg_filtered.bed --duplicate
cut -f1,2 king.con > samples_to_remove
~/plink --bfile combined_hgdp_1kg_filtered \
--remove samples_to_remove \
--make-bed \
--out combined_hgdp_1kg_unrelated


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

# get list of all ADAMS-UKB vars
rm adams_ukb_snps_to_keep
for i in {1..22};
do
	awk '{print $2}' /data/scratch/hmy117//ukb_adams_merged_genotypes_chr$i\.bim >> adams_ukb_snps_to_keep
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
--pca approx allele-wts 50 \
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

# Define the training control & tunegrid
trainControl = trainControl(method = "cv", number = 10)
tuneGrid = expand.grid(.mtry = c(2, 3, 4, 5))

rf_fit = train(superpop ~ .,
                      data=kg_hgdp,
                      method='rf',
                      metric='Accuracy',
                      trControl = trainControl,
                      tuneGrid = tuneGrid)

# Plot feature importance
importance <- varImp(rf_fit, scale = FALSE)
plot(importance)

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

adams_ukb = adams_ukb %>%
  mutate(max_prob = predicted_ancestry_confidence$max_prob)
cov = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/adams_covars.tsv") %>% dplyr::select(-FID)
adams_fam = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed.fam",col_names=F)

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

# pca plot
adams_ukb = adams_ukb %>% filter(!is.na(cohort))

p = ggplot(data = adams_ukb,aes(PC1_AVG,PC2_AVG,fill=predicted_ancestry))+
  geom_point(shape=21)+
	theme_bw()+
  facet_wrap(~cohort)+
  labs(x="PC1",y="PC2",fill="Ancestry")+
	scale_fill_brewer(palette="Set2")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/all_pcs.png",res=900,units="in",width=6,height=4)
p
dev.off()

# numbers
adams_ukb %>%
    mutate(predicted_ancestry = ifelse(predicted_ancestry %in% c("EUR","CSA","AFR"),as.character(predicted_ancestry),"other")) %>%
    group_by(cohort,MS_status) %>%
    dplyr::count(predicted_ancestry) %>%
    mutate(total = sum(n), prop = n/total)

# plot predicted probabilities vs ethnicity
# validate with self-reported ethnicity
# join with UKB pheno
adams_ukb = adams_ukb %>%
	mutate("EID" = as.character(IID)) %>%
	left_join(ukb_pheno %>% mutate(EID = as.character(EID)),by="EID")



adams_ukb_with_ethnicity = adams_ukb %>%
	filter(!is.na(ethnic_background_f21000_0_0))

p = ggplot(adams_ukb_with_ethnicity %>%
  filter(genetic_ethnic_grouping_f22006_0_0 == "Caucasian"),
  aes(PC1_AVG,PC2_AVG,fill=max_prob))+
  geom_point(shape=21)+
  theme_bw()+
  facet_wrap(~predicted_ancestry)+
  labs(x="PC1",y="PC2",fill="Prob")
p2 = ggplot(adams_ukb_with_ethnicity %>%
  filter(genetic_ethnic_grouping_f22006_0_0 == "Caucasian"),
  aes(max_prob,fill=predicted_ancestry))+
  geom_histogram(color="black")+
  theme_bw()+
  theme(legend.position="none")+
  labs(y="N",x="Probability of ancestry call")+
  geom_vline(xintercept=0.5,linetype="dashed")


png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_caucasians.png",res=900,units="in",width=10,height=6)
cowplot::plot_grid(p,p2,ncol=2,align="h")
dev.off()


ethnicities_to_plot = adams_ukb_with_ethnicity %>%
	dplyr::count(ethnic_background_f21000_0_0) %>%
	filter(n>200)

p = ggplot(adams_ukb_with_ethnicity %>%
  filter(ethnic_background_f21000_0_0 %in% ethnicities_to_plot$ethnic_background_f21000_0_0),
  aes(max_prob,fill=predicted_ancestry))+
  geom_density(alpha=0.7)+
  facet_wrap(~ethnic_background_f21000_0_0,ncol=6)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Probability of most likely ancestry call",fill="Inferred ancestry",y="Density")+
  geom_vline(xintercept=0.5,linetype="dashed")+
  theme(legend.position="top",strip.text.x=element_text(size=10))+
  scale_x_continuous(labels = c(0,0.5,1),breaks = c(0,0.5,1))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ukb_ethnicity_vs_prob.png",res=900,units="in",width=16,height=6)
p
dev.off()

p = ggplot(adams_ukb_with_ethnicity %>% filter(ethnic_background_f21000_0_0 %in% ethnicities_to_plot$ethnic_background_f21000_0_0 ),
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
	filter(!is.na(ethnicity_clean))

table(factor(adams_with_ethnicity$ethnicity_clean),
factor(adams_with_ethnicity$predicted_ancestry))
table(adams_with_ethnicity$predicted_ancestry)


p = ggplot(adams_with_ethnicity,
  aes(max_prob,fill=predicted_ancestry))+
  geom_density(alpha=0.7)+
  facet_wrap(~ethnicity_clean,ncol=4)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Probability of most likely ancestry call",fill="Inferred ancestry",y="Density")+
  geom_vline(xintercept=0.5,linetype="dashed")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_ethnicity_vs_prob.png",res=900,units="in",width=12,height=4)
p
dev.off()

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


# plot confidence
adams_ukb = adams_ukb %>%
  mutate(ambiguous = ifelse(max_prob < 0.50,"yes","no"))

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


# ambiguous counts
# numbers
adams_ukb %>%
    mutate(predicted_ancestry = ifelse(predicted_ancestry %in% c("EUR","CSA","AFR"),as.character(predicted_ancestry),"other")) %>%
    group_by(cohort,MS_status) %>%
    filter(ambiguous == "no") %>%
    dplyr::count(predicted_ancestry) %>%
    mutate(total = sum(n), prop = n/total)



counts_inc_exclusions = adams_ukb %>%
	group_by(cohort) %>%
	dplyr::count(ambiguous,MS_status,predicted_ancestry) %>%
	ungroup() %>%
	pivot_wider(id_cols = c(cohort,MS_status,predicted_ancestry),names_from = ambiguous, values_from=n)
write_csv(counts_inc_exclusions,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_ancestry_exclusions_tbl.csv")

adams_ukb %>%
	group_by(cohort) %>%
    dplyr::count(predicted_ancestry) %>%
    mutate(total = sum(n),pct = n/sum(n)*100)

adams_ukb %>%
	group_by(cohort) %>%
    filter(ambiguous=="no") %>%
    dplyr::count(predicted_ancestry) %>%
    mutate(total = sum(n),pct = n/sum(n)*100)

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

# write all ancestry calls
adams_ukb = adams_ukb %>%
	filter(!is.na(cohort))


# trim EUR controls to 10k
set.seed(123456)
random_eur_controls = adams_ukb %>%
filter(
    cohort == "UKB" & MS_status == "Control" & predicted_ancestry=="EUR" & ethnic_background_f21000_0_0 == "British" & max_prob >0.9
) %>%
sample_n(size = 10000)

adams_ukb = adams_ukb %>%
filter(
    !(cohort == "UKB" & MS_status == "Control" & predicted_ancestry=="EUR")
) %>%
bind_rows(random_eur_controls)


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
pcs = adams_ukb %>%
  dplyr::select(-IID) %>%
	dplyr::rename("IID" = EID) %>%
	dplyr::select(IID,contains("PC"))
colnames(pcs) = c("IID",paste0("PC",c(1:50)))
combo_dat = combo_dat %>%
    left_join(pcs,by="IID") %>%
    dplyr::select(c(1:9))

write_tsv(combo_dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv")

# demographics
combo_dat = combo_dat %>%
	left_join(anc_calls,by="IID")


````


#### Identify ancestry groupings (part 2)

##### Split into ancestry groups
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/split_ancestry_groups.sh
````

#### Crude Fisher tests prior to QC
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

for anc in sas afr eur;
  do

      awk '{if($3=="MS") print $1,$1}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv >   ms_cases_$anc

     # fisher's exact test
    ~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_$anc \
    --make-pheno ms_cases_$anc \* \
    --assoc \
    --maf 0.05 \
    --out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_fisher$anc \
    --allow-no-sex \
    --threads $NSLOTS
  done
Rscript ./scripts/fisher_to_manhattan.R afr
Rscript ./scripts/fisher_to_manhattan.R sas
Rscript ./scripts/fisher_to_manhattan.R eur




````

#### GWAS prep
###### Filter out PCA outliers within each ancestry

````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_prep.sh
````

###### Get counts of excluded people
````R
library(tidyverse)
dat = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv")

combo_dat = list()

for(anc in c("afr","sas","eur")){

non_outliers = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_",anc,"_covars_with_pcs.tsv"))

dat %>% filter(IID %in% non_outliers$IID) %>%
dplyr::count(cohort,MS_status) %>%
mutate(total = sum(n),prop = n/total) %>%
print()

combo_dat[[length(combo_dat)+1]] = dat %>% filter(IID %in% non_outliers$IID) %>%
mutate(ancestry = anc)

}

combo_dat = do.call("bind_rows",combo_dat)


# make basic demographics table
library(compareGroups)
combo_dat = combo_dat %>%
	mutate(gender = ifelse(
		sex == 1,"Male","Female"
		))
tbl = compareGroups(data = combo_dat %>% mutate(full_group = paste0(cohort,"_",MS_status)),
	full_group ~ age + gender + ancestry,
	method = c(2,3,3)
	) %>%
	createTable()

export2csv(tbl,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/risk_gwas_demographics.csv")

````



#### Re-imputation
````unix
# compress 1kg reference
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/prepare_kg_ref.sh

# phase with eagle
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/phase_afr.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/phase_sas.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/phase_eur.sh
````


#### Impute to 1kg with minimac4
````unix

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/impute_afr.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/impute_sas.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/impute_eur.sh

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/imputation_step2_qc.sh

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/merge_kg_imputed_data.sh

````
#### Power calcs 
````unix
Rscript power_calcs_risk.R
````

#### Get some QC stats
````unix
# Imputation quality
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/get_imputation_metrics.sh"

# Missingness, HWE, MAF  
for anc in sas afr eur;
    do
        ~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
        --freq \
        --out /data/scratch/hmy117/freqs_reimputation_$anc

        ~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
                --missing \
                --out /data/scratch/hmy117/missing_reimputation_$anc

        ~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
        --hardy \
        --out /data/scratch/hmy117/hwe_reimputation_$anc
    done

````

#### QC data
````unix
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/reimputation_metrics.sh"
````

#### QC on imputed data
````unix
for anc in sas afr eur;
  do
    ~/plink2 --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
    --maf 0.05 \
    --hwe 1e-50 \
    --mind 0.1 \
    --geno 0.1 \
    --snps-only just-acgt \
    --make-bed \
    --threads 1 \
    --out /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$anc\_ALL_CHRS
  done
````


#### PCA on imputed data
````unix

qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_reimputed_prep.sh"


````

#### Fisher test pre-QC (to see if losing signal)
````unix
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

for anc in sas afr eur;
  do
     # fisher's exact test
    ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
    --pheno ./pheno/reimputed_$anc\_pheno.tsv \
    --assoc \
    --maf 0.05 \
    --out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_fisher$anc \
    --allow-no-sex \
    --threads 1
  done
module load R
Rscript ./scripts/fisher_to_manhattan.R afr &
Rscript ./scripts/fisher_to_manhattan.R sas &
Rscript ./scripts/fisher_to_manhattan.R eur &

````

#### GWAS of re-imputed data
````unix
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_gwas.sh"

````


#### VEP
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/run_vep.sh

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

# add Z scores
imsgc_risk = imsgc_risk %>%
  mutate(
    Z = (qnorm(1 - (P / 2)  )),
    SE = abs(BETA / Z)
    )

# replace hg19 positions
imsgc_risk = imsgc_risk %>%
  dplyr::select(-BP) %>%
  dplyr::rename("BP" = BP_hg38) %>%
  mutate(rsID = SNP) %>%
  mutate(SNP = paste0(CHR,":",BP))

# write to file
write_tsv(imsgc_risk,"imsgc_ms_risk_discovery_hg38.tsv")

# just SNPs with rsIDs
write_tsv(imsgc_risk %>% filter(grepl("rs",rsID) & !is.na(P)),"imsgc_rsids_ms_risk_discovery_hg38.tsv")

````


#### Manhattan plots, beta-beta plots, QQ plots
````unix
module load R
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/make_manhattans.R sas &
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/make_manhattans.R afr &
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/make_manhattans.R eur &

````


#### Combine all GWAS results
````R

  library(tidyverse)
  setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
  eur = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_eur_MS_status_sig_hits_with_freqs.csv",col_types = cols(.default="c")) %>% mutate(ancestry="EUR")

  sas = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status_sig_hits_with_freqs.csv",col_types = cols(.default="c")) %>% mutate(ancestry="SAS")

  afr = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_afr_MS_status_sig_hits_with_freqs.csv",col_types = cols(.default="c")) %>% mutate(ancestry="AFR")

  dat = bind_rows(eur,sas,afr) %>%
  mutate(significant = ifelse(as.numeric(P_PLINK)<5e-8,"yes","no")) %>%
  mutate(significant_regenie = ifelse(as.numeric(P_REGENIE)<5e-8,"yes","no")) %>%
  mutate(near_imsgc_hit = ifelse(as.numeric(imsgc_hit_snp_p)<1e-5,"yes","no"))

  write_csv(dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/all_gwas_clumped_res_suggestive.csv")

````


#### Comparison of GWAS methods
````R
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

for(anc in c("sas","afr","eur")){
all_snps = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_allsnps_step1",anc,"_MS_status.regenie"))
pcs_1_2 = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_pcs_1_2_",anc,"_MS_status.regenie"))
age = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_age_",anc,"_MS_status.regenie"))
main = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",anc,"_MS_status.regenie"))

# read in plink results
plink_agesex = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm_age",anc,".MS_status.glm.logistic.hybrid"))
plink_nocov = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm_nocov",anc,".MS_status.glm.logistic.hybrid"))
plink = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",anc,".MS_status.glm.logistic.hybrid"))

# combine
combo = main %>%
    mutate("Primary (REGENIE)" = LOG10P) %>%
    filter(LOG10P > -log10(0.05)) %>%
    dplyr::select(ID,"Primary (REGENIE)") %>%
    inner_join(
        all_snps %>%
        mutate("All SNPs (REGENIE)" = LOG10P) %>%
        dplyr::select(ID,"All SNPs (REGENIE)"),
        by="ID"         
    ) %>%
    inner_join(
        pcs_1_2 %>%
        mutate("PCs 1-2 (REGENIE)" = LOG10P) %>%
        dplyr::select(ID,"PCs 1-2 (REGENIE)"),
        by="ID"      
    ) %>%
    inner_join(
        age %>%
        mutate("Age (REGENIE)" = LOG10P) %>%
        dplyr::select(ID,"Age (REGENIE)"),
        by="ID"          
    ) %>%
    inner_join(
        plink %>%
        mutate("Primary (PLINK)" = -log10(P)) %>%
        dplyr::select(ID,"Primary (PLINK)"),
        by="ID"          
    ) %>%
    inner_join(
        plink_agesex %>%
        mutate("Age & Sex (PLINK)" = -log10(P)) %>%
        dplyr::select(ID,"Age & Sex (PLINK)"),
        by="ID"          
    ) %>%
    inner_join(
        plink_nocov %>%
        mutate("No Covariates (PLINK)" = -log10(P)) %>%
        dplyr::select(ID,"No Covariates (PLINK)"),
        by="ID"          
    )
# plot all PCs
png(paste0("./outputs/gwas_comparison_log10p_pairs_",anc,"_plot.png"),res=900,units="in",width=12,height=12)
print(
  GGally::ggpairs(
    combo,
    columns = c(2:ncol(combo)),
    diag = list(continuous = "blankDiag"),
    upper="blank")+
    scale_color_brewer(palette="Set1")+
    theme_bw()+
    geom_abline(slope=1,intercept=0,color="red",alpha=0.5,linetype="dashed")
)
dev.off()
}

````



#### Clump IMSGC
````R
# find SNP intersection
library(tidyverse)
full_imsgc = read_table("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/discovery_metav3.0.meta")
sas = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status.regenie"
))
afr = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_afr_MS_status.regenie"
))
eur = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_eur_MS_status.regenie"
))


# filter
filtered_imsgc = full_imsgc %>%
  dplyr::rename("CHROM"=CHR,"GENPOS" =BP,"ALLELE1"=A1,"ALLELE0"=A2) %>%
  inner_join(sas,by=c("CHROM","GENPOS")) %>%
   inner_join(eur,by=c("CHROM","GENPOS")) %>%
  inner_join(afr,by=c("CHROM","GENPOS"))

# sig hits
sig_hits = filtered_imsgc %>% filter(P < 1e-5)
write_tsv(sig_hits %>% dplyr::select(SNP,P),"/data/scratch/hmy117/imsgc_for_clump.tsv")

````

##### DO clumping
````unix

cd /data/scratch/hmy117

module load miniforge
mamba activate vep_env



for i in {1..22};
  do
    ~/plink --bfile /data/Wolfson-UKBB-Dobson/1kg_reference/raw_1kg_plink_EUR_filesets/chr$i \
    --clump "/data/scratch/hmy117/imsgc_for_clump.tsv" \
    --clump-p1 1e-5 \
    --clump-p2 1 \
    --clump-r2 0.001 \
    --clump-kb 1000 \
    --out clumped_imsgc_chr$i
  done

rm imsgc_sig_rsids
for i in {1..22};
  do
    awk 'NR>1{print $3}' clumped_imsgc_chr$i\.clumped >> imsgc_sig_rsids
  done


# annotate with vep - freqs
vep \
-i imsgc_sig_rsids \
-o snp_annotations_freqs_susceptibility_imsgc \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--af_gnomadg \
--af \
--tab \
--fields "Uploaded_variation,Existing_variation,Allele,gnomADg_AFR_AF,gnomADg_NFE_AF,gnomADg_SAS_AF"

# annotate with vep
vep \
-i imsgc_sig_rsids \
-o snp_annotations_nearest_susceptibility_imsgc \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence"


````

#### Concordance
````unix
module purge
module load R

Rscript ~/ADAMS/genotypes/QMUL_Aug_23/scripts/analyse_concordance.R
````

#### PRS
````unix

# PRSice
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
awk 'NR==1{print "SNP","A1","A2","P","OR"};NR>1{print "chr"$1":"$2,$4,$5,$7,$8}' /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/discovery_metav3.0.meta > snps_for_risk_prs_all_snps

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/prs_risk.sh  
module load R

Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/prs_risk_analysis.R
````


#### HLA-TAPAS
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_imputation.sh
````

#### Prepare genotypes for HIBAG
````unix
# Liftover to hg38 for HIBAG

for ancestry_out in sas afr eur;
  do
    ## Liftover to hg38
    cd /data/scratch/hmy117
    awk '{print "chr"$1,$4-1,$4,$2}' /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS.bim > hg19_bedfile$ancestry_out

    # run liftover
    /data/Wolfson-UKBB-Dobson/liftover/liftOver \
    hg19_bedfile$ancestry_out \
    /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
    hg38_bedfile$ancestry_out \
    unmapped

    awk '{print $4,$1":"$3}' hg38_bedfile$ancestry_out > hg38_snps$ancestry_out
    awk '{print $4,$3}' hg38_bedfile$ancestry_out > hg38_snp_positions$ancestry_out

    # Update SNP positions and IDs
    ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
    --update-map hg38_snp_positions$ancestry_out \
    --make-bed \
    --chr 6 \
    --out /data/scratch/hmy117/risk_gwas_genotypes_hg38_tmp$ancestry_out

    # remove dups, just SNPs, restrict pheno
    ~/plink2 --bfile /data/scratch/hmy117/risk_gwas_genotypes_hg38_tmp$ancestry_out \
    --make-bed \
    --snps-only just-acgt \
    --rm-dup exclude-all \
    --chr 6 \
    --keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry_out\_pheno.tsv \
    --out /data/scratch/hmy117/risk_gwas_genotypes_hg38_hibag$ancestry_out

    # hg19 file
    ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$ancestry_out\_ALL_CHRS \
    --chr 6 \
    --keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/susceptibility_ALL_$ancestry_out\_pheno.tsv \
    --out /data/scratch/hmy117/risk_gwas_genotypes_hg19_hibag$ancestry_out \
    --make-bed
  done



````

#### HIBAG
````unix
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hibag_imputation_afr.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hibag_imputation_sas.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hibag_imputation_eur.sh"
````

#### HLA associations - SNP2HLA
````unix

for ancestry in eur sas afr;

do

cd /data/scratch/hmy117

# make plink file with phased vcf
~/plink2 --vcf /data/scratch/hmy117/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--make-bed \
--double-id \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink \
--threads $NSLOTS

# get HLA alleles
grep HLA /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink.bim | cut -f2 > /data/scratch/hmy117/hla_alleles_$ancestry

# filter to just hla alleles
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink \
--extract /data/scratch/hmy117/hla_alleles_$ancestry \
--make-bed \
--maf 0.01 \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla

# frequencies
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--freq case-control \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry\_pheno.tsv \
--allow-no-sex \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_freqs_$ancestry


# run fisher test
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry\_pheno.tsv \
--allow-no-sex \
--assoc \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_fisher_$ancestry


# run basic logistic regression
~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry\_pheno.tsv \
--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_$ancestry\_covars_with_pcs.tsv \
--glm hide-covar \
--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_$ancestry

# recode with no. of alleles
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--recode A \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla_genotypes

done


````

#### HIBAG QC
````unix
Rscript "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_qc_comparison.R"
````

##### Analyze HIBAG results
````unix
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts//hibag_analysis.R
````


##### HLA (SNP2HLA & HIBAG)
````R
library(tidyverse)

all_res = data.frame()

for(ancestry in c("sas","afr","eur")){
# read in SNP2HLA
in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry,".MS_status.glm.logistic.hybrid")
dat = read_table(in_file)

# read in HIBAG
hibag = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_additive.csv") %>%  
    filter(anc == ancestry) %>%
    dplyr::select(full_allele,beta,p,af_Control,af_MS) %>%
    dplyr::rename("hibag_beta"=beta,"hibag_p"=p)

# clean allele names
snp2hla = dat %>%
    mutate(full_allele = str_remove_all(ID,"HLA_")) %>%
      filter(A1 == "T") %>%
    dplyr::select(full_allele,OR,P) %>%
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>%
    mutate(full_allele = paste0(field1,":",field2)) %>%
    filter(!is.na(field2)) %>%
    dplyr::select(full_allele,OR,P) %>%
    mutate(snp2hla_beta = log(OR),snp2_hla_p = P)


# snp2hla freqs
snp2hla_freqs = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_freqs_",ancestry,".frq.cc")) %>%
  filter(A1 == "T") %>%
    mutate(full_allele = str_remove_all(SNP,"HLA_")) %>%
    dplyr::select(full_allele,MAF_U,MAF_A) %>%
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>%
    mutate(full_allele = paste0(field1,":",field2)) %>%
    filter(!is.na(field2)) %>%
    dplyr::select(full_allele,MAF_A,MAF_U)

# combine
combo_dat =
    hibag %>%
        left_join(snp2hla,by="full_allele") %>%
        left_join(snp2hla_freqs,by="full_allele")

# make fdr
combo_dat$hibag_fdr = p.adjust(combo_dat$hibag_p,method="fdr")
combo_dat$snp2hla_fdr = p.adjust(combo_dat$snp2_hla_p,method="fdr")

# beta beta
p = ggplot(combo_dat,aes(-log10(hibag_fdr),-log10(snp2hla_fdr)))+
    geom_point()+
    theme_bw()+
    geom_vline(xintercept=-log10(0.05),linetype="dashed",color="red",alpha=0.5)+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="red",alpha=0.5)+
    ggrepel::geom_text_repel(data = combo_dat %>% filter(hibag_fdr<0.05),
    aes(label = full_allele))+
    ggtitle(toupper(ancestry))+
    labs(x=bquote(-log[10]~FDR(HIBAG)),
    y=bquote(-log[10]~FDR(SNP2HLA)))
png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hibag_vs_snp2hla_p_vs_p",ancestry,".png"),res=900,units="in",width=4,height=4)
print(p)
dev.off()

# beta beta
p = ggplot(combo_dat,aes(hibag_beta,snp2hla_beta))+
    geom_vline(xintercept=0,linetype="dashed",color="red",alpha=0.5)+
    geom_hline(yintercept=0,linetype="dashed",color="red",alpha=0.5)+
    geom_abline(intercept=0,slope=1,linetype="dashed",color="red",alpha=0.5)+
    geom_point(data = combo_dat %>% filter(hibag_p>=0.05),alpha=0.2)+
    geom_point(data = combo_dat %>% filter(hibag_p<0.05),alpha=1)+
    theme_bw()+
    ggrepel::geom_text_repel(data = combo_dat %>% filter(hibag_p<0.05),
    aes(label = full_allele),size=3)+
    ggtitle(toupper(ancestry))+
    labs(x="Effect size (logOR) from HIBAG",y="Effect size (logOR) from SNP2HLA")
png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hibag_vs_snp2hla_",ancestry,".png"),res=900,units="in",width=4,height=4)
print(p)
dev.off()

# save data
write_csv(combo_dat,paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hibag_vs_snp2hla_",ancestry,".csv"))

all_res <<- bind_rows(all_res,combo_dat %>% mutate(anc = ancestry))
}

all_res = all_res %>%
distinct(anc,full_allele,.keep_all=T) %>%
mutate(snp2hla_or = exp(snp2hla_beta),hibag_or = exp(hibag_beta)) %>%
dplyr::select(-P,-OR,-hibag_fdr,-snp2hla_fdr,-hibag_beta,-snp2hla_beta) %>%
pivot_wider(id_cols = full_allele,names_from = anc,
values_from = c(2,3,4,5,6,7,9,10)) %>%
arrange(full_allele) %>%
dplyr::select(full_allele,contains("sas"),contains("afr"),contains("eur"))

write_csv(all_res,paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hibag_vs_snp2hla_both_ancs.csv"))

````


##### Conditional analysis HLA results
````R
library(tidyverse)


# conditional stepwise analysis
run_initial_analysis = function(ancestry){


  # find two-field alleles
# read in SNP2HLA
in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry,".MS_status.glm.logistic.hybrid")
dat = read_table(in_file,col_types = "ddccccccdddddc")

# clean allele names
snp2hla = dat %>%
    mutate(full_allele = str_remove_all(ID,"HLA_")) %>%
      filter(A1 == "T") %>%
    dplyr::select(full_allele,ID) %>%
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>%
    mutate(full_allele = paste0(field1,":",field2)) %>%
    filter(!is.na(field2)) %>%
    dplyr::select(full_allele,ID) %>%
    filter(!grepl("HLA_DPA1",ID)) %>%
    filter(!grepl("HLA_DQA1",ID))


# save these alleles
write_tsv(snp2hla %>% dplyr::select(ID),paste0("/data/scratch/hmy117/hla_alleles_to_keep_",ancestry,".tsv"))

    # step 1
    cmd = paste0(
        "~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_",ancestry,"_for_hla_imp_plink_just_hla ",
        "--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv ",
        "--extract /data/scratch/hmy117/hla_alleles_to_keep_",ancestry,".tsv ",
        "--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv ",
        "--glm hide-covar --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 ","--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry)
    system(cmd)

    # read in step1 result
    step1_res = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry,".MS_status.glm.logistic.hybrid"))


    # find top hit
    tophit_step1 = step1_res %>% slice_min(P,n=1)
    print(nrow(step1_res))

    all_tophits = tophit_step1

    # save to condition list
    write_tsv(all_tophits %>% dplyr::select(ID),
    paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry),
    col_names=F)
}
run_initial_analysis("sas")
run_initial_analysis("afr")
run_initial_analysis("eur")


# now condition
for(ancestry in c("sas","afr","eur")){
  # read initial hit
  all_tophits = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry),col_names=F)
  colnames(all_tophits)="ID"
  # iterate over further steps
  for(i in c(2:100)){

    message("Doing iteration ",i)

    cmd = paste0(
            "~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_",ancestry,"_for_hla_imp_plink_just_hla --extract /data/scratch/hmy117/hla_alleles_to_keep_",ancestry,".tsv ",
            "--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv ",
            "--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv ",
            "--glm hide-covar --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 ","--condition-list /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry," ",
        "--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_step",i,"_",ancestry)
    system(cmd)

    # read in step i result
    step_res = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_step",i,"_",ancestry,".MS_status.glm.logistic.hybrid"))

    # check top hit is below P < 0.05
    nhit = step_res %>% filter(P < 0.05) %>% nrow()
    if(nhit == 0){
        break
    } else {

        # find top hit
        tophit = step_res %>% slice_min(P,n=1)

        # add to overall top hits
        all_tophits <<- bind_rows(all_tophits,tophit)

        # save to condition list
        write_tsv(all_tophits %>% dplyr::select(ID),
        paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry),
        col_names=F)
      }
  }

# save overall conditional hits
 write_tsv(all_tophits,
    paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_conditional_hits_",ancestry)
 )

}

# read conditional hits
cond = bind_rows(
    read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_conditional_hits_sas"),col_types = cols(.default="c")) %>% mutate(anc = "SAS"),
    read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_conditional_hits_afr"),col_types = cols(.default="c")) %>% mutate(anc = "AFR"),
    read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_conditional_hits_eur"),col_types = cols(.default="c")) %>% mutate(anc = "EUR")    
)


cond = cond %>%
    mutate(full_allele = str_remove_all(ID,"HLA_")) %>%
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>%
    mutate(full_allele = paste0(field1,":",field2)) %>%
    filter(!is.na(field2)) %>%
    dplyr::select(full_allele,anc,OR,P) %>%
    group_by(anc) %>%
    mutate(condstep = row_number())

# save

write_csv(cond,
"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_conditioned.csv")


````


##### PAF and European alleles
````R
library(tidyverse)


##############################
# PAF
##############################

library(tidyverse)

# read in moutsianas alleles & join with HIBAG
eur_alleles= read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/moutsianas_alleles.csv")
hibag_res = read_csv ("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_additive.csv") %>%
dplyr::rename("hla_allele"=full_allele) %>%
dplyr::select(hla_allele,OR,contains("af_"),p,anc) %>%
pivot_wider(names_from = anc,id_cols = hla_allele,values_from = c(2:5)) %>%
filter(hla_allele %in% eur_alleles$hla_allele)
eur_alleles = eur_alleles %>% inner_join(hibag_res,by="hla_allele")
write_csv(eur_alleles,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/moutsianas_alleles_complete.csv")

hla_dat1 = eur_alleles %>%
  mutate(sas_paf = (af_MS_sas^2 + 2*af_MS_sas*(1-af_MS_sas) ) * (1 - 1 / OR_sas)  ) %>%
  mutate(afr_paf = (af_MS_afr^2 + 2*af_MS_afr*(1-af_MS_afr) ) * (1 - 1 / OR_afr)  ) %>%
  mutate(eur_paf = (af_MS_eur^2 + 2*af_MS_eur*(1-af_MS_eur) ) * (1 - 1 / OR_eur)  ) %>%
  mutate(eur_paf_imsgc = (eur_case_af^2 + 2*eur_case_af*(1-eur_case_af) ) * (1 - 1 / eur_or)  ) %>%
  pivot_longer(contains("paf")) %>%
  dplyr::select(hla_allele,name,value) %>%
  mutate(name = str_remove_all(name,"_paf")) %>%
  mutate(paf = 100* value) %>%
  mutate(PAF = "Miettinen's")


hla_dat2 = eur_alleles %>%
  mutate(sas_paf = ( (af_Control_sas^2 + 2*af_Control_sas*(1-af_Control_sas) ) * (OR_sas - 1) ) / (1 + (af_Control_sas^2 + 2*af_Control_sas*(1-af_Control_sas) ) * (OR_sas - 1) ) ) %>%
   mutate(afr_paf = ( (af_Control_afr^2 + 2*af_Control_afr*(1-af_Control_afr) ) * (OR_afr - 1) ) / (1 + (af_Control_afr^2 + 2*af_Control_afr*(1-af_Control_afr) ) * (OR_afr - 1) ) ) %>%
  mutate(eur_paf = ( (eur_cont_af^2 + 2*eur_cont_af*(1-eur_cont_af) ) * (eur_or - 1) ) / (1 + (eur_cont_af^2 + 2*eur_cont_af*(1-eur_cont_af) ) * (eur_or - 1) ) ) %>%
  pivot_longer(contains("paf")) %>%
  dplyr::select(hla_allele,name,value) %>%
  mutate(name = str_remove_all(name,"_paf")) %>%
  mutate(paf = 100* value) %>%
  mutate(PAF = "Levin's")

dat =  bind_rows(hla_dat1,hla_dat2)


p = ggplot(dat,aes(
  hla_allele,paf,fill=toupper(name)
)) +
geom_col(color="black",position=position_dodge())+
theme_bw()+
labs(x="HLA allele",y="PAF (%)",fill="Ancestry")+
coord_flip()+
scale_fill_brewer(palette="Set1")+
facet_wrap(~PAF,nrow=2)


````

#### LD

````unix

# get LD
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_sas_for_hla_imp_plink_just_hla \
--ld-snp "HLA_A*33:03:01:01" \
--ld-window-kb 1000 \
--r2 \
--ld-window 99999 \
--ld-window-r2 0.1

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_sas_for_hla_imp_plink_just_hla \
--ld-snp "HLA_A*26:01:01:01" \
--ld-window-kb 1000 \
--r2 \
--ld-window 99999 \
--ld-window-r2 0.1
````

#### Prep for GWAS catalogue
````R 
library(tidyverse)
for(anc in c("afr","sas","eur")){

# PLINK GWAS
# hg19 names
plink_dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
    anc,".MS_status.glm.logistic.hybrid_hg38"
)) 

# add REGENIE GWAS
dat_regenie = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    mutate(SNP = paste0(CHROM,":",GENPOS))


# format 
plink_dat = plink_dat %>% 
  dplyr::rename("chromosome"=`#CHROM`,
    "base_pair_location" = POS,
    "effect_allele"=A1) %>% 
    mutate("other_allele" = ifelse(effect_allele==REF,ALT,REF),
    "beta" = log(OR)) %>% 
    dplyr::rename("standard_error" = `LOG(OR)_SE`,
    "p_value" = P)
dat_regenie = dat_regenie %>% 
  dplyr::select(CHROM,GENPOS,A1FREQ,ALLELE1,ALLELE0) %>% 
  dplyr::rename("chromosome"=CHROM,"base_pair_location"=GENPOS)
plink_dat = plink_dat %>% 
  dplyr::select(chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,p_value) %>% 
  left_join(dat_regenie,by=c("chromosome","base_pair_location"))
plink_dat = plink_dat %>% 
  mutate(effect_allele_frequency = ifelse(
    effect_allele == ALLELE1,A1FREQ,1-A1FREQ
  )) %>% 
  filter(effect_allele==ALLELE1 | effect_allele==ALLELE0)

plink_dat = plink_dat %>% 
  dplyr::select(chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,effect_allele_frequency,p_value) 
write_tsv(plink_dat,paste0("/data/scratch/hmy117/sumstats_for_gwas_catalogue_plink_hg38_",anc,".tsv"))
}