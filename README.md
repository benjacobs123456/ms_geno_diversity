
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

fam_file = read_table("./outputs/ADAMS_geno_fid_iid.fam",col_names=F)
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
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.7.zip
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
  -Oz -o ../imputation_raw_files/sorted_chr$i\.vcf.gz &
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
wget https://imputationserver.sph.umich.edu/share/results/b0d354441e0952fa0e1bb107f6077518afbd68a41a7cfba0ac43448cf8d2080b/chr_1.zip &
wget https://imputationserver.sph.umich.edu/share/results/5bda742f3d72b9506d5f5e1456e1c79ea6c3a615a4341969bf146b751b32c83f/chr_10.zip &
wget https://imputationserver.sph.umich.edu/share/results/cdba779ba6f8d5da8a7005bd77d3c5be10b676db5a5bbaf80426e4a6b6685453/chr_11.zip &
wget https://imputationserver.sph.umich.edu/share/results/1deb259b781f37d19067699d77ae41cf33ec6a07c8f049bdf514b3ccc96135a0/chr_12.zip &
wget https://imputationserver.sph.umich.edu/share/results/153dc72b1e8e099eaf6ab67412f6edd001722d41ebdfdacac64d2a105967d9b6/chr_13.zip &
wget https://imputationserver.sph.umich.edu/share/results/5e050fbb9cb28458569958a1dc78537227765ed5ecde9a56fd44e76f965a662a/chr_14.zip &
wget https://imputationserver.sph.umich.edu/share/results/be1891f065300e9f6341b0ba3601d3e52ff8692b2c44761245ab52dade24aeec/chr_15.zip &
wget https://imputationserver.sph.umich.edu/share/results/b643ac3c26d21f60557d4de99f76e6b85e28eb32675e8ea8b934c9ed3edc75a2/chr_16.zip &
wget https://imputationserver.sph.umich.edu/share/results/97da216cb395d314e95c5c05c355127a99fc0073d7dab4c6bc6913d7a6bd709c/chr_17.zip &
wget https://imputationserver.sph.umich.edu/share/results/f21c3cd391643e22072f13f9abd9dd03914fe69388fcf09cf5a84c8a8814ca37/chr_18.zip &
wget https://imputationserver.sph.umich.edu/share/results/566b87b721c5a85480fd2dffa7cd7912f1f930b3cf75a5d9d06fa92f0e335486/chr_19.zip &
wget https://imputationserver.sph.umich.edu/share/results/5f6fa46ccdfbd059294d52a794099ef1f56e6ce17d451fb17206d70fd3ed32ff/chr_2.zip &
wget https://imputationserver.sph.umich.edu/share/results/54072888bc7ec6b462c62aa3b2e34b5de1139b2a0669b9fb13558a34938ae4da/chr_20.zip &
wget https://imputationserver.sph.umich.edu/share/results/64165cee0a1a69b42f23f654f92e1682e03dfeb7a31c06d8f7a0038ed7ecb9ce/chr_21.zip &
wget https://imputationserver.sph.umich.edu/share/results/c55d806681680138ce5005e66e4e5b23abdeaa8f610ff73259524f2145e67525/chr_22.zip &
wget https://imputationserver.sph.umich.edu/share/results/40f1b15bc547d2544028c6f0b5d92f18f3385c5044bc0ecc434d10015f158741/chr_3.zip &
wget https://imputationserver.sph.umich.edu/share/results/c33da06e02d31938b30f6fb7a7f707bb3c1df399b588a98d3788c9f2372b256e/chr_4.zip &
wget https://imputationserver.sph.umich.edu/share/results/046470c79d5f43df67e8157c9c46dbb992f5940517f68e6039c561e2078ab32d/chr_5.zip &
wget https://imputationserver.sph.umich.edu/share/results/ebb8436cb1053452331e9ce7548ddaeae5507d67432b6206ea483b6970c3b4f2/chr_6.zip &
wget https://imputationserver.sph.umich.edu/share/results/b5c824944158c4d7ce1f578ec0856f6401e75c1517f9edd27aa545c8a608afb0/chr_7.zip &
wget https://imputationserver.sph.umich.edu/share/results/3b39a5e809c4523399bc994d5c07c9a00d9ae2242d12700882060c48cf0a3216/chr_8.zip &
wget https://imputationserver.sph.umich.edu/share/results/073a612251df28f6d43ae8b9fe2cf8474fc9e7ca469ad2a71199c50e5704b32a/chr_9.zip &
wget https://imputationserver.sph.umich.edu/share/results/5200331737fcda9315637611e5433003245aa6f63d97812353588a956898e2f8/qc_report.txt
wget https://imputationserver.sph.umich.edu/share/results/d540798db143488c80a02009150fe01eebafc12768bffe597f2d13a22530e9d5/quality-control.html
wget https://imputationserver.sph.umich.edu/share/results/759b452647dd867bd75898cb868a2fa2370c61596c9f3afe30138f315f56d957/statistics/chunks-excluded.txt
wget https://imputationserver.sph.umich.edu/share/results/6e6614b31d77df44a663613c9348528a3c78af38591e7ff6b7afcf60af158770/statistics/snps-excluded.txt

# unzip
for i in {1..22};
do
  unzip -P joGUiWvc=?C57B -o chr_$i\.zip &
done

````



## Ancestry inference


### Download reference data
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/download_hgdp_1kg.sh
````

### Filter to ADAMS variants
````R
library(tidyverse)
df = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/adams_hg19_cpra.bim",col_names=F) %>%
  dplyr::select(2)
message(nrow(df)," SNPs")
df$X2 = str_remove_all(df$X2,"GSA-")
df = df %>%
  filter(grepl("^rs",X2))
message(nrow(df)," SNPs with rsids")

write_tsv(df,"/data/scratch/hmy117/hgdp_1kg_genomes/adams_snps_for_filtering.tsv")
````



````
### Filter HGDP genomes to ADAMS genotyped variants
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/filter_hgdp_to_adams_vars.sh
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

### Download metadata
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv



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
ggplot(adams_snps,aes(MAF,R2))+
  geom_point(alp)
adams_snps$MAF_bin = Hmisc::cut2(adams_snps$MAF,cuts = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5))
levels(adams_snps$MAF_bin) = c("<0.01","0.01-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5")
p = ggplot(adams_snps,aes(MAF_bin,R2))+
  geom_boxplot()+
  theme_bw()+
  labs(x="MAF bin",y=bquote(Imputation~R^2))
p2 = ggplot(adams_snps %>% filter(ER2 != ".") ,aes(MAF_bin,as.numeric(ER2)))+
  geom_boxplot()+
  theme_bw()+
  labs(x="MAF bin",y=bquote(Empirical~R^2))

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
dups = read_table("king.con")
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

### Remove duplicates
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

#### Explore batch effects 

````R

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

# pca plot
p = ggplot(data = adams_ukb %>% filter(!is.na(cohort)),aes(PC1_AVG,PC2_AVG,fill=predicted_ancestry))+
  geom_point(shape=21)+
	theme_bw()+
  facet_wrap(~cohort)+
  labs(x="PC1",y="PC2",fill="Ancestry")+
	scale_fill_brewer(palette="Set2")
png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/all_pcs.png",res=900,units="in",width=6,height=4)
p
dev.off()

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
  geom_vline(xintercept=0.8,linetype="dashed")


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
  geom_vline(xintercept=0.8,linetype="dashed")+
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



p = ggplot(adams_with_ethnicity,
  aes(max_prob,fill=predicted_ancestry))+
  geom_density(alpha=0.7)+
  facet_wrap(~ethnicity_clean,ncol=4)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Probability of most likely ancestry call",fill="Inferred ancestry",y="Density")+
  geom_vline(xintercept=0.8,linetype="dashed")

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
  mutate(ambiguous = ifelse(max_prob < 0.80,"yes","no"))

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
## RUNNIG FOR EUR 25-09
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/split_ancestry_groups.sh 
````

#### Crude Fisher tests prior to QC 
````unix 
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

for anc in sas afr;
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




````

#### GWAS prep 
###### Filter out PCA outliers within each ancestry

````unix 
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_risk_gwas_all_within_ancestry_prep.sh
````

#### Re-imputation
````unix 
# compress 1kg reference
qsub ~/ADAMS/genotypes/QMUL_Aug_23/scripts/prepare_kg_ref.sh

# phase with eagle
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/phase_afr.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/phase_sas.sh

# impute to 1kg with minimac4
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/impute_afr.sh
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/impute_sas.sh

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
for anc in sas afr;
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
````R 
library(tidyverse)

for(anc in c("sas","afr")){
snps = list()
for(i in c(1:22)){
    message(i)
    in_file = paste0("/data/scratch/hmy117/snp_qc_chr",i,"_anc_",anc)
    dat = read_table(in_file,col_names=F) 
    colnames(dat) = c("CHR","POS","REF","ALT","MAF","R2")
    snps[[i]] = dat
}
snps = do.call("bind_rows",snps)

hardy = read_table(paste0("/data/scratch/hmy117/hwe_reimputation_",anc,".hardy"))
missing = read_table(paste0("/data/scratch/hmy117/missing_reimputation_",anc,".vmiss"))

# combine 
snps = snps %>%    
    filter(MAF >= 0.01 & R2 >= 0.3) %>%
    mutate(ID = paste0(CHR,":",POS,":",REF,":",ALT)) %>% 
    inner_join(hardy,by="ID") %>% 
    inner_join(missing,by="ID")

# get high qual SNPs for step1 regenie 
snps = snps %>% 
    mutate(`maf_0.05` = ifelse(MAF >= 0.05,"yes","no")) %>% 
    mutate(`hwe_1e-50` = ifelse(P > 1e-50,"yes","no")) %>% 
    mutate(`info_0.9` = ifelse(R2 > 0.9,"yes","no")) %>% 
    mutate(`miss_0.1` = ifelse(F_MISS < 0.1,"yes","no"))

# get high qual snps  
step1_snps = snps %>% filter(MAF >= 0.05 & P > 1e-10 & F_MISS < 0.01 & R2 > 0.97)
write_tsv(step1_snps,paste0("/data/scratch/hmy117/info_stats_for_step1_",anc,".tsv"))
write_tsv(snps,paste0("/data/scratch/hmy117/info_stats_all_snps_",anc,".tsv"))
}


````

#### QC on imputed data
````unix 
for anc in sas afr;
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

for anc in sas afr;
  do
     # fisher's exact test
    ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_qc_$anc\_ALL_CHRS \
    --pheno ./pheno/reimputed_$anc\_pheno.tsv \
    --assoc \
    --maf 0.05 \
    --out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_fisher$anc \
    --allow-no-sex \
    --threads $NSLOTS
  done
Rscript ./scripts/fisher_to_manhattan.R afr &
Rscript ./scripts/fisher_to_manhattan.R sas &

````

#### GWAS of re-imputed data 
````unix 
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/regenie_gwas.sh"
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/plink_gwas.sh"

````

#### GWAS of matched case-control dataset 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

for(anc in c("sas","afr")){
# fx for eulidean distances
euclidean_dist = function(x,y){
	sqrt(rowSums((x - y)^2 ))
}

# read in data
pcs = read_tsv(paste0("./pheno/reimputed_",anc,"_covars_with_pcs.tsv"))
pheno = read_tsv(paste0("./pheno/reimputed_",anc,"_pheno.tsv"))
pcs = pcs %>% left_join(pheno,by=c("FID","IID"))


# define cases & controls
cases = pcs %>% filter(MS_status=="2")
controls = pcs %>% filter(MS_status=="1")

# loop through each case and find x nearest neighbours
matched_controls = list()
n_controls = 1
for(i in c(1:nrow(cases))){
this_case = cases[i,]
message("matching case", i)
matched_controls_df = do.call("bind_rows",matched_controls)
message("There are now ", nrow(matched_controls_df)," matched controls")

these_matched_controls = controls %>%
  filter(!IID %in% matched_controls_df$IID)

# make numeric matrix
these_matched_controls_mat = as.matrix(these_matched_controls %>%
dplyr::select(PC1,PC2))
this_case_mat = as.matrix(
	data.frame(
        PC1 = rep(this_case$PC1,nrow(these_matched_controls)),
	    PC2 = rep(this_case$PC2,nrow(these_matched_controls))	))

# calculate euclidean distance
euclidean_distances = euclidean_dist(this_case_mat,these_matched_controls_mat)

# add back in
these_matched_controls$euclidean = euclidean_distances

# select top x matches
these_matched_controls = these_matched_controls %>%
  slice_min(euclidean,n=n_controls,with_ties=F) %>%
	mutate(matched_case = this_case$IID)


# combine with existing controls
matched_controls[[i]] = these_matched_controls

}
matched_controls_df = do.call("bind_rows",matched_controls)

# recombine
all_dat = bind_rows(cases,matched_controls_df)

# save 
write_tsv(all_dat %>% dplyr::select(FID,IID),paste0("/data/scratch/hmy117/matched_fid_iid_anc_",anc,".tsv"))
}
````

#### Matched case-control GWAS
````unix 
echo doing $ancestry

cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

for anc in sas afr;
    do

        # fisher 
        # glm
        ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$anc\_ALL_CHRS \
        --pheno ./pheno/reimputed_$anc\_pheno.tsv \
        --keep /data/scratch/hmy117/matched_fid_iid_anc_$anc\.tsv \
        --assoc \
        --allow-no-sex \
        --maf 0.05 \
        --out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_fisher_ALL_ANCESTRY_glm_MATCHED_$anc \
        --threads $NSLOTS
    done

````

#### Basic Manhattan 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

for(anc in c("sas","afr")){
fisher = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_fisher_ALL_ANCESTRY_glm_MATCHED_",anc,".assoc"))

qqman::manhattan(
    fisher %>% 
    filter(P < 0.05))

}

#### Comparison of GWAS methods 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

for(anc in c("sas","afr")){
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


#### Manhattan plots, beta-beta plots, QQ plots
````unix 
module load R/4.2.2
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/quick_manhattan_regenie.R sas &
Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/quick_manhattan_regenie.R afr &
````

#### Clump 
````unix 
cd /data/scratch/hmy117

for anc in sas afr;
  do
    awk 'NR==1{print "SNP","P"};NR>1{if($16 > 4) print $3,10^-$16}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_$anc\_MS_status.regenie > snps_for_clump_$anc

    ~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_$anc\_ALL_CHRS \
    --clump snps_for_clump_$anc \
    --clump-p1 1e-4 \
    --clump-r2 0.01 \
    --out clumps_$anc
  done

````

#### Process top hits 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

all_res = data.frame()
for(anc in c("sas","afr")){

  # annotations
nearest = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_nearest_susceptibility_",anc),skip=30, col_types = cols(.default = "c"))
colnames(nearest)[1] = "SNP"
nearest = nearest %>% distinct(SNP, NEAREST,.keep_all=T)

# regenie

dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    left_join(nearest,by="SNP")
print("N SNPs:")
print(nrow(dat))


# read in frequencies
freqs = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_freqs_susceptibility_",anc),skip=56, col_types = cols(.default = "c"))
colnames(freqs)[1] = "SNP"
freqs = freqs %>% distinct(SNP, .keep_all=T)
freqs = freqs %>% dplyr::select(SNP,contains("AF"),Allele)
dat = dat %>%  left_join(freqs,by="SNP")

gnomad_anc_code = ifelse(anc == "afr","AFR","SAS")
freq_code = paste0("gnomADg_",gnomad_anc_code,"_AF")


# compare with imsgc
imsgc = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv",col_types = "dccdddddddc")

# sig hits
imsgc_sighits = imsgc %>% filter(P < 1e-5)

# combine
dat = dat %>% mutate(SNP = paste0(CHROM,":",GENPOS))
dat = dat %>%
    mutate(sig_imsgc = ifelse(SNP %in% imsgc_sighits$SNP,"sig","not_sig"))

# add OR
dat = dat %>% 
    mutate(OR = exp(BETA)) %>% 
    mutate(lower_ci = exp(BETA - 1.96*SE)) %>%
    mutate(upper_ci = exp(BETA + 1.96*SE)) 

# add QC data
snp_qc_data = read_tsv(paste0("/data/scratch/hmy117/info_stats_all_snps_",anc,".tsv"))

# combine
dat = dat %>%
    left_join(snp_qc_data %>% dplyr::select(CHR,POS,R2) %>%
    dplyr::rename("GENPOS_hg19"=POS,"CHROM"=CHR),by=c("CHROM","GENPOS_hg19"))

# filter to clumps 
# get clumps 
clumps = read_table(paste0("/data/scratch/hmy117/clumps_",anc,".clumped"))
clumps$chrpos = paste0(clumps$CHR,":",clumps$BP)
dat$chrpos = paste0(dat$CHR,":",dat$GENPOS_hg19)

dat = dat %>% 
  filter(chrpos %in% clumps$chrpos)

# see if near IMSGC hit
near_imsgc_hits = list()
for(i in c(1:nrow(dat))){
    this_snp = dat[i,]
    this_imsgc_locus = imsgc %>% filter(CHR == this_snp$CHR & BP > this_snp$BP - 1e6 & BP < this_snp$BP + 1e6)
    hits = this_imsgc_locus %>% filter( P < 1e-5)
    near_imsgc_hits[[i]] = nrow(hits) > 0
}
dat$near_imsgc_hits = unlist(near_imsgc_hits)

sig_hits_with_imsgc = dat %>% 
    left_join(imsgc %>% 
    dplyr::rename("OR_IMSGC" = OR,"P_IMSGC"=P), 
        by="SNP") %>% 
        mutate(ancestry=anc)
all_res <<- bind_rows(all_res,sig_hits_with_imsgc)

}


all_res_selected_cols = all_res %>% 
  dplyr::select(ancestry,ID,ALLELE0,ALLELE1,A1FREQ_CASES,A1FREQ_CONTROLS,N_CASES,N_CONTROLS,OR,lower_ci,upper_ci,
  P,NEAREST,Consequence,AF,gnomADg_NFE_AF,
  gnomADg_AFR_AF,gnomADg_SAS_AF,sig_imsgc,R2,near_imsgc_hits,
  A1,P_IMSGC,OR_IMSGC)

# flip IMSGC OR 
all_res_selected_cols = all_res_selected_cols %>% 
  mutate(OR_IMSGC = ifelse(A1 == ALLELE0,1/OR,OR))


write_csv(all_res_selected_cols,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/gwas_hits_both_ancestries.csv")

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


# filter 
filtered_imsgc = full_imsgc %>% 
  dplyr::rename("CHROM"=CHR,"GENPOS" =BP,"ALLELE1"=A1,"ALLELE0"=A2) %>%
  inner_join(sas,by=c("CHROM","GENPOS")) %>% 
  inner_join(afr,by=c("CHROM","GENPOS"))

# sig hits 
sig_hits = filtered_imsgc %>% filter(P < 5e-8)
write_tsv(sig_hits %>% dplyr::select(SNP,P),"/data/scratch/hmy117/imsgc_for_clump.tsv")

````

##### DO clumping
````unix 

cd /data/scratch/hmy117
module load R/4.2.2
module load ensembl-vep


# prepare for vep
awk 'NR>1{if($5<5e-8) print $1,$7,$7,$2"/"$3,"+",$11}' /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv > susceptibility_vep_input_ms_risk
    
# annotate with vep
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/ensembl-vep/vep \
-i susceptibility_vep_input_ms_risk \
-o snp_annotations_nearest_susceptibility_imsgc \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Consequence"



for i in {1..22};
  do
    ~/plink --bfile /data/Wolfson-UKBB-Dobson/1kg_reference/raw_1kg_plink_EUR_filesets/chr$i \
    --clump "/data/scratch/hmy117/imsgc_for_clump.tsv" \
    --clump-p1 5e-8 \
    --clump-p2 1 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --out clumped_imsgc_chr$i
  done

rm imsgc_sig_rsids
for i in {1..22};
  do 
    awk 'NR>1{print $3}' clumped_imsgc_chr$i\.clumped >> imsgc_sig_rsids
  done 


# annotate with vep - freqs 
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/ensembl-vep/vep \
-i imsgc_sig_rsids \
-o snp_annotations_freqs_susceptibility_imsgc \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--af_gnomadg \
--af \
--tab \
--fields "Uploaded_variation,Existing_variation,Allele,gnomADg_AFR_AF,gnomADg_NFE_AF,gnomADg_SAS_AF"

````

#### Concordance 
````unix 
Rscript ~/ADAMS/genotypes/QMUL_Aug_23/scripts/analyse_concordance.R
````

#### PRS 
````unix 

# PRSice
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/
awk 'NR==1{print "SNP","A1","A2","P","OR"};NR>1{print "chr"$1":"$2,$4,$5,$7,$8}' /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/discovery_metav3.0.meta > snps_for_risk_prs_all_snps

qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/prs_risk.sh  

Rscript /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/prs_risk_analysis.R 
````


#### HLA-TAPAS
````unix
qsub /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_imputation.sh 
````

#### Prepare genotypes for HIBAG
````unix 
# Liftover to hg38 for HIBAG 

for ancestry_out in sas afr;
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

````

#### HLA associations - SNP2HLA
````unix 

for ancestry in sas afr;
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
--covar-name age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_$ancestry

# recode with no. of alleles
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--recode A \
--out /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_omnibus_just_hla_genotypes

done


````

##### Define haplotypes
````unix

for ancestry in sas afr;
do

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--r2 in-phase dprime \
--ld-window-kb 999999999 \
--ld-window-r2 0 \
--out /data/scratch/hmy117/hla_ld_$ancestry

~/plink2 --vcf /data/scratch/hmy117/imputed_hla_$ancestry\.bgl.phased.vcf.gz \
--make-bed \
--double-id \
--extract /data/scratch/hmy117/hla_alleles_$ancestry \
--export haps \
--out /data/scratch/hmy117/haps_$ancestry

done




````


#### HLA associations with SNP2HLA
````unix
qsub "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_association.sh"
````


#### HIBAG QC 
````unix 
Rscript "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts/hla_qc_comparison.R" &
````

##### Analyze HIBAG results 
````unix 
Rscript ./scripts/hibag_analysis.R
````


##### HLA (SNP2HLA & HIBAG)
````R 
library(tidyverse)

all_res = data.frame()

for(ancestry in c("sas","afr")){
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
dat = read_table(in_file)

# clean allele names 
snp2hla = dat %>% 
    mutate(full_allele = str_remove_all(ID,"HLA_")) %>% 
      filter(A1 == "T") %>%
    dplyr::select(full_allele,ID) %>% 
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>% 
    mutate(full_allele = paste0(field1,":",field2)) %>% 
    filter(!is.na(field2)) %>% 
    dplyr::select(full_allele,ID) 

# save these alleles 
write_tsv(snp2hla %>% dplyr::select(ID),paste0("/data/scratch/hmy117/hla_alleles_to_keep_",ancestry,".tsv"))

    # step 1
    cmd = paste0(
        "~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_",ancestry,"_for_hla_imp_plink_just_hla ",
        "--pheno /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv ",
        "--extract /data/scratch/hmy117/hla_alleles_to_keep_",ancestry,".tsv ",
        "--covar /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv ",
        "--glm hide-covar --covar-name age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 ","--out /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry)
    system(cmd)

    # read in step1 result 
    step1_res = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_",ancestry,".MS_status.glm.logistic.hybrid"))


    # find top hit 
    tophit_step1 = step1_res %>% slice_min(P,n=1)
    all_tophits = tophit_step1

    # save to condition list 
    write_tsv(all_tophits %>% dplyr::select(ID),
    paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry),
    col_names=F)
}
run_initial_analysis("sas")
run_initial_analysis("afr")


# now condition
for(ancestry in c("sas","afr")){
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
            "--glm hide-covar --covar-name age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 ","--condition-list /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_condlist_",ancestry," ",
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
    read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_logistic_conditional_hits_afr"),col_types = cols(.default="c")) %>% mutate(anc = "AFR")
)



# save 

  
write_csv(cond,
"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputed_hla_associations_conditioned.csv")


````


##### Explore haplotypes

````R 
library(tidyverse)
library(ggsankey)


make_sankey = function(ancestry,allele){
dat = read_table(paste0("/data/scratch/hmy117/hla_ld_",ancestry,".ld")) %>% 
    mutate(full_allele = str_remove_all(SNP_A,"HLA_")) %>% 
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>% 
    mutate(full_allele = paste0(field1,":",field2)) %>% 
    filter(!is.na(field2)) %>% 
    dplyr::select(full_allele,SNP_A,everything()) %>% 
    filter(full_allele == allele)

dat = dat %>% dplyr::select(full_allele,SNP_B,PHASE,R2,DP) %>% 
  filter(PHASE %in% c("TT/AA","AA/TT"))


haps = read_table(paste0("/data/scratch/hmy117/haps_",ancestry,".haps"),col_names=F)
haps_sample = read_table(paste0("/data/scratch/hmy117/haps_",ancestry,".sample")) %>% 
filter(ID_1 != "0") 

# get two field SNP ids 
haps = haps %>% 
    mutate(full_allele = str_remove_all(X2,"HLA_")) %>% 
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>% 
    mutate(full_allele = paste0(field1,":",field2)) %>% 
    filter(!is.na(field2)) %>% 
    dplyr::select(full_allele,X2,everything()) %>% 
    dplyr::select(-X1,-X3,-X4,-X5,-contains("field")) 


ids = expand.grid(haps_sample$ID_1,c("_a","_b")) 
names = c("full_allele","SNP",paste0(paste0(ids$Var1,ids$Var2)))
colnames(haps) = names

# flip (as all orientated to absence)
flip = function(x){
  ifelse(x == 0,1,0)
}
haps = haps %>% mutate_at(.vars = -c(1:2),.funs = flip)

# make long 
haps_long = haps %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(gene = full_allele) %>% 
  separate(gene,sep="\\*",into=c("locus","other"))

# spread per gene
 haps_long = haps_long %>% 
 filter(value==1) %>% 
 filter(locus %in% c("A","B","C","DRB1","DQB1","DPB1"))

# add 'next locus'
haps_long = haps_long %>% 
  mutate(next_locus = case_when(
    locus == "A" ~ "B",
    locus == "B" ~ "C",
    locus == "C" ~ "DRB1",
    locus == "DRB1" ~ "DQB1",
    locus == "DQB1" ~ "DPB1"    
  ))

# add 'next node'
dat_for_join = haps_long %>% 
  dplyr::select(name,next_locus) %>% 
  dplyr::rename("locus"=next_locus) %>% 
  left_join(haps_long,by = c("name","locus")) %>% 
  dplyr::select(name,locus,full_allele) %>% 
  dplyr::rename("next_locus" = locus,"next_allele" = full_allele)

haps_long = haps_long %>% left_join(dat_for_join,
by=c("name","next_locus"))




# get freq
ac_af = haps_long %>% 
  group_by(locus) %>% 
  dplyr::count(full_allele) %>% 
  mutate(prop = n/sum(n)) %>% 
  ungroup() %>% 
  filter(full_allele == allele)


allele_carriers = haps_long %>% filter(full_allele == allele)

plot_dat = haps_long %>% filter(name %in% allele_carriers$name) 

plot_dat = plot_dat %>% 
left_join(plot_dat %>% group_by(locus) %>% dplyr::count(full_allele) %>% mutate(prop = n/sum(n)) %>% ungroup() %>% dplyr::select(full_allele,prop),by="full_allele")

# set labels
plot_dat = plot_dat %>% arrange((prop))
plot_dat$full_allele = factor(plot_dat$full_allele,ordered=T,levels=unique(plot_dat$full_allele))

plot_dat$locus = factor(plot_dat$locus,ordered=T,levels = c("A","B","C","DRB1","DQB1","DPB1"))
p = ggplot(plot_dat, aes(x = locus, 
               next_x = next_locus, 
               node = full_allele, 
               next_node = next_allele,
               label = paste0(full_allele, " (",round(prop,3)*100,"%)"),
               group=name,
               fill = factor(locus))) +
  geom_sankey(color="black") +
  geom_sankey_text(size = 3, color = 1, fill = "white",alpha=1,hjust = -0.2) +
  theme_sankey(base_size = 16)+
  theme(legend.position="none")+
  scale_fill_brewer(palette="Set1")+
  labs(x="HLA gene")+
  expand_limits(x = c(1, 7))+
  ggtitle(paste0(toupper(ancestry),"\n",allele,"\nAC: ",ac_af$n,"\nAF (%): ",round(ac_af$prop*100,1),"%"))

allele_outname = str_replace_all(
  str_replace_all(allele,"\\*","_"),
  "\\:","_")


png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/sankey_plot_",ancestry,"_",allele_outname,".png"),res=900,units="in",width=12,height=6)
print(p)
dev.off()

print(dat)
}

make_sankey("sas","DPB1*10:01")
make_sankey("sas","B*37:01")
make_sankey("sas","A*26:01")
make_sankey("sas","DQB1*06:03")
make_sankey("sas","DRB1*15:01")
make_sankey("sas","A*23:01")
make_sankey("afr","A*66:01")



# experimental - all haps 
library(ggsankey)
library(tidyverse)
make_sankey = function(ancestry){

haps = read_table(paste0("/data/scratch/hmy117/haps_",ancestry,".haps"),col_names=F)
haps_sample = read_table(paste0("/data/scratch/hmy117/haps_",ancestry,".sample")) %>% 
filter(ID_1 != "0") 

# get two field SNP ids 
haps = haps %>% 
    mutate(full_allele = str_remove_all(X2,"HLA_")) %>% 
    separate(full_allele,sep=":",into=c("field1","field2","field3")) %>% 
    mutate(full_allele = paste0(field1,":",field2)) %>% 
    filter(!is.na(field2)) %>% 
    dplyr::select(full_allele,X2,everything()) %>% 
    dplyr::select(-X1,-X3,-X4,-X5,-contains("field")) 


ids = expand.grid(haps_sample$ID_1,c("_a","_b")) 
names = c("full_allele","SNP",paste0(paste0(ids$Var1,ids$Var2)))
colnames(haps) = names

# flip (as all orientated to absence)
flip = function(x){
  ifelse(x == 0,1,0)
}
haps = haps %>% mutate_at(.vars = -c(1:2),.funs = flip)


# make long 
haps_long = haps %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(gene = full_allele) %>% 
  separate(gene,sep="\\*",into=c("locus","other"))

# spread per gene
 haps_long = haps_long %>% 
 filter(value==1) %>% 
 filter(locus %in% c("A","B","C","DRB1","DQB1","DPB1"))

# add 'next locus'
haps_long = haps_long %>% 
  mutate(next_locus = case_when(
    locus == "A" ~ "B",
    locus == "B" ~ "C",
    locus == "C" ~ "DRB1",
    locus == "DRB1" ~ "DQB1",
    locus == "DQB1" ~ "DPB1"    
  ))

# add 'next node'
dat_for_join = haps_long %>% 
  dplyr::select(name,next_locus) %>% 
  dplyr::rename("locus"=next_locus) %>% 
  left_join(haps_long,by = c("name","locus")) %>% 
  dplyr::select(name,locus,full_allele) %>% 
  dplyr::rename("next_locus" = locus,"next_allele" = full_allele)

haps_long = haps_long %>% left_join(dat_for_join,
by=c("name","next_locus"))

# find all six-gene haplotypes
unique_alleles = haps_long %>% 
  distinct(full_allele,locus) 

# restrict to complete data 
complete_dat = haps_long %>% 
  distinct(name,locus) %>% 
  dplyr::count(name) %>% 
  filter(n == 6)
haps_long = haps_long %>% filter(name %in% complete_dat$name)

# define six gene haplotype 

samples = unique(haps_long$name)
all_haps = list()
for(i in c(1:length(samples))){
message("doing ",i," of ",length(samples))
sample = samples[i]
  hap = paste0(
    haps_long[haps_long$name==sample & haps_long$locus=="A",]$full_allele,"-",
    haps_long[haps_long$name==sample & haps_long$locus=="B",]$full_allele,"-",
    haps_long[haps_long$name==sample & haps_long$locus=="C",]$full_allele,"-",
    haps_long[haps_long$name==sample & haps_long$locus=="DRB1",]$full_allele,"-",
    haps_long[haps_long$name==sample & haps_long$locus=="DQB1",]$full_allele,"-",
    haps_long[haps_long$name==sample & haps_long$locus=="DPB1",]$full_allele
)
  all_haps[[length(all_haps)+1]] = c("name"=sample,"haplotype"=hap)
}

# add to table 
hap_dat = do.call("bind_rows",all_haps)

# get freqs 
props = hap_dat %>% 
  dplyr::count(haplotype) %>% 
  filter(!is.na(haplotype)) %>%
  arrange(desc(n)) %>% 
  mutate(prop_hap = n/sum(n),total_haps = sum(n))
write_csv(props,paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/sankey_plot_",ancestry,"_allhaps.csv"))

# add disease disease 
pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv"))

all_haps_with_pheno = hap_dat %>% 
  separate(name,sep="_",into=c("IID","Other")) %>% 
  mutate(IID = as.numeric(IID)) %>%
  left_join(pheno,by="IID")

# loop through haps 
hap_assoc_res = list()
common_haps =  props %>% filter(prop_hap> 0.001)
unique_haps = unique(common_haps$haplotype)
for(i in c(1:length(unique_haps))){
  this_hap = unique_haps[i]
  this_dat = all_haps_with_pheno %>% 
    filter(haplotype == this_hap)
  model_dat = pheno %>% 
    mutate(has_hap = ifelse(IID %in% this_dat$IID,1,0))
  tbl = table(model_dat$has_hap,model_dat$MS_status)
  fish = fisher.test(tbl)
  
  hap_assoc_res[[i]] = c("Haplotype" = this_hap,"P" = as.numeric(fish$p.value),"OR" = as.numeric(fish$estimate))
}
hap_assoc_res = do.call("bind_rows",hap_assoc_res) %>% 
mutate(P = as.numeric(P),OR = as.numeric(OR))
write_csv(hap_assoc_res,paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/haplotype_associations_",ancestry,"_allhaps.csv"))


# add back in 
haps_long = haps_long %>% 
  left_join(hap_dat %>% dplyr::select(name,haplotype),by="name")

# set labels
plot_dat = haps_long
plot_dat$full_allele = factor(plot_dat$full_allele,ordered=T,levels=unique(plot_dat$full_allele))

plot_dat$locus = factor(plot_dat$locus,ordered=T,levels = c("A","B","C","DRB1","DQB1","DPB1"))

abundant =  props %>% 
head(n=30)

plot_dat = plot_dat %>% filter(!is.na(haplotype))
plot_dat$haplotype = factor(plot_dat$haplotype )

ggplot(plot_dat %>% filter(haplotype %in% abundant$haplotype), aes(axis1 = locus, axis2 = next_locus, y = value)) +
  geom_alluvium(aes(fill = full_allele), width = 0.25, alpha = 0.7) +
  geom_stratum(width = 0.25, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("locus", "next_locus"), expand = c(0.15, 0.15)) +
  labs(title = "Allele Transitions by Locus", x = "Locus", y = "Flow Value") +
  theme_minimal() +
  theme(legend.position = "none")


p = ggplot(plot_dat %>% filter(haplotype %in% abundant$haplotype), aes(x = locus, 
               next_x = next_locus, 
               node = full_allele, 
               next_node = next_allele,
               label = full_allele,
               fill = factor(haplotype))) +
  geom_sankey(color="black") +
  geom_sankey_text(size = 3, color = 1,alpha=1,hjust = -0.2) +
  theme_sankey(base_size = 16)+
  theme(legend.position="none")+
  labs(x="HLA gene")+
  expand_limits(x = c(1, 7))+
  ggtitle(paste0(toupper(ancestry)))


png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/sankey_plot_",ancestry,"_allhaps.png"),res=900,units="in",width=12,height=6)
print(p)
dev.off()




}
make_sankey("sas")

make_sankey("afr")

````

#### Inspect specific allele pairs

````unix 


grep "HLA_DPB1\\*10:01:01:01" /data/scratch/hmy117/hla_ld_sas.ld
grep "HLA_B\\*37:01" /data/scratch/hmy117/hla_ld_sas.ld
grep "HLA_A\\*26:01" /data/scratch/hmy117/hla_ld_sas.ld
grep "HLA_DQB1\\*06:01" /data/scratch/hmy117/hla_ld_sas.ld
grep "HLA_DRB1\\*15:01" /data/scratch/hmy117/hla_ld_sas.ld

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_DQB1*05:01:01:01"




~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_DQB1*05:01:01:01"

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_DRB1*10:01:01:01"

ancestry="sas"
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_C*06:02:01:01"

ancestry="sas"
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_B*37:01:01:01"

ancestry="sas"
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_$ancestry\_for_hla_imp_plink_just_hla \
--ld "HLA_DPB1*10:01:01:01" "HLA_A*02:06:01:01"

# AFR
~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_afr_for_hla_imp_plink_just_hla \
--ld "HLA_A*66:01:01:01" "HLA_B*53:01:01"

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_afr_for_hla_imp_plink_just_hla \
--ld "HLA_A*66:01:01:01" "HLA_C*04:01:01:01"

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_afr_for_hla_imp_plink_just_hla \
--ld "HLA_A*66:01:01:01" "HLA_DRB1*15:03:01:01"

~/plink --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_afr_for_hla_imp_plink_just_hla \
--ld "HLA_A*66:01:01:01" "HLA_DQB1*03:01:01:01"

````
#### PLINK random effects meta 
````unix 
# data prep
awk 'NR==1{print "SNP","BETA","SE","P","CHR","BP","A1","A2"};NR>1{print $1":"$18,$12,$13,10^-$15,$1,$18,$4,$3}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status.regenie_hg38 > /data/scratch/hmy117/sas_gwas_meta
awk 'NR==1{print "SNP","BETA","SE","P","CHR","BP","A1","A2"};NR>1{print  $1":"$18,$12,$13,10^-$15,$1,$18,$4,$3}' /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_afr_MS_status.regenie_hg38 > /data/scratch/hmy117/afr_gwas_meta

awk 'NR==1{print "SNP","BETA","SE","P","CHR","BP","A1","A2"};NR>1{print $11,$8,sqrt($10^2),$5,$1,$7,$2,$3}'  /data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv > /data/scratch/hmy117/imsgc_gwas_meta

# run meta
~/plink \
--meta-analysis \
/data/scratch/hmy117/imsgc_gwas_meta \
/data/scratch/hmy117/afr_gwas_meta \
/data/scratch/hmy117/sas_gwas_meta + logscale

````

````R 
library(tidyverse)

imsgc = read_table("/data/scratch/hmy117/imsgc_gwas_meta",col_types = "cdddddcc")

meta = read_table("~/plink.meta",col_types = "ddcccddddddd")

# filter meta 
meta = meta %>%
filter(N==3)

# filter imsgc 
imsgc = imsgc %>% 
filter(SNP %in% meta$SNP)

# look at novelty 
meta_hits = meta %>% 
  filter(`P(R)` < 5e-8) %>% 
  left_join(imsgc %>% dplyr::rename("IMSGC_P"=P),by="SNP")

#### MR-MEGA 

````unix 
# Prep sum stats 
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/scripts

module load R/4.2.2
Rscript mr_mega_prep.R /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status.regenie_hg38

Rscript mr_mega_prep.R /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_afr_MS_status.regenie_hg38
````

````unix 
# get allele freqs for MS risk 
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/

module load R/4.2.2
module load ensembl-vep

grep rs "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/discovery_metav3.0.meta" | awk '{print $3}' > ms_rsids

  # frequencies
/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/ensembl-vep/vep -i ms_rsids \
-o "/data/scratch/imsgc_ms_risk_discovery_hg38_for_vep_annotated" \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--af_1kg \
--tab \
--fields "Uploaded_variation,EUR_AF,Allele"

# prepare for vep
awk 'NR>1{print $1,$7,$7,$3"/"$2,"+",$11}' "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv" | sort -n -k1,2 > "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38_for_vep"



````

````R 
library(tidyverse)

in_file = "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv"
dat = read_tsv(in_file)

# format
dat$z = qnorm(1-dat$P/2)
dat$BETA = log(dat$OR)
dat$SE = abs(dat$BETA / dat$z)

dat = dat %>%
    mutate(MARKERNAME = paste0(CHR,":",BP),
    EA = A1,
    NEA = A2,
    OR_95L = exp(BETA - 1.96* SE),
    OR_95U = exp(BETA + 1.96* SE),
    EAF = NA,
    CHROMOSOME = CHR,
    POSITION = BP,
    N = 50000) %>%
    dplyr::select(MARKERNAME,EA,NEA,OR,OR_95L,OR_95U,EAF,CHROMOSOME,POSITION,N)

# add EAFs 

anno = read_table("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38_for_vep_annotated",skip=27,col_types = cols(.default = "c")) %>% dplyr::rename("MARKERNAME" = `#Uploaded_variation`) %>% 
filter(EUR_AF != "-") %>% 
dplyr::rename("EAF" = EUR_AF,"EA" = Allele) %>%
distinct() 

# join 
dat = dat %>% 
  dplyr::select(-EAF) %>% 
  inner_join(anno,by=c("EA","MARKERNAME"))

outfile = paste0(in_file,"_mrmega_input.tsv")
write_tsv(dat,outfile)
````

````unix 
cd /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs

echo "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv_mrmega_input.tsv" > mr_mega_input

echo "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_afr_MS_status.regenie_hg38_mrmega_input.tsv" >> mr_mega_input

echo "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status.regenie_hg38_mrmega_input.tsv" >> mr_mega_input

# RUN MR-MEGA
~/MR-MEGA -i mr_mega_input \
-o mr_mega_risk_res \
--pc 1
