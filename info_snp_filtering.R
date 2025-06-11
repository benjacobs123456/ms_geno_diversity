library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")

# read in SNPs
all_snps = list()
for(i in c(1:22)){
  file = paste0("chr",i,".info.gz")
  all_snps[[i]] = read_table(file)

}
snps = do.call("bind_rows",all_snps)


# filter adams snps
adams_snps = snps %>%
  filter(Rsq > 0.3 & MAF >= 0.01)

# write ADAMS SNPs to keep
adams_snps_to_keep = adams_snps %>%
  dplyr::select(SNP)
write_tsv(adams_snps,"/data/scratch/hmy117/adams_imputed/adams_snps_to_keep.tsv")

# write ADAMS SNPs to keep for REGENIE step 1
adams_snps_to_keep_for_step1 = adams_snps %>%
  filter(Genotyped=="Genotyped" & MAF >= 0.05) %>%
  dplyr::select(SNP)
write_tsv(adams_snps_to_keep_for_step1,"/data/scratch/hmy117/adams_imputed/adams_snps_to_keep_for_regenie_step1.tsv")
