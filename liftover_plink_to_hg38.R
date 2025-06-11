library(tidyverse)

args = commandArgs(trailingOnly=T)


# read in regenie file
infile = args[1]
gwas_hg19 = read_table(infile)

# prepare bed file for liftover
gwas_hg19_bed = gwas_hg19 %>%
  mutate(chr = paste0("chr",`#CHROM`),bp1=POS-1,bp2=POS,snp = ID) %>%
  dplyr::select(chr,bp1,bp2,snp)

# save to file
write_tsv(gwas_hg19_bed,"/data/scratch/hmy117/gwas_hg19_bedfile_tmp",col_names=F)

# run liftover
cmd = paste0("/data/Wolfson-UKBB-Dobson/liftover/liftOver ",
"/data/scratch/hmy117/gwas_hg19_bedfile_tmp ",
"/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz ",
"/data/scratch/hmy117/gwas_hg19_bedfile_tmp_hg38 unmapped")
system(cmd)

# read back in
gwas_hg38 = read_table("/data/scratch/hmy117/gwas_hg19_bedfile_tmp_hg38",col_names=F) %>%
dplyr::select(X4,X3) %>%
dplyr::rename("ID" = X4,"POS_hg38" = X3)

# merge
gwas_hg19_hg38 = gwas_hg19 %>%
  mutate("POS_hg19" = POS) %>%
  filter(ID %in% gwas_hg38$ID) %>%
  left_join(gwas_hg38,by="ID") %>%
  dplyr::select(-POS) %>%
  dplyr::rename("POS" = POS_hg38) %>%
  mutate(ID = paste0(`#CHROM`,":",POS,":",REF,":",ALT))

# write to file
write_tsv(gwas_hg19_hg38,paste0(infile,"_hg38"))
