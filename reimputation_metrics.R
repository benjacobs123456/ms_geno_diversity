library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
anc = args[1]

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
