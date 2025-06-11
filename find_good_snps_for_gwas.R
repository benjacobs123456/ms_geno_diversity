library(tidyverse)

args = commandArgs(trailingOnly=T)

in_file = paste0("/data/scratch/hmy117/freqs_",args[1],".frq.cc")
in_file2 = paste0("/data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_ALL_independent_snps_",args[1],".prune.in")
out_file = paste0("/data/scratch/hmy117/snps_to_keep_for_gwas",args[1],".tsv")
out_file2 = paste0("/data/scratch/hmy117/snps_to_keep_for_gwas_step1",args[1],".tsv")

freqs = read_table(in_file,col_names=T)
step1_snps = read_table(in_file2,col_names=F)

message("N snps read in: ",nrow(freqs))
# calculate delta
freqs = freqs %>%
	mutate(delta_freq = abs(MAF_A - MAF_U))

# filter to MAF 0.01 in both datasets
freqs = freqs %>%
	filter(MAF_A > 0.05 & MAF_U > 0.05)
message("N snps MAF > 0.05 in both datasets: ",nrow(freqs))

# filter to missingness <10% in both datasets
miss = read_table(paste0("/data/scratch/hmy117/missingness_",args[1],".missing"))

# find snps missing at >10% in either dataset of P diff miss < 1e-5
missing_snps = miss %>% filter(F_MISS_A > 0.1 | F_MISS_U > 0.1 | P < 1e-5)

# remove missing snps
freqs = freqs %>%
	filter(!SNP %in% missing_snps$SNP)

write_tsv(freqs %>% dplyr::select(SNP),out_file,col_names=F)

# step 1 snps
freqs = freqs %>%
	filter(MAF_A > 0.1 & MAF_U > 0.1)

step1_snps = step1_snps %>%
	filter(X1 %in% freqs$SNP) %>%
	dplyr::select(X1)
write_tsv(step1_snps,out_file2,col_names=F)
