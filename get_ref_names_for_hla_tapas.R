library(tidyverse)

args = commandArgs(trailingOnly=T)
in_file = args[1]

ref_bim = read_table("/data/home/hmy117/HLA-TAPAS/resources/1000G.bglv4.bim",col_names=F)
target_bim = read_table(in_file,col_names=F)

# merge
snp_names_file = target_bim %>%
  inner_join(ref_bim,by=c("X1","X4","X5","X6")) %>%
  dplyr::select(X2.x,X2.y)

# save
write_tsv(snp_names_file,paste0(in_file,"_snp_name_updates.tsv"),col_names=F)
message("written file successfully")
