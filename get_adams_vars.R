library(tidyverse)
df = read_table2("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno.bim",col_names=F) %>%
  dplyr::select(2)
message(nrow(df)," SNPs")
df$X2 = str_remove_all(df$X2,"GSA-")
df = df %>%
  filter(grepl("^rs",X2))
message(nrow(df)," SNPs with rsids")

write_tsv(df,"/data/scratch/hmy117/hgdp_1kg_genomes/adams_snps_for_filtering.tsv")

