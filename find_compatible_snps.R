library(tidyverse)

args = commandArgs(trailingOnly=T)
file1 = args[1]
file2 = args[2]

# read in files 
message("Reading in bim files")
file1_data = read_table(file1,col_names=F) %>% mutate(chrpos = paste0(X1,":",X4))
file2_data = read_table(file2,col_names=F) %>% mutate(chrpos = paste0(X1,":",X4))


# find intersection of chr:pos 
compatible_positions = file1_data %>%
  inner_join(file2_data,by="chrpos")
message("There are ",nrow(compatible_positions)," compatible chr:pos SNPs")

# filter to compatible alleles & filter out palindromes
compatible_positions_and_alleles = compatible_positions %>%
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


write_tsv(compatible_positions_and_alleles %>% dplyr::select(X2.x),"snps_to_keep_file1.tsv",col_names=F)
write_tsv(compatible_positions_and_alleles %>% dplyr::select(X2.y),"snps_to_keep_file2.tsv",col_names=F)

