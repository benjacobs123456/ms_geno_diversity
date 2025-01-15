library(tidyverse)

args = commandArgs(trailingOnly=T)
# args = c("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_sas_MS_status.regenie")
dat = read_table(args[1])

# format
dat = dat %>%
    mutate(MARKERNAME = paste0(CHROM,":",GENPOS),
    EA = ALLELE1,
    NEA = ALLELE0,
    OR = exp(BETA),
    OR_95L = exp(BETA - 1.96* SE),
    OR_95U = exp(BETA + 1.96* SE),
    EAF = A1FREQ,
    CHROMOSOME = CHROM,
    POSITION = GENPOS) %>%
    dplyr::select(MARKERNAME,EA,NEA,OR,OR_95L,OR_95U,EAF,CHROMOSOME,POSITION,N)
outfile = paste0(args[1],"_mrmega_input.tsv")
write_tsv(dat,outfile)
