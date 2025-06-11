
args = commandArgs(trailingOnly=T)
library(tidyverse)

infile = paste0( "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",args[1],".MS_status.glm.logistic.hybrid_hg38")
dat = read_table(infile)
regenie_infile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",args[1],"_MS_status.regenie_hg38")
dat2 = read_table(regenie_infile) %>% 
dplyr::select(ID,A1FREQ)
dat = dat %>%
inner_join(dat2,by="ID") %>%
distinct(ID,.keep_all=T)

# format
dat = dat %>%
    mutate(MARKERNAME = paste0(`#CHROM`,":",POS),
    EA = A1,
    NEA = ifelse(A1 == REF,ALT,REF),
    BETA = log(OR),
    N = OBS_CT,
    SE = `LOG(OR)_SE`,
    OR_95L = exp(BETA - 1.96* SE),
    OR_95U = exp(BETA + 1.96* SE),
    EAF = A1FREQ,
    CHROMOSOME = `#CHROM`,
    POSITION = POS) %>%
    dplyr::select(MARKERNAME,EA,NEA,OR,OR_95L,OR_95U,EAF,CHROMOSOME,POSITION,N)
outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/",args[1],"_mrmega_input.tsv")
write_tsv(dat,outfile)
