library(tidyverse)
library(qqman)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)

# read in 
anc = args[1]


# regenie 

dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie"
)) %>% 
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID)


write_csv(sig_res,
paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_regenie_sigres.csv"
)
)
# plot 
outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.png"
)

png(outplot,units="in",res=900,width=8,height=4)
print(manhattan(dat %>% filter(P < 0.5), col = c("blue4", "orange3"), main = toupper(anc)))
dev.off()


