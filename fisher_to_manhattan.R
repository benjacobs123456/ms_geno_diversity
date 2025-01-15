library(qqman)
library(tidyverse)

args = commandArgs(trailingOnly=T)

# read in 
in_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_fisher",args[1],".assoc")
dat = read_table(in_file)
out_file = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_fisher",args[1],".png")

# manhattan 
png(out_file,res=900,units="in",width=8,height=4)
print(
qqman::manhattan(dat %>% filter(P < 0.05),annotatePval=5e-8,annotateTop=T)
)
dev.off()
