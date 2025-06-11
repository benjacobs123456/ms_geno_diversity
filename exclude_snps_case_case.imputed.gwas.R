library(tidyverse)
library(qqman)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)

# read in 
anc = args[1]
cc_gwas = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/fisher_case_case_",anc,".assoc"))
cc_gwas_counts = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/fisher_case_case_counts_",anc,".assoc"))

cc_gwas = cc_gwas %>% 
    left_join(cc_gwas_counts %>% 
        dplyr::select(SNP,C_A,C_U),
        by="SNP")

# get AF diff
cc_gwas = cc_gwas %>% mutate(delta_af = abs(F_U - F_A)) 

# exclusions
cc_gwas = cc_gwas %>% 
    mutate(exclude = ifelse
        (delta_af>0.2 | P < 0.05 | (C_U == 0 & C_A == 0) | (F_A < 0.01 & F_U < 0.01),
        "exclude",
        "keep")
    )

n_excl = sum(cc_gwas$exclude=="exclude") 
pct_excl = round(100* n_excl / nrow(cc_gwas),1)
lab = paste0(n_excl," SNPs excluded (",pct_excl,"%)")


png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/fisher_case_case_delta_af_plot_imputed",anc,".png"),res=900,units="in",width=6,height=4)
ggplot(cc_gwas,aes(F_A,F_U))+
    geom_point(data = cc_gwas %>% filter(exclude == "keep"), shape=21,fill="orange",stroke=0.1)+
    geom_point(data = cc_gwas %>% filter(exclude == "exclude"), shape=21,alpha=0.01,fill="purple",stroke=0.1)+    
    theme_bw()+
    labs(x = "AF (UKB)", y = "AF (ADAMS)", fill ="Excluded")+
    geom_abline(intercept=0,slope=1,lwd=3,linetype="dashed",alpha=0.1)+
    ggtitle(lab)
dev.off()

write_tsv(
    cc_gwas %>% filter(exclude == "keep") %>% 
    dplyr::select(SNP),
    paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/snps_to_keep_case_case_gwas_imputed",anc)
    )


