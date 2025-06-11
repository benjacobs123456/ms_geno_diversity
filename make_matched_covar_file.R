library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
# args = c("CSA","sas")
# args = c("AFR","afr")

set.seed(123456)

# read in new pcs
pcs = read_table(paste0("/data/scratch/hmy117/matched_case_control_pcs_",args[2]))

# read current covars and pheno
covars = read_tsv(paste0("./pheno/susceptibility_ALL_",args[2],"_covars_with_pcs.tsv"))
pheno = read_tsv(paste0("./pheno/susceptibility_ALL_",args[2],"_pheno.tsv"))

# combine all of the above
all_cov = pheno %>%
    inner_join(covars,by=c("FID","IID")) %>%
    dplyr::select(-contains("PC")) %>%
    inner_join(pcs,by=c("FID","IID"))


# standardize pcs
z_score = function(x){
	( as.numeric(x) - mean( as.numeric(x)) ) / sd(as.numeric(x))
}
all_cov = all_cov %>%
	mutate_at(.vars = vars(contains("PC")),z_score)


p = ggplot(all_cov,aes(PC1,PC2,col=factor(MS_status)))+
    geom_point(alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")


# save plot
png(paste0("./outputs/susceptibility_ALL_",args[2],"_pc_plot_round2.png"),res=900,units="in",width=6,height=6)
p
dev.off()

# prepare pheno & covar file for whole ancestry
covars_all = all_cov %>%
filter(!is.na(PC1)) %>%
dplyr::select(FID,IID,age,sex,contains("PC")) %>%
na.omit()

write_tsv(covars_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_covars_with_pcs_matched.tsv"))

# make pheno file
pheno_all = all_cov %>%
filter(!is.na(PC1)) %>%
dplyr::select(FID,IID,MS_status) %>%
na.omit() %>%
write_tsv(pheno_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_pheno_matched.tsv"))
