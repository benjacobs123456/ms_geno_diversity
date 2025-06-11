library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed/")
snps = list()
for(i in c(1:22)){
		snps[[i]] = read_table2(paste0("chr",i,".info.gz"))
}
snps = do.call("bind_rows",snps)
nrow(snps)

# cut into MAF bins
snps$maf_bin = Hmisc::cut2(snps$MAF,cuts = c(0,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.5))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_maf.png",res = 300, units="in",width=4, height=4)
ggplot(snps,aes(maf_bin,Rsq,fill=maf_bin))+geom_boxplot()+
theme_minimal()+
labs(x="Minor allele frequency",y="Imputation quality")+
scale_x_discrete(labels = c("<0.001","<0.01","<0.05","<0.1","<0.5"))+
scale_fill_brewer(palette="Set1")+
theme(legend.position = "none")
dev.off()

# genotyped snps
genotyped = snps %>% filter(Genotyped =="Genotyped")

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/imputation_quality_r.png",res = 300, units="in",width=4, height=4)
ggplot(genotyped,aes(maf_bin,as.numeric(EmpR),fill=maf_bin))+
geom_boxplot()+
theme_minimal()+
labs(x="Minor allele frequency",y="Imputation-Genotyping correlation (R)")+
scale_x_discrete(labels = c("<0.01","<0.05","<0.1","<0.5"))+
scale_fill_brewer(palette="Set1")+
theme(legend.position = "none")
dev.off()
