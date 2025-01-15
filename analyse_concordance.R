library(tidyverse)

setwd("/data/scratch/hmy117")
# read in clumped imsgc snps 
clumps = list() 
for(i in c(1:22)){
    clumps[[i]] = read_table(paste0("clumped_imsgc_chr",i,".clumped"))
}
imsgc = do.call("bind_rows",clumps)
full_imsgc = read_table("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/discovery_metav3.0.meta")
full_imsgc = full_imsgc %>% filter(SNP %in% imsgc$SNP)
full_imsgc$BETA = log(full_imsgc$OR)

nearest = read_table("snp_annotations_nearest_susceptibility_imsgc",skip=30,col_types = cols(.default="c")) %>% 
    dplyr::select(`#Uploaded_variation`,NEAREST) %>% 
    distinct()

overall_res = data.frame()
# beta-beta plots
for(anc in c("sas","afr")){
dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID)

# join with IMSGC
combo_dat = dat %>%
    dplyr::select(-BP) %>%
    dplyr::rename("BP"= GENPOS_hg19) %>%
    inner_join(full_imsgc,by=c("CHR","BP")) %>%
    filter(
        (ALLELE1 == A1 & ALLELE0 == A2 ) |
        (ALLELE1 == A2 & ALLELE0 == A1 ))

# join with nearest gene (which is hg38)
combo_dat = combo_dat %>%
    mutate(`#Uploaded_variation` = paste0(CHR,":",GENPOS)) %>%
    inner_join(nearest,by=c("#Uploaded_variation"))

# flip betas to orient to ALLELE1
combo_dat = combo_dat %>%
    mutate(BETA.y = ifelse(
        (ALLELE1 == A1 & ALLELE0 == A2 ),
        BETA.y,
        BETA.y*-1))

# look at concordance
# beta-beta plot for non-MHC SNPs SNPs achieving P < 5e-8 in IMSGC

combo_dat = combo_dat %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))


overall_res <<- bind_rows(overall_res,combo_dat %>% mutate(ancestry = anc))
concordance = combo_dat %>% 
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value
cor = cor.test(combo_dat$BETA.y,combo_dat$BETA.x,method="spearman")
label = paste0(
    "One-tailed binomial test P: ",binom,"\nAncestry: ",anc,"\nConcordant SNPs: ",
    concordance[concordance$concordant=="yes",]$n,"/",concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
    "\nRho = ",round(cor$estimate,2)," P = ",cor$p.value)

print(label)

p = ggplot(
    combo_dat,
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.1)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.1)+
    geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.1)+
    geom_point()+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size",size=bquote(-log[10]~P))+
    ggtitle(label)+
    scale_color_brewer(palette="Set1")+
    ggrepel::geom_text_repel(mapping = aes(label =NEAREST),color="black",size=3)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot.png"
)


png(outplot,units="in",res=900,width=8,height=8)
print(p)
dev.off()

# repeat with only P < 0.5 SNPs
full_data = combo_dat
combo_dat = combo_dat %>%
    filter(P.x < 0.5) %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))

concordance = combo_dat %>% 
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value
cor = cor.test(combo_dat$BETA.y,combo_dat$BETA.x,method="spearman")
label = paste0(
    "One-tailed binomial test P: ",binom,"\nAncestry: ",anc,"\nConcordant SNPs: ",
    concordance[concordance$concordant=="yes",]$n,"/",concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
    "\nRho = ",round(cor$estimate,2)," P = ",cor$p.value)

print(label)

p = ggplot(
    combo_dat,
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.1)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.1)+
    geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.1)+
    geom_point()+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size",size=bquote(-log[10]~P))+
    ggtitle(label)+
    scale_color_brewer(palette="Set1")+
    ggrepel::geom_text_repel(mapping = aes(label =NEAREST),color="black",size=3)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot_precisely_estimated_snps.png"
)


png(outplot,units="in",res=900,width=8,height=8)
print(p)
dev.off()


# repeat with only P < 0.1 SNPs
combo_dat = full_data
combo_dat = combo_dat %>%
    filter(P.x < 0.1) %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))

concordance = combo_dat %>% 
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value
cor = cor.test(combo_dat$BETA.y,combo_dat$BETA.x,method="spearman")
label = paste0(
    "One-tailed binomial test P: ",binom,"\nAncestry: ",anc,"\nConcordant SNPs: ",
    concordance[concordance$concordant=="yes",]$n,"/",concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
    "\nRho = ",round(cor$estimate,2)," P = ",cor$p.value)

print(label)

p = ggplot(
    combo_dat,
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.1)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.1)+
    geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.1)+
    geom_point()+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size",size=bquote(-log[10]~P))+
    ggtitle(label)+
    scale_color_brewer(palette="Set1")+
    ggrepel::geom_text_repel(mapping = aes(label =NEAREST),color="black",size=3)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot_very_precisely_estimated_snps.png"
)


png(outplot,units="in",res=900,width=8,height=8)
print(p)
dev.off()



# repeat without MHC 

combo_dat = full_data %>%
    filter(!(
        CHR == 6 & BP > 25e6 & BP < 35e6
    )) %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))

concordance = combo_dat %>% 
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value
cor = cor.test(combo_dat$BETA.y,combo_dat$BETA.x,method="spearman")
label = paste0(
    "One-tailed binomial test P: ",binom,"\nAncestry: ",anc,"\nConcordant SNPs: ",
    concordance[concordance$concordant=="yes",]$n,"/",concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
    "\nRho = ",round(cor$estimate,2)," P = ",cor$p.value)

print("NO MHC")
print(label)

p = ggplot(
    combo_dat,
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.1)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.1)+
    geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.1)+
    geom_point()+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size",,size=bquote(-log[10]~P))+
    ggtitle(label)+
    scale_color_brewer(palette="Set1")+
    ggrepel::geom_text_repel(mapping = aes(label =NEAREST),color="black",size=3)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot_NOMHC.png"
)


png(outplot,units="in",res=900,width=8,height=8)
print(p)
dev.off()


# no MHC, precisely estimated 

# repeat without MHC 

combo_dat = full_data %>%
    filter(P.x < 0.5) %>%
    filter(!(
        CHR == 6 & BP > 25e6 & BP < 35e6
    )) %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))

concordance = combo_dat %>% 
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value
cor = cor.test(combo_dat$BETA.y,combo_dat$BETA.x,method="spearman")
label = paste0(
    "One-tailed binomial test P: ",binom,"\nAncestry: ",anc,"\nConcordant SNPs: ",
    concordance[concordance$concordant=="yes",]$n,"/",concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
    "\nRho = ",round(cor$estimate,2)," P = ",cor$p.value)

print("NO MHC")
print(label)

p = ggplot(
    combo_dat,
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.1)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.1)+
    geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.1)+
    geom_point()+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size",,size=bquote(-log[10]~P))+
    ggtitle(label)+
    scale_color_brewer(palette="Set1")+
    ggrepel::geom_text_repel(mapping = aes(label =NEAREST),color="black",size=3)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot_NOMHC_precise.png"
)


png(outplot,units="in",res=900,width=8,height=8)
print(p)
dev.off()



}


wide_dat = overall_res %>% 
    mutate(OR_ADAMS = exp(BETA.x)) %>%
    mutate(OR_IMSGC = exp(BETA.y)) %>%
    dplyr::select(CHROM,GENPOS,SNP.y,ALLELE1,OR_ADAMS,P.x,OR_IMSGC,P.y,NEAREST,ancestry,contains("FREQ"),concordant) %>% 
    pivot_wider(id_cols = c(CHROM,GENPOS,ALLELE1,SNP.y,NEAREST,OR_IMSGC,P.y),
    values_from = c(OR_ADAMS,P.x,A1FREQ_CONTROLS,A1FREQ_CASES,A1FREQ,concordant),
    names_from = ancestry) %>% 
    dplyr::rename("SNP"=SNP.y)

# add reference frequencies 
ref_freqs = read_table("/data/scratch/hmy117/snp_annotations_freqs_susceptibility_imsgc",skip=30) %>% 
    dplyr::rename("SNP"=`#Uploaded_variation`,"ALLELE1"=Allele) %>% 
    dplyr::select(-2) %>% 
    filter(gnomADg_NFE_AF != "-") %>% 
    distinct(SNP,gnomADg_NFE_AF,ALLELE1,.keep_all=T)

wide_dat %>% 
    inner_join(ref_freqs,by=c("SNP","ALLELE1"))

write_csv(wide_dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/concordance_imsgc_eur_hits.csv")
