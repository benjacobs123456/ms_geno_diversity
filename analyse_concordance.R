library(tidyverse)

setwd("/data/scratch/hmy117")
# read in clumped imsgc snps
clumps = list()
for(i in c(1:22)){
    file = paste0("clumped_imsgc_chr",i,".clumped")
    clumps[[i]] = if(file.exists(file)){
    read_table(file)
    }
}
imsgc = do.call("bind_rows",clumps)

full_imsgc = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv",col_types = "dcccddddddc")
full_imsgc = full_imsgc %>% filter(rsID %in% imsgc$SNP)

nearest = read_table("snp_annotations_nearest_susceptibility_imsgc",skip=31,col_types = cols(.default="c")) %>%
dplyr::select(`#Uploaded_variation`,NEAREST) %>%
distinct() %>% 
dplyr::rename("rsID" = `#Uploaded_variation`)
full_imsgc = full_imsgc %>% left_join(nearest,by="rsID")
overall_res = data.frame()

# beta-beta plots
for(anc in c("sas","afr","eur")){

            dat = read_table(paste0(
                "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
                anc,".MS_status.glm.logistic.hybrid_hg38"
            )) %>%
                mutate(A2 = ifelse(A1==REF,ALT,REF))%>%
                dplyr::rename("CHR" = `#CHROM`,"BP"=POS,"SNP"=ID,"ALLELE1"=A1,"ALLELE0"=A2) %>%
                mutate(BETA=log(OR))

            # add freqs from regenie
            dat_regenie = read_table(paste0(
                "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
                anc,"_MS_status.regenie_hg38"
            )) %>%
                mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>% 
                dplyr::select(SNP,A1FREQ_CASES,A1FREQ_CONTROLS)


            # join 
            dat = dat %>% left_join(dat_regenie,by="SNP")

            # join with IMSGC
                combo_dat = dat %>%
                    inner_join(full_imsgc,by=c("CHR","BP")) %>%
                    filter(
                        (ALLELE1 == A1 & ALLELE0 == A2 ) |
                        (ALLELE1 == A2 & ALLELE0 == A1 ))


            # flip betas to orient to ALLELE1
            combo_dat = combo_dat %>%
                mutate(BETA.y = ifelse(
                    (ALLELE1 == A1 & ALLELE0 == A2 ),
                    BETA.y,
                    BETA.y*-1))

        # look at concordance
        # beta-beta plot for 

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
png(outplot,res=900,units="in",width=8,height=8)
print(p)
dev.off()

# w/o mhc 
combo_dat = combo_dat %>%
    filter(!
    (CHR==6 & BP > 25e6 & BP <35e6)
    )

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
    anc,"_MS_status_beta_beta_plot_nomhc.png"
)
png(outplot,res=900,units="in",width=8,height=8)
print(p)
dev.off()

# sensitivity analysis 

}

# heatmap 
overall_res = overall_res %>%
mutate(MHC = ifelse(
    CHR == 6 & BP > 25e6 & BP < 35e6,"MHC","Non-MHC"
)) %>%
mutate(snpgene = paste0(rsID," (",NEAREST,")"))
overall_res = overall_res %>% arrange(CHR,BP)
overall_res$snpgene = factor(overall_res$snpgene,levels = unique(overall_res$snpgene),ordered=T)

# p val
overall_res = overall_res %>% mutate(
    pval_simple = case_when(
        P.x < 5e-8 ~ "***",
        P.x < 1e-5 ~ "**",
        P.x < 0.05 ~ "*",
        P.x >= 0.05 ~ " "
    )
)
p = ggplot(overall_res%>%
filter(P.y < 5e-8),aes(toupper(ancestry),snpgene,fill=concordant,label=pval_simple))+
geom_tile(color="black") +
scale_fill_brewer(palette="Set1")+
geom_text(alpha=0.7)+
labs(fill="Concordance",x="Ancestry",y="IMSGC risk SNP and nearest gene")+
theme_bw()

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/concordance_summary.png",res=900,units="in",width=6,height=8)
print(p)
dev.off()

# edit this


wide_dat = overall_res %>%
    mutate(OR_ADAMS = exp(BETA.x)) %>%
    mutate(OR_IMSGC = exp(BETA.y)) %>%
    dplyr::select(CHR,BP,SNP.x,rsID,ALLELE1,OR_ADAMS,P.x,OR_IMSGC,P.y,NEAREST,ancestry,contains("FREQ"),concordant) %>%
    distinct(rsID,ancestry,.keep_all=T) %>%
    pivot_wider(id_cols = c(CHR,BP,ALLELE1,SNP.x,NEAREST,OR_IMSGC,P.y,rsID),
    values_from = c(OR_ADAMS,P.x,A1FREQ_CONTROLS,A1FREQ_CASES,concordant),
    names_from = ancestry) %>%
    dplyr::rename("SNP"=SNP.x)

# add reference frequencies
ref_freqs = read_table("/data/scratch/hmy117/snp_annotations_freqs_susceptibility_imsgc",skip=31) %>%
    dplyr::rename("SNP"=`#Uploaded_variation`,"ALLELE1"=Allele) %>%
    dplyr::select(-2) %>%
    filter(gnomADg_NFE_AF != "-") %>%
    distinct(SNP,.keep_all=T)

wide_dat = wide_dat %>%
    inner_join(ref_freqs %>%
    dplyr::rename("rsID"=SNP),by=c("rsID"))

write_csv(wide_dat,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/concordance_imsgc_eur_hits.csv")
