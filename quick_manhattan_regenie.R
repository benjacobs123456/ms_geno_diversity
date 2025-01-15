library(tidyverse)
library(qqman)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
anc = args[1]


# annotations
nearest = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_nearest_susceptibility_",anc),skip=30, col_types = cols(.default = "c"))
colnames(nearest)[1] = "SNP"
nearest = nearest %>% distinct(SNP, NEAREST,.keep_all=T)

# regenie

dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    left_join(nearest,by="SNP")
print("N SNPs:")
print(nrow(dat))


# read in frequencies
freqs = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_freqs_susceptibility_",anc),skip=56, col_types = cols(.default = "c"))
colnames(freqs)[1] = "SNP"
freqs = freqs %>% distinct(SNP, .keep_all=T)
freqs = freqs %>% dplyr::select(SNP,contains("AF"),Allele)
dat = dat %>%  left_join(freqs,by="SNP")

gnomad_anc_code = ifelse(anc == "afr","AFR","SAS")
freq_code = paste0("gnomADg_",gnomad_anc_code,"_AF")


# plot
outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.png"
)

# compare with imsgc
imsgc = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv",col_types = "dccdddddddc")

# sig hits
imsgc_sighits = imsgc %>% filter(P < 1e-5)

# combine
dat = dat %>% mutate(SNP = paste0(CHROM,":",GENPOS))
dat = dat %>%
    mutate(sig_imsgc = ifelse(SNP %in% imsgc_sighits$SNP,"sig","not_sig"))

# add OR
dat = dat %>% 
    mutate(OR = exp(BETA)) %>% 
    mutate(lower_ci = exp(BETA - 1.96*SE)) %>%
    mutate(upper_ci = exp(BETA - 1.96*SE)) 

# manhattan
chrcoords = dat %>% group_by(CHROM) %>% summarise(min_bp = min(GENPOS), max_bp = max(GENPOS)) %>%
    mutate(CHROM = CHROM + 1) %>%
    mutate(cumbp_total = cumsum(max_bp))
dat = dat %>%
    left_join(chrcoords,by="CHROM") %>%
    mutate(cumbp = ifelse(is.na(cumbp_total),GENPOS,GENPOS+cumbp_total))
midpoints = dat %>%
    group_by(CHROM) %>%
    summarise(midpoint = median(cumbp))


# define windows
dat = dat %>% mutate(window = round(cumbp / 1e6,0))
sig_imsgc_windows = dat %>% group_by(window) %>% dplyr::count(sig_imsgc) %>% filter(n>1) %>%
    filter(sig_imsgc=="sig")
dat = dat %>%
    mutate(colcode = case_when(
        window %in% sig_imsgc_windows$window ~ "sig",
        CHROM %%2 == 0 ~ "even",
        CHROM %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

p = ggplot(dat %>% filter(P < 0.5),aes(cumbp,LOG10P,color=colcode))+
    geom_point(data = dat %>% filter(colcode != "sig"),alpha=0.6)+
    geom_point(data = dat %>% filter(colcode == "sig"),alpha=0.8)+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHROM)+
    scale_color_manual(values = pal)+
    ggrepel::geom_text_repel(data = dat %>%
        filter(P < 1e-4 & colcode == "sig") %>%
        filter(!is.na(NEAREST)) %>%
        group_by(CHROM) %>%
        slice_min(P,with_ties=F,n=3) %>%
        distinct(NEAREST,.keep_all=T),
    mapping = aes(label = NEAREST),min.segment.length=0,nudge_y=0.5,color="black")+
    theme_bw()+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")

png(outplot,units="in",res=900,width=10,height=4)
p
dev.off()

# save to file
outfile = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_sig_hits_with_freqs.csv"
)

# add QC data
snp_qc_data = read_tsv(paste0("/data/scratch/hmy117/info_stats_all_snps_",anc,".tsv"))

# combine
dat = dat %>%
    left_join(snp_qc_data %>% dplyr::select(CHR,POS,R2) %>%
    dplyr::rename("GENPOS_hg19"=POS,"CHROM"=CHR),by=c("CHROM","GENPOS_hg19"))

median(dat[dat$R2==1,]$CHISQ)
sig_hits_with_imsgc = dat %>% 
filter(P < 1e-4) %>%
    left_join(imsgc %>% 
    dplyr::rename("OR_IMSGC" = OR,"P_IMSGC"=P), 
        by="SNP")

write_csv(sig_hits_with_imsgc,
outfile)



# AF plot
# plot
outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_gnomad_maf_",
    anc,"_MS_status.png"
)
sig_hits = dat %>% filter(P < 1e-4)

png(outplot,units="in",res=900,width=4,height=4)
ggplot(sig_hits,aes(A1FREQ_CONTROLS,as.numeric(.data[[freq_code]])))+
    geom_point()+
    labs(x="Control AF",y="gnomAD AF")+
    theme_bw()+
    geom_abline(slope=1,intercept=0)
dev.off()

# Beta vs MAF
outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_beta_maf_",
    anc,"_MS_status.png"
)

png(outplot,units="in",res=900,width=4,height=4)
ggplot(sig_hits,aes(A1FREQ,abs(BETA)))+
    geom_point()+
    labs(x="AF",y="|Beta|")+
    theme_bw()
dev.off()

# Info plot

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_info_maf_p_",
    anc,"_MS_status.png"
)

png(outplot,units="in",res=900,width=4,height=4)
ggplot(sig_hits,aes(R2,LOG10P))+
    geom_point()+
    labs(x=bquote(Imputation~R^2),y=bquote(-log[10]~P))+
    theme_bw()
dev.off()


# manhattan
dat = dat %>%
    mutate(colcode = case_when(
        P < 1e-4 ~ "sig",
        P >=1e-4 & CHROM %%2 == 0 ~ "even",
        P >=1e-4 & CHROM %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

# define windows
dat = dat %>% mutate(window = round(cumbp / 1e6,0))
sig_windows = dat %>% group_by(window) %>% 
    dplyr::count(sighits = colcode=="sig") %>% 
    filter(sighits==T)
dat = dat %>% mutate(colcode = ifelse(window %in% sig_windows$window,"sig",colcode))

dat_for_gene_labels = data = dat %>%
        filter(P < 1e-4) %>%
        filter(!is.na(NEAREST)) %>%
        group_by(window) %>%
        slice_min(P,with_ties=F) %>%
        distinct(NEAREST,.keep_all=T)

near_imsgc_hist = list()
for(i in c(1:nrow(dat_for_gene_labels))){
    this_snp = dat_for_gene_labels[i,]
    this_imsgc_locus = imsgc %>% filter(CHR == this_snp$CHR & BP > this_snp$BP - 1e6 & BP < this_snp$BP + 1e6)
    hits = this_imsgc_locus %>% filter( P < 1e-5)
    near_imsgc_hist[[i]] = nrow(hits) > 0
}
dat_for_gene_labels$near_imsgc_hit = unlist(near_imsgc_hist)
p = ggplot(dat %>% filter(P < 0.1),aes(cumbp,LOG10P,color=colcode))+
    geom_point(data = dat %>% filter(P < 0.1 & colcode != "sig"),alpha=0.7)+
    geom_point(data = dat %>% filter(P < 0.1 & colcode == "sig" & P >1e-4),alpha=0.7)+
    geom_point(data = dat %>% filter(P < 0.1 & colcode == "sig" & P <=1e-4),alpha=1)+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHROM)+
    scale_color_manual(values = pal)+
    ggrepel::geom_label_repel(data = dat_for_gene_labels,
    mapping = aes(label = NEAREST,fill=near_imsgc_hit),min.segment.length=0,nudge_y=1,color="black",segment.linetype = 5,segment.alpha=0.6)+
    theme_bw()+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")+
    geom_hline(yintercept=4,linetype="dashed",alpha=0.5,color="blue")+
    geom_hline(yintercept=-log10(5e-8),linetype="dashed",alpha=0.5,color="red")

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_main.png"
)


png(outplot,units="in",res=900,width=12,height=4)
p
dev.off()

# beta-beta plots
dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    left_join(nearest,by="SNP")

# join
combo_dat = dat %>%
    inner_join(imsgc,by=c("CHR","BP")) %>%
    filter(
        (ALLELE1 == A1 & ALLELE0 == A2 ) |
        (ALLELE1 == A2 & ALLELE0 == A1 ))

# flip betas
combo_dat = combo_dat %>%
    mutate(BETA.x = ifelse(
        (ALLELE1 == A1 & ALLELE0 == A2 ),
        BETA.x,
        BETA.x*-1))

# look at concordance
# beta-beta plot for non-MHC SNPs SNPs achieving P < 5e-8 in IMSGC

combo_dat = combo_dat %>%
    mutate(concordant = ifelse(sign(BETA.x) == sign(BETA.y),
    "yes","no"))

concordance = combo_dat %>% filter(!(CHROM == 6 & GENPOS > 25e6 & GENPOS <35e6) &
    P.y < 5e-8 & P.x < 0.05) %>%
    dplyr::count(concordant)


binom = binom.test(concordance[concordance$concordant=="yes",]$n,
concordance[concordance$concordant=="yes",]$n + concordance[concordance$concordant=="no",]$n,
0.5,alternative="greater")$p.value

p = ggplot(
    combo_dat %>% filter(!(CHROM == 6 & GENPOS > 25e6 & GENPOS <35e6) & P.y < 5e-8 & P.x < 0.05),
    aes(BETA.y,BETA.x,size=-log10(P.x),col=concordant))+
    geom_point()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    geom_errorbarh(mapping = aes(xmin = BETA.y - 1.96* SE.y,xmax = BETA.y + 1.96* SE.y,BETA.x),height=0.1,linewidth=0.1)+
    geom_errorbar(mapping = aes(ymin = BETA.x - 1.96* SE.x,ymax = BETA.x + 1.96* SE.x,BETA.y),width=0.1,linewidth=0.1)+
    theme_bw()+
    labs(x="IMSGC effect size",y="ADAMS effect size")+
    ggtitle("One-tailed binomial test P: ",binom)

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_beta_beta_plot.png"
)

png(outplot,units="in",res=900,width=6,height=6)
p
dev.off()

# PLINK glm vs regenie

dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    left_join(nearest,by="SNP")

plink_dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
    anc,".MS_status.glm.logistic.hybrid"
))
plink_dat = plink_dat %>% dplyr::rename("CHR" = `#CHROM`,"BP"=POS) %>%
    mutate(BETA = log(OR))

# join
combo_dat = dat %>%
    left_join(plink_dat,by="ID")

# beta-beta plot & P-P plot vs plink
p = ggplot(
    combo_dat %>% filter(P.x < 0.5 & P.y < 0.5),
    aes(BETA.y,BETA.x))+
    geom_point()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme_bw()+
    geom_abline(intercept=0,slope=1)+
    labs(x="PLINK effect size",y="REGENIE effect size")
p1 = ggplot(
    combo_dat %>% filter(P.x < 0.5 & P.y < 0.5),
    aes(-log10(P.y),-log10(P.x),col=factor(CHROM)))+
    geom_point()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme_bw()+
    geom_abline(intercept=0,slope=1)+
    labs(x="PLINK -log10(P)",y="REGENIE -log10(P)",col="Chromosome")

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_plink_vs_regenie.png"
)

png(outplot,units="in",res=900,width=10,height=4)
cowplot::plot_grid(p,p1,ncol=2)
dev.off()

# QQ plots
# inflation
dat$P = 10^-dat$LOG10P
lambda = median(qchisq(dat$P,df=1,lower.tail=F)) /  qchisq(0.5,lower.tail=F,df=1)
message(lambda)
pval_dat = data.frame(observed = -log10(dat$P)) %>%
    arrange(observed)

pval_dat = pval_dat %>%
    mutate(expected = sort(-log10(runif(n = nrow(pval_dat)))) )

# plot
p = ggplot(pval_dat,aes(expected,observed))+
    geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+
    geom_line()+
    labs(x="Expected -log10(P)",y="Observed -log10(P)")+
    ggtitle(paste0("Ancestry: ",toupper(anc),"\nLambda =  ",round(lambda,2)))+
    theme_bw()

png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/qq_plot_susceptibility_",anc,".png"),res=900,units="in",width=4,height=4)
print(p)
dev.off()

# manhattan for plink glm
plink_dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
    anc,".MS_status.glm.logistic.hybrid"
))
plink_dat = plink_dat %>% dplyr::rename("CHR" = `#CHROM`,"BP"=POS,"SNP"=ID)

# manhattan
chrcoords = plink_dat %>% group_by(CHR) %>% summarise(min_bp = min(BP), max_bp = max(BP)) %>%
    mutate(CHR = CHR + 1) %>%
    mutate(cumbp_total = cumsum(max_bp))
plink_dat = plink_dat %>%
    left_join(chrcoords,by="CHR") %>%
    mutate(cumbp = ifelse(is.na(cumbp_total),BP,BP+cumbp_total))
midpoints = plink_dat %>%
    group_by(CHR) %>%
    summarise(midpoint = median(cumbp))

plink_dat = plink_dat %>%
    mutate(colcode = case_when(
        P < 1e-5 ~ "sig",
        P >=1e-5 & CHR %%2 == 0 ~ "even",
        P >=1e-5 & CHR %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

# define windows
plink_dat = plink_dat %>% mutate(window = round(cumbp / 1e6,0))
p = ggplot(plink_dat %>% filter(P < 0.1),aes(cumbp,-log10(P),color=colcode))+
    geom_point()+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHR)+
    scale_color_manual(values = pal)+
    ggrepel::geom_text_repel(data = plink_dat %>%
        filter(P < 1e-5) %>%
        group_by(CHR) %>%
        slice_min(P,with_ties=F),
    mapping = aes(label = SNP),min.segment.length=0,nudge_y=0.5,color="black")+
    theme_bw()+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_plink_",
    anc,"_MS_status.png"
)

png(outplot,units="in",res=900,width=10,height=4)
p
dev.off()

# inflation
dat= plink_dat %>% filter(!is.na(P))
lambda = median(qchisq(dat$P,df=1,lower.tail=F)) /  qchisq(0.5,lower.tail=F,df=1)
message(lambda)
pval_dat = data.frame(observed = -log10(dat$P)) %>%
    arrange(observed)

pval_dat = pval_dat %>%
    mutate(expected = sort(-log10(runif(n = nrow(pval_dat)))) )

# plot
p = ggplot(pval_dat,aes(expected,observed))+
    geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+
    geom_line()+
    labs(x="Expected -log10(P)",y="Observed -log10(P)")+
    ggtitle(paste0("Ancestry: ",toupper(anc),"\nLambda =  ",round(lambda,2)))+
    theme_bw()

png(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/qq_plot_susceptibility_PLINK_",anc,".png"),res=900,units="in",width=4,height=4)
print(p)
dev.off()
