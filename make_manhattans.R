library(tidyverse)
library(qqman)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
anc = args[1]

# annotations
nearest = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_nearest_susceptibility_",anc),skip=31, col_types = cols(.default = "c"))
colnames(nearest)[1] = "SNP"
nearest = nearest %>% distinct(SNP, NEAREST,.keep_all=T)

# PLINK GWAS
plink_dat = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
    anc,".MS_status.glm.logistic.hybrid_hg38"
)) %>%
    dplyr::rename("CHROM" = `#CHROM`,"GENPOS"=POS,"SNP"=ID) %>%
    left_join(nearest,by="SNP")
print("N SNPs:")
print(nrow(plink_dat))

# calculate lambda
lambda = median(qchisq(plink_dat$P,df=1,lower.tail=F)) /  qchisq(0.5,lower.tail=F,df=1)
message(lambda)

# QQ plots
# inflation
pval_dat = data.frame(observed = -log10(plink_dat$P)) %>%
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

# adjust P value for GC
plink_dat$P_unadjusted = plink_dat$P
plink_dat$P = pchisq(plink_dat$Z_STAT^2/lambda, df = 1, lower = F)

# hg19 names
plink_dat_hg19 = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_glm",
    anc,".MS_status.glm.logistic.hybrid"
)) %>%
    dplyr::rename("CHROM" = `#CHROM`,"POS_hg19"=POS,"SNP_hg19"=ID) %>%
    dplyr::select(CHROM,POS_hg19,SNP_hg19)


# clump
plink_dat = plink_dat %>%
  left_join(plink_dat_hg19,by=c("CHROM","POS_hg19"))


# write SNP list for clumping
readr::write_tsv(
plink_dat %>% filter(P_unadjusted<1e-5) %>%
dplyr::select(SNP_hg19,P_unadjusted) %>%
dplyr::rename("SNP"=SNP_hg19,"P"=P_unadjusted),
paste0("/data/scratch/hmy117/snps_for_clump_",anc))

cmd = paste0(
"~/plink --bfile /data/scratch/hmy117/imputed_risk_gwas_genotypes_step2_imputation_post_qc_",anc,"_ALL_CHRS ",
"--clump /data/scratch/hmy117/snps_for_clump_",anc," ",
"--clump-p1 1e-5 --clump-r2 0.001 --clump-kb 1000 --out /data/scratch/hmy117/clumps_",anc)
system(cmd)

# read back in
clumps = read_table(paste0("/data/scratch/hmy117/clumps_",anc,".clumped"))
plink_dat = plink_dat %>%
  mutate(tophit = ifelse(SNP_hg19 %in% clumps$SNP,"tophit"," "))
message(nrow(clumps))

# read in frequencies
freqs = read_table(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/snp_annotations_freqs_susceptibility_",anc),skip=57, col_types = cols(.default = "c"))
colnames(freqs)[1] = "SNP"
freqs = freqs %>% distinct(SNP, .keep_all=T)
freqs = freqs %>% dplyr::select(SNP,contains("AF"),Allele)
plink_dat = plink_dat %>%  left_join(freqs,by="SNP")

gnomad_anc_code = if(anc == "afr"){
    "AFR"
} else if(anc=="sas") {
    "SAS"
} else {
    "NFE"
}

freq_code = paste0("gnomADg_",gnomad_anc_code,"_AF")

# compare with imsgc
imsgc = read_tsv("/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_ms_risk_discovery_hg38.tsv",col_types = "dcccddddddc")

# sig hits
imsgc_sighits = imsgc %>% filter(P < 1e-5)

# define chr:pos ID in hg38
dat = plink_dat %>% mutate(SNP = paste0(CHROM,":",GENPOS))

# define whether SNP is a suggestive hit in IMSGC
dat = dat %>%
    mutate(sig_imsgc = ifelse(SNP %in% imsgc_sighits$SNP,"sig","not_sig"))

# add OR
dat = dat %>%
    mutate(BETA = log(OR)) %>%
    mutate(lower_ci = exp(BETA - 1.96*`LOG(OR)_SE` )) %>%
    mutate(upper_ci = exp(BETA + 1.96*`LOG(OR)_SE` ))


# add QC data
snp_qc_data = read_tsv(paste0("/data/scratch/hmy117/info_stats_all_snps_",anc,".tsv"))

# combine
dat = dat %>%
    left_join(snp_qc_data %>% dplyr::select(CHR,POS,R2) %>%
    dplyr::rename("POS_hg19"=POS,"CHROM"=CHR),by=c("CHROM","POS_hg19"))


sig_hits_with_imsgc = dat %>%
    filter(P_unadjusted < 1e-5) %>%
    left_join(imsgc %>%
    dplyr::rename("OR_IMSGC" = OR,"P_IMSGC"=P),
        by="SNP")

# add REGENIE GWAS
dat_regenie = read_table(paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status.regenie_hg38"
)) %>%
    mutate(P = 10^-LOG10P, CHR = CHROM, BP = GENPOS, SNP = ID) %>%
    mutate(SNP = paste0(CHROM,":",GENPOS))

# join
sig_hits_with_imsgc = sig_hits_with_imsgc %>%
  dplyr::select(-`FIRTH?`,-OBS_CT,-TEST,-Z_STAT,-ERRCODE,-POS_hg19,-Location,-Allele.x,-Gene,-AF,-gnomADg_AF,-gnomADg_AMI_AF,-gnomADg_AMR_AF,-gnomADg_EAS_AF,-gnomADg_ASJ_AF,-gnomADg_FIN_AF,-gnomADg_MID_AF,-gnomADg_REMAINING_AF,-BP,-BETA.y,-rsID,-`LOG(OR)_SE`) %>%
  dplyr::rename("P_PLINK"=P,"OR_PLINK"=OR) %>%
  left_join(
    dat_regenie %>%
    dplyr::select(SNP,ALLELE0,ALLELE1,A1FREQ,A1FREQ_CASES,A1FREQ_CONTROLS,N_CASES,N_CONTROLS,BETA,SE,P) %>%
    dplyr::rename("P_REGENIE" = P),
    by="SNP") %>%
  dplyr::select(-Allele.y,-BETA.x,-lower_ci,-upper_ci,-CHR,-N,-Z,-SE.x) %>%
  mutate(OR_REGENIE = exp(BETA)) %>%
  dplyr::select(-BETA,-SE.y)

# filter to top hits
sig_hits_with_imsgc = sig_hits_with_imsgc %>% filter(tophit=="tophit") %>%
dplyr::select(-tophit)

# scan through and look within window to see if any sig imsgc hits
imsgc_near_hit_snp = list()
imsgc_near_hit_dist = list()
imsgc_near_hit_p = list()

for(i in c(1:nrow(sig_hits_with_imsgc))){
message(i)
this_snp = sig_hits_with_imsgc[i,]

# filter IMSGC to this window
near_imsgc_hit = imsgc %>%
  filter(CHR == this_snp$CHROM & BP > (this_snp$GENPOS - 1e6) & BP < (this_snp$GENPOS + 1e6) ) %>%
  slice_min(P,with_ties=F)

imsgc_near_hit_snp[[i]] = near_imsgc_hit$SNP[1]
imsgc_near_hit_dist[[i]] = abs(near_imsgc_hit$BP[1] - this_snp$GENPOS)
imsgc_near_hit_p[[i]] = near_imsgc_hit$P[1]
}
sig_hits_with_imsgc$imsgc_hit_snp = unlist(imsgc_near_hit_snp)
sig_hits_with_imsgc$imsgc_hit_snp_dist = unlist(imsgc_near_hit_dist)
sig_hits_with_imsgc$imsgc_hit_snp_p = unlist(imsgc_near_hit_p)

# save to file
outfile = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_sig_hits_with_freqs.csv"
)

write_csv(sig_hits_with_imsgc,
outfile)



sig_hits = dat %>% filter(P < 5e-8)


# manhattan
dat = dat %>%
    mutate(colcode = case_when(
        P < 5e-8 ~ "sig",
        P >=5e-8 & CHROM %%2 == 0 ~ "even",
        P >=5e-8 & CHROM %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

# define windows
chrcoords = dat %>% group_by(CHROM) %>% summarise(min_bp = min(GENPOS), max_bp = max(GENPOS)) %>%
    mutate(CHROM = CHROM + 1) %>%
    mutate(cumbp_total = cumsum(max_bp))
dat = dat %>%
    left_join(chrcoords,by="CHROM") %>%
    mutate(cumbp = ifelse(is.na(cumbp_total),GENPOS,GENPOS+cumbp_total))
midpoints = dat %>%
    group_by(CHROM) %>%
    summarise(midpoint = median(cumbp))

dat = dat %>% mutate(window = round(cumbp / 5e5,0))
sig_windows = dat %>% group_by(window) %>%
    dplyr::count(sighits = colcode=="sig") %>%
    filter(sighits==T)
dat = dat %>% mutate(colcode = ifelse(window %in% sig_windows$window,"sig",colcode))

dat_for_gene_labels = dat %>%
        filter(P < 5e-8) %>%
        filter(!is.na(NEAREST)) %>%
        group_by(CHROM,NEAREST) %>%
        slice_min(P,with_ties=F) %>%
        ungroup() %>%
        group_by(CHROM) %>%
        distinct(NEAREST,.keep_all=T)


p = ggplot(dat %>% filter(P < 0.1),aes(cumbp,-log10(P),color=colcode))+
    geom_point(data = dat %>% filter(P < 0.1 & colcode != "sig"),alpha=0.3)+
    geom_point(data = dat %>% filter(P < 0.1 & colcode == "sig" & P >5e-8),alpha=0.5)+
    geom_point(data = dat %>% filter(P < 0.1 & colcode == "sig" & P <=5e-8),alpha=1)+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHROM)+
    scale_color_manual(values = pal)+
    ggrepel::geom_label_repel(data = dat_for_gene_labels,
    mapping = aes(label = NEAREST),min.segment.length=0,nudge_y=1,color="black",segment.linetype = 5,segment.alpha=0.6)+
    theme_bw()+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")+
    geom_hline(yintercept=-log10(1e-5),linetype="dashed",alpha=0.5,color="blue")+
    geom_hline(yintercept=-log10(5e-8),linetype="dashed",alpha=0.5,color="red")

outplot = paste0(
    "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/case_control_gwas_ALL_ANCESTRY_",
    anc,"_MS_status_main.png"
)


png(outplot,units="in",res=900,width=12,height=4)
p
dev.off()
