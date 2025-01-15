
library(tidyverse)


#######################
# TWO DIGIT RESOLUTION
#######################

make_hla_forest = function(ancestry,allele,model="additive"){
  message(ancestry)
  message(allele)

  # nantes
  in_file = paste0("/data/scratch/hmy117/hla_allele_nantes_",allele,"_",ancestry,".tsv")
  # 1KG
  # in_file = paste0("/data/scratch/hmy117/hla_allele_kg_",allele,"_",ancestry,".tsv")

  # step 1 - compare frequencies & get global freqs
  allele_calls = read_tsv(in_file,col_types="cccdd")
  colnames(allele_calls) = c("IID","allele1","allele2","prob","match")

  # truncate all calls to 1st field 
  allele_calls = allele_calls %>% 
    separate(allele1,sep=":",into=c("allele1","extra")) %>% 
    separate(allele2,sep=":",into=c("allele2","extra2")) %>% 
    dplyr::select(-extra,-extra2)

  allele_counts = allele_calls %>%
    dplyr::select(c(1:3)) %>%
    pivot_longer(cols = c(2,3)) %>%
    dplyr::select(-name) %>%
    group_by(IID) %>%
    dplyr::count(value)

    # join with phenot & covars

  pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv"),col_types = cols(.default="c"))
  cov = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv"),col_types = "ccdddddddddddd")

  allele_counts = allele_counts %>%
    ungroup() %>%
    left_join(pheno,by="IID") %>%
    left_join(cov,by="IID")

  # get disease-stratified AFs
  all_allele_freqs = allele_counts %>%
    filter(!is.na(MS_status)) %>%
    group_by(value,MS_status) %>%
    summarise(ac = sum(n)) %>%
    ungroup() %>%
    mutate(MS_status = ifelse(MS_status == 1, "Control","MS")) %>%
    group_by(MS_status) %>%
    mutate(total = sum(ac)) %>%
    mutate(af = ac/total) %>%
    ungroup() %>%
    pivot_wider(id_cols = value,names_from  = MS_status,values_from = c(ac,total,af)) %>%
    mutate(delta_freq = abs(af_MS - af_Control)) %>%
    mutate(log2_fold_enrichment = log2(af_MS/af_Control)) %>%
    arrange(desc(delta_freq))


  all_allele_freqs$gene = allele
  all_allele_freqs$anc = ancestry

  # add total af 
  all_allele_freqs  = all_allele_freqs %>% 
    mutate(AF = (ac_Control + ac_MS) / (total_Control + total_MS))


  all_allele_freqs = all_allele_freqs %>%
      dplyr::rename("allele"= value) 

  common_alleles = all_allele_freqs %>% filter(AF >= 0.01)

  # regress vs MS status for each allele
  alleles = unique(common_alleles$allele)

  all_res = list()
  for(i in c(1:length(alleles))){
    message(i)
    this_allele = alleles[i]

    # filter
    this_allele_counts = allele_counts %>% filter(value == this_allele)

    
    # get people with 0 counts
    this_allele_counts_zeroes = allele_counts %>% filter(!IID %in% this_allele_counts$IID)  %>%
    distinct(IID,.keep_all=T) %>%
    mutate(n = 0)

    all_dat = this_allele_counts %>%
      bind_rows(this_allele_counts_zeroes)

    all_dat$MS_status = relevel(factor(all_dat$MS_status),ref="1")

    # fisher's exact test (dominanat coding)
    tbl = table(all_dat$MS_status,all_dat$n>=1)
    fish = fisher.test(tbl)
    fish_p = fish$p.value
    names(fish_p) = "fish_p"
    fish_or = fish$estimate 

    # dominant & recessive coding
    all_dat = all_dat %>%
        mutate(n_rec = ifelse(n <=1 ,0, 1)) %>%
        mutate(n_dom = ifelse(n >=1 ,1, 0))

    # make all models and compare

    hla_model = if(model == "additive"){
      glm(data=all_dat,factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + n, family=binomial(link="logit"))
    } else if(model=="recessive"){
      glm(data=all_dat, factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  n_rec, family=binomial(link="logit"))
    } else if(model=="dominant"){
      glm(data=all_dat, factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + n_dom, family=binomial(link="logit"))
    }


    if(nrow(summary(hla_model)$coefficients)==13){
      all_res[[i]] = c(summary(hla_model)$coefficients[13,],fish_p,fish_or)
    } else {
      na_vec = c(NA,NA,NA,NA,NA,NA)
      names(na_vec) = colnames(summary(hla_model)$coefficients)
      all_res[[i]] = na_vec
    }
  }
  all_res = do.call("bind_rows",all_res)
  colnames(all_res) = c("beta","se","z","p","fisher_p","fisher_or")
  all_res$allele = alleles

  all_res = all_res %>%
    inner_join(all_allele_freqs,
    by="allele")

  overall_res_table[[length(overall_res_table)+1]] <<- all_res

  all_res = all_res %>% arrange(desc(beta))

  all_res$hla_allele = factor(all_res$allele,ordered=T,levels = all_res$allele)
  outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_allele_two_digit_",allele,"_",ancestry,".png")
  p = ggplot(all_res,aes(beta,hla_allele))+
    geom_errorbarh(mapping = aes(xmin = beta - 1.96*se,xmax = beta +1.96*se,y=hla_allele),height=0.1)+
    geom_point(size=3,shape=21,fill="orange")+
    geom_vline(xintercept=0,linetype="dashed",color="pink")+
    labs(y=paste0(allele," allele"),x="Beta (log Odds Ratio)")+
    theme_bw()+
    ggtitle(paste0(toupper(ancestry),"\n",allele))


  png(outfile,res=900,units="in",width=4,height=4)
  print(p)
  dev.off()
}

# make parameter table
overall_res_table = list()
overall_allele_freqs = list()
alleles = c("A","B","C","DPB1","DQB1","DRB1")
ancestries = c("sas","afr")
param_tbl = expand.grid(alleles,ancestries)

## ADDITIVE
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i])
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",allele))

# n unique genes

overall_res_table_additive = overall_res_table %>%
  mutate(OR = exp(beta), lower_ci = exp(beta - 1.96*se),upper_ci = exp(beta+1.96*se)) 

overall_res_table_additive = 
bind_rows(
  overall_res_table_additive %>% 
  filter(anc=="sas") %>%
  mutate(fdr = p.adjust(p,method="fdr")),
overall_res_table_additive %>% 
  filter(anc=="afr") %>%
  mutate(fdr = p.adjust(p,method="fdr"))
)


overall_res_table_additive$direction = ifelse(overall_res_table_additive$beta>0,"up","down")
p =  ggplot(overall_res_table_additive,aes(full_allele,-log10(p),fill=gene,shape=direction))+
    geom_point(data = overall_res_table_additive %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table_additive %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="HLA locus",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table_additive %>% filter(p<0.05),mapping = aes(label = full_allele))+
    theme(axis.text.x = element_blank())+
    facet_wrap(~toupper(anc),nrow=2)

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all_TWODIGIT.png",res=900,units="in",width=10,height=6)
p
dev.off()
write_csv(overall_res_table_additive,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_additive_TWODIGIT.csv")



#######################
# TWO FIELD RESOLUTION
#######################

make_hla_forest = function(ancestry,allele,model="additive"){
  message(ancestry)
  message(allele)

  # nantes
  in_file = paste0("/data/scratch/hmy117/hla_allele_nantes_",allele,"_",ancestry,".tsv")
  # 1KG
  # in_file = paste0("/data/scratch/hmy117/hla_allele_kg_",allele,"_",ancestry,".tsv")

  # step 1 - compare frequencies & get global freqs
  allele_calls = read_tsv(in_file,col_types="cccdd")
  colnames(allele_calls) = c("IID","allele1","allele2","prob","match")
  allele_counts = allele_calls %>%
    dplyr::select(c(1:3)) %>%
    pivot_longer(cols = c(2,3)) %>%
    dplyr::select(-name) %>%
    group_by(IID) %>%
    dplyr::count(value)

    # join with phenot & covars

  pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv"),col_types = cols(.default="c"))
  cov = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv"),col_types = "ccdddddddddddd")

  allele_counts = allele_counts %>%
    ungroup() %>%
    left_join(pheno,by="IID") %>%
    left_join(cov,by="IID")

  # get disease-stratified AFs
  all_allele_freqs = allele_counts %>%
    filter(!is.na(MS_status)) %>%
    group_by(value,MS_status) %>%
    summarise(ac = sum(n)) %>%
    ungroup() %>%
    mutate(MS_status = ifelse(MS_status == 1, "Control","MS")) %>%
    group_by(MS_status) %>%
    mutate(total = sum(ac)) %>%
    mutate(af = ac/total) %>%
    ungroup() %>%
    pivot_wider(id_cols = value,names_from  = MS_status,values_from = c(ac,total,af)) %>%
    mutate(delta_freq = abs(af_MS - af_Control)) %>%
    mutate(log2_fold_enrichment = log2(af_MS/af_Control)) %>%
    arrange(desc(delta_freq))


  all_allele_freqs$gene = allele
  all_allele_freqs$anc = ancestry

  # add total af 
  all_allele_freqs  = all_allele_freqs %>% 
    mutate(AF = (ac_Control + ac_MS) / (total_Control + total_MS))

  # add global frequencies
  this_gene = allele
  global_afs = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/afnd.tsv.txt") %>%
    filter(grepl("USA NMDP",population)) %>%
    mutate(alleles_over_2n = as.numeric(alleles_over_2n)) %>%
    filter(!is.na(alleles_over_2n)) %>%
    filter(group == "hla" & gene == this_gene) %>%
    mutate(allele = str_remove(allele,"\\*")) %>%
    mutate(allele = str_remove(allele,paste0("^",this_gene))) %>%
    filter(allele %in% all_allele_freqs$value) %>%
    dplyr::select(allele,population,alleles_over_2n) %>%
    pivot_wider(id_cols = allele,names_from = population,values_from = alleles_over_2n)


  all_allele_freqs = all_allele_freqs %>%
      dplyr::rename("allele"= value) %>%
      left_join(global_afs,by="allele")

  common_alleles = all_allele_freqs %>% filter(AF >= 0.01)

  # regress vs MS status for each allele
  alleles = unique(common_alleles$allele)

  all_res = list()
  for(i in c(1:length(alleles))){
    message(i)
    this_allele = alleles[i]

    # filter
    this_allele_counts = allele_counts %>% filter(value == this_allele)

    
    # get people with 0 counts
    this_allele_counts_zeroes = allele_counts %>% filter(!IID %in% this_allele_counts$IID)  %>%
    distinct(IID,.keep_all=T) %>%
    mutate(n = 0)

    all_dat = this_allele_counts %>%
      bind_rows(this_allele_counts_zeroes)

    all_dat$MS_status = relevel(factor(all_dat$MS_status),ref="1")

    # fisher's exact test (dominanat coding)
    tbl = table(all_dat$MS_status,all_dat$n>=1)
    fish = fisher.test(tbl)
    fish_p = fish$p.value
    names(fish_p) = "fish_p"
    fish_or = fish$estimate 

    # dominant & recessive coding
    all_dat = all_dat %>%
        mutate(n_rec = ifelse(n <=1 ,0, 1)) %>%
        mutate(n_dom = ifelse(n >=1 ,1, 0))

    # make all models and compare

    hla_model = if(model == "additive"){
      glm(data=all_dat,factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + n, family=binomial(link="logit"))
    } else if(model=="recessive"){
      glm(data=all_dat, factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  n_rec, family=binomial(link="logit"))
    } else if(model=="dominant"){
      glm(data=all_dat, factor(MS_status) ~  sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + n_dom, family=binomial(link="logit"))
    }


    if(nrow(summary(hla_model)$coefficients)==13){
      all_res[[i]] = c(summary(hla_model)$coefficients[13,],fish_p,fish_or)
    } else {
      na_vec = c(NA,NA,NA,NA,NA,NA)
      names(na_vec) = colnames(summary(hla_model)$coefficients)
      all_res[[i]] = na_vec
    }
  }
  all_res = do.call("bind_rows",all_res)
  colnames(all_res) = c("beta","se","z","p","fisher_p","fisher_or")
  all_res$allele = alleles

  all_res = all_res %>%
    inner_join(all_allele_freqs,
    by="allele")

  overall_res_table[[length(overall_res_table)+1]] <<- all_res

  all_res = all_res %>% arrange(desc(beta))

  all_res$hla_allele = factor(all_res$allele,ordered=T,levels = all_res$allele)
  outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_allele_",allele,"_",ancestry,".png")
  p = ggplot(all_res,aes(beta,hla_allele))+
    geom_errorbarh(mapping = aes(xmin = beta - 1.96*se,xmax = beta +1.96*se,y=hla_allele),height=0.1)+
    geom_point(size=3,shape=21,fill="orange")+
    geom_vline(xintercept=0,linetype="dashed",color="pink")+
    labs(y=paste0(allele," allele"),x="Beta (log Odds Ratio)")+
    theme_bw()+
    ggtitle(paste0(toupper(ancestry),"\n",allele))


  png(outfile,res=900,units="in",width=4,height=4)
  print(p)
  dev.off()
}

# make parameter table
overall_res_table = list()
overall_allele_freqs = list()
alleles = c("A","B","C","DPB1","DQB1","DRB1")
ancestries = c("sas","afr")
param_tbl = expand.grid(alleles,ancestries)

## ADDITIVE
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i])
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",allele))

# n unique genes

overall_res_table_additive = overall_res_table %>%
  mutate(OR = exp(beta), lower_ci = exp(beta - 1.96*se),upper_ci = exp(beta+1.96*se)) 

overall_res_table_additive = 
bind_rows(
  overall_res_table_additive %>% 
  filter(anc=="sas") %>%
  mutate(fdr = p.adjust(p,method="fdr")),
overall_res_table_additive %>% 
  filter(anc=="afr") %>%
  mutate(fdr = p.adjust(p,method="fdr"))
)


overall_res_table_additive$direction = ifelse(overall_res_table_additive$beta>0,"up","down")
p =  ggplot(overall_res_table_additive,aes(full_allele,-log10(p),fill=gene,shape=direction))+
    ggrepel::geom_text_repel(data = overall_res_table_additive %>% filter(p<0.05),mapping = aes(label = full_allele),size=3,min.segment.length=0)+
    geom_point(data = overall_res_table_additive %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table_additive %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="HLA locus",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    theme(axis.text.x = element_blank())+
    facet_wrap(~toupper(anc),nrow=2)

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all.png",res=900,units="in",width=12,height=4)
p
dev.off()


write_csv(overall_res_table_additive,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_additive.csv")

# more formatting 
## get significant associations
overall_res_table_additive_simple = overall_res_table_additive %>% 
  dplyr::select(anc,full_allele,contains("af_"),p,fdr,fisher_p,fisher_or,OR,lower_ci,upper_ci,`USA NMDP African`,contains("South Asian"),contains("European")) %>% 
  mutate(sas_enriched = ifelse(
    `USA NMDP South Asian Indian` > 0.05 & `USA NMDP European Caucasian` < 0.01,
     "enriched",
     " "
  )) %>% 
  mutate(afr_enriched = ifelse(
    `USA NMDP African` > 0.05  & `USA NMDP European Caucasian` < 0.01,
     "enriched",
     " "
  )) %>% 
  mutate(or_ci = paste0(round(OR,2)," (",round(lower_ci,1)," - ",round(upper_ci,1),")")) %>% 
  pivot_wider(id_cols = c(full_allele,contains("NMDP"),contains("enriched")),
  names_from = anc,
  values_from = c(contains("af"),p,fdr,or_ci,fisher_p,fisher_or))

  write_csv(overall_res_table_additive_simple,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_additive_simple.csv")





## RECESSIVE
overall_res_table = list()
overall_allele_freqs = list()
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i],model="recessive")
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",allele))

# n unique genes

overall_res_table_additive = overall_res_table %>%
  mutate(OR = exp(beta), lower_ci = exp(beta - 1.96*se),upper_ci = exp(beta+1.96*se)) 

overall_res_table_additive = 
bind_rows(
  overall_res_table_additive %>% 
  filter(anc=="sas") %>%
  mutate(fdr = p.adjust(p,method="fdr")),
overall_res_table_additive %>% 
  filter(anc=="afr") %>%
  mutate(fdr = p.adjust(p,method="fdr"))
)


overall_res_table_additive$direction = ifelse(overall_res_table_additive$beta>0,"up","down")
p =  ggplot(overall_res_table_additive,aes(full_allele,-log10(p),fill=gene,shape=direction))+
    ggrepel::geom_text_repel(data = overall_res_table_additive %>% filter(p<0.05),mapping = aes(label = full_allele),size=3,min.segment.length=0)+
    geom_point(data = overall_res_table_additive %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table_additive %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="HLA locus",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    theme(axis.text.x = element_blank())+
    facet_wrap(~toupper(anc),nrow=2)

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_recessive.png",res=900,units="in",width=12,height=4)
p
dev.off()


write_csv(overall_res_table_additive,"/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_recessive.csv")






## DOMINANT

overall_res_table = list()
overall_allele_freqs = list()
# loop through and do associations
for(i in c(1:nrow(param_tbl))){
  make_hla_forest(param_tbl$Var2[i],param_tbl$Var1[i],model="dominant")
}

# combine all results
overall_res_table = do.call("bind_rows",overall_res_table)

# make full allele name
overall_res_table = overall_res_table %>%
  mutate(full_allele = paste0(gene,"*",allele))

bonf = 0.05/6
overall_res_table$direction = ifelse(overall_res_table$beta>0,"up","down")
p =  ggplot(overall_res_table,aes(full_allele,-log10(p),fill=toupper(anc),shape=direction))+
    geom_point(data = overall_res_table %>% filter(direction=="up"),size=2,color="black",shape=24)+
    geom_point(data = overall_res_table %>% filter(direction=="down"),size=2,color="black",shape=25,show.legend=F)+
    scale_fill_brewer(palette="Set1")+
    geom_hline(yintercept=-log10(bonf),linetype="dashed",color="pink")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="pink")+
    labs(fill="Ancestry",x="HLA allele",y=bquote(-log[10]~P))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggrepel::geom_text_repel(data = overall_res_table %>% filter(p<0.05),mapping = aes(label = full_allele))

png("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_plot_alleles_all_dominant.png",res=900,units="in",width=12,height=6)
p
dev.off()



##############################
# PAF 
##############################

library(tidyverse)
hla_dat1 = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/moutsianas_alleles.csv") %>% 
  mutate(sas_paf = (sas_case_af^2 + 2*sas_case_af*(1-sas_case_af) ) * (1 - 1 / sas_or)  ) %>% 
  mutate(afr_paf = (afr_case_af^2 + 2*afr_case_af*(1-afr_case_af) ) * (1 - 1 / afr_or)  ) %>% 
  mutate(eur_paf = (eur_case_af^2 + 2*eur_case_af*(1-eur_case_af) ) * (1 - 1 / eur_or)  ) %>% 
  pivot_longer(contains("paf")) %>% 
  dplyr::select(hla_allele,name,value) %>% 
  mutate(name = str_remove_all(name,"_paf")) %>% 
  mutate(paf = 100* value) %>% 
  mutate(PAF = "Miettinen's")
  

hla_dat2 = read_csv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/moutsianas_alleles.csv") %>% 
  mutate(sas_paf = ( (sas_cont_af^2 + 2*sas_cont_af*(1-sas_cont_af) ) * (sas_or - 1) ) / (1 + (sas_cont_af^2 + 2*sas_cont_af*(1-sas_cont_af) ) * (sas_or - 1) ) ) %>% 
  mutate(afr_paf = ( (afr_cont_af^2 + 2*afr_cont_af*(1-afr_cont_af) ) * (afr_or - 1) ) / (1 + (afr_cont_af^2 + 2*afr_cont_af*(1-afr_cont_af) ) * (afr_or - 1) ) ) %>% 
  mutate(eur_paf = ( (eur_cont_af^2 + 2*eur_cont_af*(1-eur_cont_af) ) * (eur_or - 1) ) / (1 + (eur_cont_af^2 + 2*eur_cont_af*(1-eur_cont_af) ) * (eur_or - 1) ) ) %>% 
  pivot_longer(contains("paf")) %>% 
  dplyr::select(hla_allele,name,value) %>% 
  mutate(name = str_remove_all(name,"_paf")) %>% 
  mutate(paf = 100* value) %>% 
  mutate(PAF = "Levin's")
  
dat =  bind_rows(hla_dat1,hla_dat2)

  
p = ggplot(dat,aes(
  hla_allele,paf,fill=toupper(name)
)) +
geom_col(color="black",position=position_dodge())+
theme_bw()+
labs(x="HLA allele",y="PAF (%)",fill="Ancestry")+
coord_flip()+
scale_fill_brewer(palette="Set1")+
facet_wrap(~PAF,nrow=2)


# prep pheno & covar 
ancestry = "sas"

pheno = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_pheno.tsv"),col_types = cols(.default="c"))
cov = read_tsv(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/pheno/reimputed_",ancestry,"_covars_with_pcs.tsv"),col_types = "ccdddddddddddd")
pheno = pheno %>% 
  left_join(cov,by="IID")

# step 1 - read in all alleles and get a wide-format dataset 
alleles = c("A","B","C","DPB1","DQB1","DRB1")

for(allele in alleles){

  # nantes
  in_file = paste0("/data/scratch/hmy117/hla_allele_nantes_",allele,"_",ancestry,".tsv")
  # 1KG
  # in_file = paste0("/data/scratch/hmy117/hla_allele_kg_",allele,"_",ancestry,".tsv")

  # step 1 - compare frequencies & get global freqs
  allele_calls = read_tsv(in_file,col_types="cccdd")
  colnames(allele_calls) = c("IID","allele1","allele2","prob","match")

  allele_counts = allele_calls %>%
    dplyr::select(c(1:3)) %>%
    pivot_longer(cols = c(2,3)) %>%
    dplyr::select(-name) %>%
    group_by(IID) %>%
    dplyr::count(value)

  # cast wide 
  allele_counts = allele_counts %>%
    ungroup() %>% 
    mutate(allele = paste0(allele,"*",value)) %>% 
    pivot_wider(id_cols = IID,values_from = n, names_from = allele)
  wide_allele_counts = allele_counts %>% replace(is.na(.),0)

  # add to phenotype 
  pheno <<- pheno %>% 
  left_join(wide_allele_counts,by="IID")
}

# do regression 1 
## get col names 
alleles_to_regress = colnames(pheno)[17:ncol(pheno)]
model_dat = pheno %>% 
  mutate(MS_status = ifelse(MS_status == 2,1,0))

# loop through all alleles 
res = list()
for(i in c(1:length(alleles_to_regress))){
  message(i)
  model = glm(data = model_dat, 
  MS_status ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
  PC8 + PC9 + PC10 + model_dat[[alleles_to_regress[i]]], family=binomial(link="logit"))

  coefs = broom::tidy(model)
  res[[i]] = c("allele"=alleles_to_regress[i],"p"=coefs[13,]$p.value)

}
res = do.call("bind_rows",res)
res$p = as.numeric(res$p)
res = res %>% arrange(p)
tophit = res %>% slice_min(p)

# step 2
res = list()
for(i in c(1:length(alleles_to_regress))){
  message(i)
  model = glm(data = model_dat, 
  MS_status ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
  PC8 + PC9 + PC10 + model_dat[[alleles_to_regress[i]]] + model_dat[[tophit$allele[1]]], family=binomial(link="logit"))

  coefs = broom::tidy(model)
  res[[i]] = c("allele"=alleles_to_regress[i],"p"=coefs[13,]$p.value)

}
res = do.call("bind_rows",res)
res$p = as.numeric(res$p)
res = res %>% filter(!allele %in% tophit$allele)
res = res %>% arrange(p)



