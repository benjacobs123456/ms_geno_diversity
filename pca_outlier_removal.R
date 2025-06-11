library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
# args = c("CSA","sas")
# args = c("AFR","afr")
# args = c("CSA","sas")
set.seed(123456)


# read in combined pheno-covar data & ancestry calls
pheno_cov_data = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/merged_ukb_adams_pheno_cov.tsv")
anc_calls = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_all_ukb_adams.tsv")

# combine
all_cov = anc_calls %>% left_join(pheno_cov_data,by="IID")

ukb_pheno_ethnicity = readRDS("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/78867r672482/78867r672482_FO.rds") %>%
	tibble %>%
	dplyr::select(eid,contains("ethnic"))
id_bridge = read_csv("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/59138_78867/Bridge_eids_59138_78867.csv")

ukb_pheno_ethnicity = ukb_pheno_ethnicity %>%
  dplyr::rename("eid_78867" = eid) %>%
  left_join(id_bridge,by="eid_78867") %>%
  dplyr::rename("EID" = eid_59138)
ukb_pheno_ethnicity = ukb_pheno_ethnicity %>% dplyr::select(EID,contains("ethnic"))

# set parameters
ancestry = args[1]
out_ancestry=args[2]

# filter covars to this ancestry 
all_cov = all_cov %>% 
    filter(predicted_ancestry == ancestry)

# add in within-ancestry PCs (vs PCs computed on whole cohort )
pcs = read_table(paste0("/data/scratch/hmy117/pcs_",out_ancestry))

all_cov = all_cov %>% 
    dplyr::select(-contains("PC"),-contains("FID")) %>% 
    left_join(pcs,by="IID")

# plot all PCs 
png(paste0("./outputs/susceptibility_ALL_",out_ancestry,"_pc_pairs_plot.png"),res=900,units="in",width=12,height=12)
GGally::ggpairs(
    all_cov,
    columns = paste0("PC",c(1:10)),
    aes(colour = cohort))+
    scale_color_brewer(palette="Set1")+
    theme_bw()
dev.off()

# plot cohort vs PCs 
pcs_long = all_cov %>% 
    pivot_longer(cols = c(PC1:PC10)) %>% 
    mutate(PC = as.numeric(str_remove_all(name,"PC"))) %>% 
    mutate(Eigenvector = as.numeric(value))

eigenvals = read_table(paste0("/data/scratch/hmy117/eigenvals_",args[2]),col_names=F) %>% 
    mutate(PC = c(1:10),Eigenvalue = X1)

# check association of PCs with cohort of origin
p_list = list()
for(i in c(1:10)){
    col = paste0("PC",i)

    model = glm(
        data = all_cov,
        factor(cohort) ~ all_cov[[col]],
        family=binomial(link = "logit")
    )
    p_list[[i]] = summary(model)$coefficients[2,4]

}
p_dat = data.frame( 
  PC = c(1:10),
  p = unlist(p_list)
) %>% mutate(sig = case_when(
    p < 0.05/10000 ~ "****",
    p < 0.05/1000 ~ "***",
    p < 0.05/100 ~ "**",
    p < 0.05/10 ~ "*"
))
    
    

p = ggplot(pcs_long,aes(factor(PC),Eigenvector,fill=cohort))+
    geom_boxplot()+
    theme_bw()+
    scale_fill_brewer(palette="Set1")+
    labs(x="PC",y="Eigenvector",fill="Cohort") + 
    theme(legend.position = "top")+
    geom_text(data = p_dat, mapping = aes(x = PC,y = 0.5, label = sig, fill = NULL))

p2 = ggplot(eigenvals,aes(factor(PC),Eigenvalue))+
    geom_point(size = 3)+
    theme_bw()+
    labs(x="PC",y="Eigenvalue")

png(paste0("./outputs/susceptibility_ALL_",out_ancestry,"_pc_eigenvals.png"),res=900,units="in",width=8,height=8)
cowplot::plot_grid(p,p2,ncol=1,align="v")
dev.off()


# standardize pcs
z_score = function(x){
	( as.numeric(x) - mean( as.numeric(x)) ) / sd(as.numeric(x))
}
all_cov = all_cov %>%
	filter(!is.na(cohort))	 %>%
    filter(!is.na(PC1)) %>%
	mutate_at(.vars = vars(contains("PC")),z_score)

# fx for eulidean distances
euclidean_dist = function(x,y){
	sqrt(rowSums((x - y)^2 ))
}


# define cases & controls
cases = all_cov %>% filter(MS_status=="MS")
controls = all_cov %>% filter(MS_status=="Control")

# initialise 
pre_filtering_dat = all_cov 

# exclude gross outliers in PCs 1 - 10
for(i in c(1:10)){
    this_pc = paste0("PC",i)
    lower = mean(cases[[this_pc]]) - 5*sd( cases[[this_pc]])
    upper = mean(cases[[this_pc]]) + 5*sd( cases[[this_pc]])
    pre = nrow(all_cov %>% filter(MS_status=="MS"))
    pre_cont = nrow(all_cov %>% filter(MS_status!="MS"))

    all_cov <<- all_cov %>% filter(.data[[this_pc]] > lower & .data[[this_pc]] < upper)
    post = nrow(all_cov %>% filter(MS_status=="MS"))
    post_cont = nrow(all_cov %>% filter(MS_status!="MS"))

    message("Removed ",pre-post," cases")
    message("Removed ",pre_cont-post_cont," controls")
    
}

# plot  
pre_filtering_dat = pre_filtering_dat %>% 
    mutate(outlier = ifelse(!IID %in% all_cov$IID,"outlier","non_outlier"))
p = ggplot(pre_filtering_dat,aes(PC1,PC2,col=cohort))+
    geom_point(data = pre_filtering_dat %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = pre_filtering_dat %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(legend.position="none")
p2 = ggplot(pre_filtering_dat,aes(PC3,PC4,col=cohort))+
    geom_point(data = pre_filtering_dat %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = pre_filtering_dat %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(legend.position="none")
p3 = ggplot(pre_filtering_dat,aes(PC5,PC6,col=cohort))+
    geom_point(data = pre_filtering_dat %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = pre_filtering_dat %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(legend.position="none")
p4 = ggplot(pre_filtering_dat,aes(PC7,PC8,col=cohort))+
    geom_point(data = pre_filtering_dat %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = pre_filtering_dat %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(legend.position="none")
p5 = ggplot(pre_filtering_dat,aes(PC9,PC10,col=cohort))+
    geom_point(data = pre_filtering_dat %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = pre_filtering_dat %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(legend.position="none")

# save plot 
png(paste0("./outputs/susceptibility_ALL_",out_ancestry,"_pc_plot.png"),res=900,units="in",width=10,height=4)
print(cowplot::plot_grid(p,p2,p3,p4,p5,ncol=5,align="h"))
dev.off()

# impute missing ages 

prop_f = sum(all_cov$sex==2,na.rm=T) / ( sum(all_cov$sex==2,na.rm=T) + sum(all_cov$sex==1,na.rm=T) ) 
covars_all = all_cov %>%
    filter(!is.na(PC1)) %>%
    dplyr::select(FID,IID,age,sex,contains("PC")) %>% 
        mutate(age = ifelse(is.na(age),median(age,na.rm=T),age)) %>% 
        mutate(sex = ifelse(is.na(sex),rbinom(size = 1,n=1,prob=prop_f) + 1,sex)) 

write_tsv(covars_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_covars_with_pcs.tsv"))

# make pheno file
pheno_all = all_cov %>%
filter(!is.na(PC1)) %>%
dplyr::select(FID,IID,MS_status) %>%
na.omit() %>%
mutate(MS_status = ifelse(MS_status == "Control",1,2))
write_tsv(pheno_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_pheno.tsv"))

