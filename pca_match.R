library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)
# args = c("CSA","sas")
# args = c("AFR","afr")

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

# standardize pcs
z_score = function(x){
	( as.numeric(x) - mean( as.numeric(x)) ) / sd(as.numeric(x))
}
all_cov = all_cov %>%
	filter(!is.na(cohort))	 %>%
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

# exclude gross outliers in PCs 1 - 4
for(i in c(1:4)){
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
    scale_color_brewer(palette="Set1")


# loop through each case and find x nearest neighbours
matched_controls = list()
n_controls = 20
for(i in c(1:nrow(cases))){
this_case = cases[i,]
message("matching case", i)
matched_controls_df = do.call("bind_rows",matched_controls)
message("There are now ", nrow(matched_controls_df)," matched controls")

these_matched_controls = controls %>%
  filter(!IID %in% matched_controls_df$IID)

# make numeric matrix
these_matched_controls_mat = as.matrix(these_matched_controls %>%
dplyr::select(PC1,PC2))
this_case_mat = as.matrix(
	data.frame(
        PC1 = rep(this_case$PC1,nrow(these_matched_controls)),
	    PC2 = rep(this_case$PC2,nrow(these_matched_controls))	))

# calculate euclidean distance
euclidean_distances = euclidean_dist(this_case_mat,these_matched_controls_mat)

# add back in
these_matched_controls$euclidean = euclidean_distances

# select top x matches
these_matched_controls = these_matched_controls %>%
  slice_min(euclidean,n=n_controls,with_ties=F) %>%
	mutate(matched_case = this_case$IID)


# combine with existing controls
matched_controls[[i]] = these_matched_controls

}
matched_controls_df = do.call("bind_rows",matched_controls)

# recombine
all_dat = bind_rows(cases,matched_controls_df)

# set outliers 
all_cov = all_cov %>% 
    mutate(outlier = ifelse(
        !(IID %in% all_dat$IID),
        "outlier",
        "not_outlier"
    ))
p = ggplot(all_cov,aes(PC1,PC2,col=cohort))+
    geom_point(data = all_cov %>% filter(outlier == "outlier"),alpha=0.05)+
    geom_point(data = all_cov %>% filter(outlier != "outlier"),alpha=1)+
    theme_bw()+
    scale_color_brewer(palette="Set1")

# save plot 
png(paste0("./outputs/susceptibility_ALL_",out_ancestry,"_pc_plot.png"),res=900,units="in",width=6,height=6)
print(p)
dev.off()

# filter 
all_cov = all_cov %>% 
    filter(outlier !="outlier")

# prepare pheno & covar file for whole ancestry
covars_all = all_cov %>%
filter(!is.na(PC1)) %>%
dplyr::select(FID,IID,age,sex,contains("PC")) %>%
na.omit()

write_tsv(covars_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_covars_with_pcs.tsv"))

# make pheno file
pheno_all = all_cov %>%
filter(!is.na(PC1)) %>%
dplyr::select(FID,IID,MS_status) %>%
na.omit() %>%
mutate(MS_status = ifelse(MS_status == "Control",1,2))
write_tsv(pheno_all,paste0("./pheno/susceptibility_ALL_",out_ancestry,"_pheno.tsv"))

