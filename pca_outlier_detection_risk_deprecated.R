library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

args = commandArgs(trailingOnly=T)

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

# find outliers
ancestry = args[1]
out_ancestry=args[2]
sd_num=3
n_controls = 10

# read in pcs
infile = paste0("/data/scratch/hmy117/pcs_vec_",out_ancestry)
pcs = read_tsv(infile,col_types = cols(.default="c"))


# combine with covars & standardize pcs
z_score = function(x){
	( as.numeric(x) - mean( as.numeric(x)) ) / sd(as.numeric(x))
}
pcs = pcs %>%
left_join(all_cov,by="IID") %>%
	filter(!is.na(cohort))	 %>%
	mutate_at(.vars = vars(contains("PC")),z_score)

# fx for eulidean distances
euclidean_dist = function(x,y){
	sqrt(rowSums((x - y)^2 ))
}


# define cases & controls
cases = if(ancestry=="EUR"){
	pcs %>% filter(MS_status=="MS" & cohort != "CAM")
	} else {
	pcs %>% filter(MS_status=="MS")
	}

controls = pcs %>% filter(MS_status=="Control")

matched_controls = list()
# loop through each case and find 20 nearest neighbours
for(i in c(1:nrow(cases))){
this_case = cases[i,]
message("matching case", i)
matched_controls_df = do.call("bind_rows",matched_controls)
message("There are now ", nrow(matched_controls_df)," matched controls")

these_matched_controls = controls %>%
  filter(!IID %in% matched_controls_df$IID)

# make numeric matrix
these_matched_controls_mat = as.matrix(these_matched_controls %>%
dplyr::select(PC1,PC2,PC3,PC4))
this_case_mat = as.matrix(
	data.frame(PC1 = rep(this_case$PC1,nrow(these_matched_controls)),
	PC2 = rep(this_case$PC2,nrow(these_matched_controls)),
	PC3 = rep(this_case$PC3,nrow(these_matched_controls)),
	PC4 = rep(this_case$PC4,nrow(these_matched_controls))
	))

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

all_dat = bind_rows(cases,matched_controls_df)

# plot UKB ethnicities
plot_dat = all_dat %>%
        left_join(  ukb_pheno_ethnicity %>%
      mutate("IID" = as.character(EID)),
      by="IID")
top_ethnicities = plot_dat %>%
  dplyr::count(ethnic_background_f21000_0_0) %>%
  filter(!is.na(ethnic_background_f21000_0_0)) %>%
  slice_max(n,n=8)
plot_dat  = plot_dat %>%
      mutate(ethnicity = ifelse(
			ethnic_background_f21000_0_0 %in% top_ethnicities$ethnic_background_f21000_0_0,
			as.character(ethnic_background_f21000_0_0),
			"Other/missing"))

p1 = ggplot(data = plot_dat,
      aes(PC1,PC2,fill=ethnicity))+
      facet_wrap(~MS_status,nrow=3)+
      geom_point(shape=21,color="black")+
      theme_bw()+
      labs(fill="Ethnicity")+
      scale_fill_brewer(palette="Set1")

png(paste0("./outputs/ukb_cam_adams_pca_ukb_ethnicity",ancestry,".png"),res=900,units="in",width=8,height=8)
print(p1)
dev.off()


# first define outliers using euclidean distance from population mean
means = all_dat %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	summarise_all(mean)

mean_mat = as.matrix(data.frame(PC1 = rep(means$PC1,nrow(all_dat)),
PC2 = rep(means$PC2,nrow(all_dat)),
PC3 = rep(means$PC3,nrow(all_dat)),
PC4 = rep(means$PC4,nrow(all_dat))
))

# calculate euclidean distance
all_dat_mat = all_dat %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	as.matrix()

all_dat$euclidean_distance =  euclidean_dist(all_dat_mat,mean_mat)

threshold = 3
outlier_ids = all_dat %>%
	mutate(outlier = ifelse(euclidean_distance > threshold,"yes","no")) %>%
	filter(outlier=="yes")

all_dat = all_dat %>%
mutate(outlier = ifelse(IID %in% outlier_ids$IID,
  "outlier",
  "keep"))
message("PCA outliers:")
print(all_dat %>% dplyr::count(outlier,cohort) %>% mutate(prop = n/sum(n)))

p = ggplot(all_dat,aes(PC1,PC2,fill=euclidean_distance))+
	geom_point(data = all_dat %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = all_dat %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")
p1 = ggplot(all_dat,aes(PC3,PC4,fill=euclidean_distance))+
	geom_point(data = all_dat %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = all_dat %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")

png(paste0("./outputs/ukb_cam_adams_pca_outliers_",ancestry,".png"),res=900,units="in",width=8,height=3)
print(cowplot::plot_grid(p,p1,align="h",ncol=2))
dev.off()

# save non-outliers to file
all_dat = all_dat %>% filter(outlier=="keep")
non_outliers = all_dat %>% dplyr::select(FID,IID)
write_tsv(non_outliers,paste0("./outputs/pca_non_outliers_susceptibility",ancestry,".tsv"))

p1 = ggplot(data = all_dat,
aes(PC1,PC2,fill=cohort))+
geom_point(shape=21,color="black")+
theme_bw()
p2 = ggplot(all_dat,
aes(PC1,PC2,fill=MS_status))+
geom_point(shape=21,color="black")+
theme_bw()+
labs(fill="MS status")

png(paste0("./outputs/ukb_cam_adams_pca_no_outliers_by_ms_status_and_cohort",ancestry,".png"),res=900,units="in",width=8,height=4)
print(cowplot::plot_grid(p1,p2,align="h",ncol=2))
dev.off()

# re-run pca with these people
# run PCA
cmd = paste0("~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_pca_",out_ancestry," ",
"--pca approx 4 ",
"--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pca_non_outliers_susceptibility",ancestry,".tsv ",
"--out /data/scratch/hmy117/round_2_pcs_",out_ancestry)

system(cmd)

# read in refreshed pcs
pcs_2 = read_table(paste0("/data/scratch/hmy117/round_2_pcs_",out_ancestry,".eigenvec"),col_types = cols(.default = "c")) %>%
mutate_at(.vars = vars(contains("PC")),z_score) %>%
left_join(all_cov %>% dplyr::select(-`#FID`),by="IID")


# first define outliers using euclidean distance from population mean
means = pcs_2 %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	summarise_all(mean)

mean_mat = as.matrix(data.frame(PC1 = rep(means$PC1,nrow(pcs_2)),
PC2 = rep(means$PC2,nrow(pcs_2)),
PC3 = rep(means$PC3,nrow(pcs_2)),
PC4 = rep(means$PC4,nrow(pcs_2))
))

# calculate euclidean distance
all_dat_mat = pcs_2 %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	as.matrix()

pcs_2$euclidean_distance =  euclidean_dist(all_dat_mat,mean_mat)

threshold = 3
outlier_ids = pcs_2 %>%
	mutate(outlier = ifelse(euclidean_distance > threshold,"yes","no")) %>%
	filter(outlier=="yes")

pcs_2 = pcs_2 %>%
mutate(outlier = ifelse(IID %in% outlier_ids$IID,
  "outlier",
  "keep"))
message("PCA outliers:")
print(pcs_2 %>% dplyr::count(outlier,cohort) %>% mutate(prop = n/sum(n)))


p = ggplot(pcs_2,aes(PC1,PC2,fill=euclidean_distance))+
	geom_point(data = pcs_2 %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = pcs_2 %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")
p1 = ggplot(pcs_2,aes(PC3,PC4,fill=euclidean_distance))+
	geom_point(data = pcs_2 %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = pcs_2 %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")

png(paste0("./outputs/ukb_cam_adams_pca_outliers_round_2",ancestry,".png"),res=900,units="in",width=8,height=3)
print(cowplot::plot_grid(p,p1,align="h",ncol=2))
dev.off()

# filter to non-outliers
pcs_2 = pcs_2 %>% filter(outlier == "keep")
write_tsv(pcs_2 %>% dplyr::select(1,2),
paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pca_non_outliers_susceptibility_round_2",ancestry,".tsv"))

# final round of re-computing PCs
# re-run pca with these people
# run PCA
cmd = paste0("~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_pca_",out_ancestry," ",
"--pca 4 approx ",
"--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pca_non_outliers_susceptibility_round_2",ancestry,".tsv ",
"--out /data/scratch/hmy117/round_3_pcs_",out_ancestry)

system(cmd)

# read in refreshed pcs

# read in refreshed pcs
pcs_3 = read_table(paste0("/data/scratch/hmy117/round_3_pcs_",out_ancestry,".eigenvec"),col_types = cols(.default = "c")) %>%
mutate_at(.vars = vars(contains("PC")),z_score) %>%
left_join(all_cov %>% dplyr::select(-`#FID`),by="IID")


# first define outliers using euclidean distance from population mean
means = pcs_3 %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	summarise_all(mean)

mean_mat = as.matrix(data.frame(PC1 = rep(means$PC1,nrow(pcs_3)),
PC2 = rep(means$PC2,nrow(pcs_3)),
PC3 = rep(means$PC3,nrow(pcs_3)),
PC4 = rep(means$PC4,nrow(pcs_3))
))

# calculate euclidean distance
all_dat_mat = pcs_3 %>%
	dplyr::select(PC1,PC2,PC3,PC4) %>%
	as.matrix()

pcs_3$euclidean_distance =  euclidean_dist(all_dat_mat,mean_mat)

threshold = 3
outlier_ids = pcs_3 %>%
	mutate(outlier = ifelse(euclidean_distance > threshold,"yes","no")) %>%
	filter(outlier=="yes")

pcs_3 = pcs_3 %>%
mutate(outlier = ifelse(IID %in% outlier_ids$IID,
  "outlier",
  "keep"))
message("PCA outliers:")
print(pcs_3 %>% dplyr::count(outlier,cohort) %>% mutate(prop = n/sum(n)))


p = ggplot(pcs_3,aes(PC1,PC2,fill=euclidean_distance))+
	geom_point(data = pcs_3 %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = pcs_3 %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")
p1 = ggplot(pcs_3,aes(PC3,PC4,fill=euclidean_distance))+
	geom_point(data = pcs_3 %>% filter(outlier =="outlier"),shape=21,color="black",alpha=0.3)+
	geom_point(data = pcs_3 %>% filter(outlier =="keep"),shape=21,color="black")+
	scale_fill_gradient(low = "purple",high="orange")+
	theme_bw()+
	labs(fill="Euclidean distance\nfrom mean")

png(paste0("./outputs/ukb_cam_adams_pca_outliers_round_3",ancestry,".png"),res=900,units="in",width=8,height=3)
print(cowplot::plot_grid(p,p1,align="h",ncol=2))
dev.off()

# filter to non-outliers
pcs_3 = pcs_3 %>% filter(outlier == "keep")
write_tsv(pcs_3 %>% dplyr::select(1,2),
paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pca_non_outliers_susceptibility_round_3",ancestry,".tsv"))

# final round of re-computing PCs
# re-run pca with these people
# run PCA
cmd = paste0("~/plink2 --bfile /data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_pca_",out_ancestry," ",
"--pca 10 approx ",
"--keep /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/pca_non_outliers_susceptibility_round_3",ancestry,".tsv ",
"--out /data/scratch/hmy117/round_4_pcs_",out_ancestry)

system(cmd)

# read in refreshed pcs
pcs_4 = read_table(paste0("/data/scratch/hmy117/round_4_pcs_",out_ancestry,".eigenvec"),col_types = cols(.default = "c")) %>%
mutate_at(.vars = vars(contains("PC")),as.numeric) %>%
left_join(all_cov %>% dplyr::select(-`#FID`),by="IID")


# plot
p1 = ggplot(data = pcs_4,
  aes(PC1,PC2,fill=cohort))+
  geom_point(shape=21,color="black")+
  theme_bw()
p2 = ggplot(pcs_4,
  aes(PC1,PC2,fill=MS_status))+
  geom_point(shape=21,color="black")+
  theme_bw()+
  labs(fill="MS status")

png(paste0("./outputs/ukb_cam_adams_pca_without_outliers_",ancestry,".png"),res=900,units="in",width=8,height=4)
print(cowplot::plot_grid(p1,p2,align="h",ncol=2))
dev.off()


# write covar file
all_cov = pcs_4 %>%
filter(!is.na(PC1)) %>%
dplyr::rename("FID" = `#FID`) %>%
dplyr::select(FID,IID,age,sex,contains("PC")) %>%
na.omit()

write_tsv(all_cov,paste0("./pheno/susceptibility_",out_ancestry,"_covars_with_pcs.tsv"))


# make pheno file
all_pheno = pcs_4 %>%
filter(!is.na(PC1)) %>%
dplyr::rename("FID" = `#FID`) %>%
dplyr::select(FID,IID,MS_status) %>%
na.omit() %>%
mutate(MS_status = ifelse(MS_status == "Control",1,2))

write_tsv(all_pheno,paste0("./pheno/susceptibility_",out_ancestry,"_pheno.tsv"))

# plot pcs vs cohort
pcs_long = pcs_4 %>% pivot_longer(cols = contains("PC"))
pcs_long$PC = as.numeric(str_remove_all(pcs_long$name,"PC"))

pvals = list()
for(i in c(1:10)){
pc_colname = paste0("PC",i)
pvals[[i]] = summary(lm(data = pcs_4,pcs_4[[pc_colname]] ~ age + sex + MS_status))$coefficients[4,4]
}
pvals_df = data.frame(p = unlist(pvals),PC = c(1:10))
pvals_df = pvals_df %>% mutate(p_ind = ifelse(p<0.05,"*","NS"))

p = ggplot()+
geom_boxplot(data = pcs_long,aes(factor(PC),value,fill=MS_status),outlier.shape=NA)+
theme_bw()+
scale_fill_brewer(palette="Set1")+
labs(x="PC",y="Eigenvector")+
scale_y_continuous(limits = c(-0.25,0.25))+
geom_text(data = pvals_df,aes(x = PC,y = 0.2,label = p_ind))

png(paste0("./outputs/ukb_cam_adams_pca_boxplots_vs_ms_status",ancestry,".png"),res=900,units="in",width=10,height=4)
print(p)
dev.off()
