library(tidyverse)
library(RNOmni)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

ancestry_calls_hgdp_1kg = read_tsv("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ancestry_calls_detailed.tsv") %>% dplyr::select(IID,contains("ancestry"))
cov = read_tsv("./pheno/adams_covars_fixed.tsv") %>% dplyr::select(-FID)
all_cov = ancestry_calls_hgdp_1kg %>% left_join(cov,by="IID")
pheno = read_tsv("./pheno/adams_pheno.tsv") %>% dplyr::select(-FID)

# count ancestries of people with a gARMSS
pheno %>%
	filter(!is.na(gARMSS)) %>%
	left_join(ancestry_calls_hgdp_1kg,by="IID") %>%
	dplyr::count(predicted_ancestry)

# loop through each ancestry to find outliers & process phenotype data
# find outliers
find_outliers_process_pheno = function(ancestry,sd_num=2){


	pcs = read_table(paste0("./outputs/",ancestry,"_pcs.eigenvec")) %>% dplyr::select(-1)

	# define limits
	limits = pcs %>%
		dplyr::select(-1) %>%
		summarise_all(.funs = c("mean","sd")) %>%
		pivot_longer(cols = everything()) %>%
		separate(name,sep="_",into=c("PC","fx")) %>%
		pivot_wider(id_cols = PC,values_from = value,names_from = fx) %>%
		mutate(lower = mean - sd_num * sd, upper = mean + sd_num * sd)

	outlier_ids = 	pcs %>%
			pivot_longer(cols = c(2:5)) %>%
			dplyr::rename("PC" = name) %>%
			left_join(limits,by="PC") %>%
			mutate(outlier = ifelse(value < lower | value > upper,"yes","no")) %>%
			filter(outlier=="yes")

	pcs = pcs %>%
		mutate(outlier = ifelse(IID %in% outlier_ids$IID,
			"outlier",
			"keep"))
	message("PCA outliers:")
	print(pcs %>% dplyr::count(outlier) %>% mutate(prop = n/sum(n)))
	# save outliers to file
	outliers = pcs %>% filter(outlier=="outlier") %>% dplyr::select(IID) %>% mutate(FID = IID)
	write_tsv(outliers,paste0("./outputs/pca_outliers_",ancestry,".tsv"))

	# combine with main covar file
	pcs = pcs  %>% left_join(all_cov,by="IID")
	p = ggplot(pcs %>% filter(outlier=="keep"),aes(PC1,PC2,fill = bigsnpr_ancestry))+
		geom_point(size=3,data = pcs %>% filter(outlier!="keep"),alpha=0.5,shape=13,show.legend=F)+
		geom_point(size=3,color="black",shape=21,alpha=0.8)+
		theme_bw()+
		scale_fill_brewer(palette="Paired")+
		labs(fill="Ancestry group")+
		ggtitle(toupper(ancestry))
	all_plots[[length(all_plots)+1]] <<- p
	png(paste0("./outputs/pca_outliers_",ancestry,".png"),res=900,units="in",width=6,height=4)
	print(p)
	dev.off()

	# write covar file
	pcs = pcs %>% filter(!is.na(PC1)) %>% dplyr::select(-contains("ancestry")) %>% mutate(FID = IID) %>% dplyr::select(FID,IID,contains("age"),contains("sex"),contains("PC"),Site) %>% na.omit()
	write_tsv(pcs,paste0("./pheno/",ancestry,"_covars_with_pcs.tsv"))

	# make pheno file
	pheno = pheno %>%
		filter(IID %in% pcs$IID)

	# join with pcs
	pheno = pheno %>%
		mutate(FID = IID) %>%
		dplyr::select(FID,IID,gARMSS,edss,contains("msis"),eq5d_vas,age_at_dx)

 	# define edss >6
	pheno = pheno %>%
		mutate(edss6 = ifelse(edss >=6,"1","0"))

  write_tsv(pheno,paste0("./pheno/",ancestry,"_pheno.tsv"))

}

all_plots = list()
find_outliers_process_pheno(ancestry="sas",sd_num=3)
find_outliers_process_pheno(ancestry="afr",sd_num=3)
find_outliers_process_pheno(ancestry="eur",sd_num=2)

png("./outputs/pca_outliers_all_ancestries.png",res=900,units="in",width=6,height=8)
cowplot::plot_grid(plotlist=all_plots,align="v",nrow=3)
dev.off()

# write combined fixed covars file

covars = purrr::map(c("sas","afr","eur"),function(x){
	read_tsv(paste0("./pheno/",x,"_covars_with_pcs.tsv"))
	})
covars = do.call("bind_rows",covars)
write_tsv(covars,"./outputs/covars_for_prs.tsv")

# normalise phenos for PRS

# for each phenotype, apply RINT (for PRS) but add NAs back in after
for(ancestry in c("sas","afr","eur")){
pheno_for_prs = read_tsv(paste0("./pheno/",ancestry,"_pheno.tsv"))
message("N cols", ncol(pheno_for_prs))

for(pheno_to_normalise in c("gARMSS","edss","msis_physical_normalised","eq5d_vas","age_at_dx")){
pheno_no_nas = pheno_for_prs %>%
    filter(!is.na(.data[[pheno_to_normalise]])) %>%
    mutate(normalised_pheno = RNOmni::RankNorm(.data[[pheno_to_normalise]]))
pheno_nas = pheno_for_prs %>%
    filter(is.na(.data[[pheno_to_normalise]])) %>%
    mutate(normalised_pheno = NA)
pheno_all = bind_rows(pheno_no_nas,pheno_nas)
pheno_all[[paste0("rint_",pheno_to_normalise)]] = pheno_all$normalised_pheno

# save to main data frame
pheno_for_prs <<- pheno_all
}

message("N cols", ncol(pheno_for_prs))
write_tsv(pheno_for_prs,paste0("./pheno/",ancestry,"_prs_pheno.tsv"))
}
