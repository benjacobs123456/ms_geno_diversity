library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")

# read summary of all prs at r2 of 0.1
all_prs = data.frame()
for(r2 in c(0.001,0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","afr","eur")){
	prs = read_table(paste0("./outputs/prs_risk_",r2,"_",ancestry,".prsice")) %>%
	mutate(anc=ancestry,clump_r2 = r2)
	all_prs <<- bind_rows(all_prs,prs)
}
}

# plot with MHC @ r2=0.1
all_prs = all_prs %>%
	mutate(P_summ = ifelse(P < 1e-3,"<0.001",round(P,3)))
p = ggplot(all_prs %>% filter(clump_r2==0.1),
  aes(factor(Threshold),Coefficient,fill=toupper(anc),label=paste0("SNPs: ",Num_SNP,"\nP: ",P_summ)))+
  geom_errorbar(mapping = aes(x = factor(Threshold),
  ymin = Coefficient - 1.96 * `Standard.Error`,
  ymax = Coefficient + 1.96 * `Standard.Error`),
  width=0.1)+
  geom_point(size=3,shape=21,color="black")+
  facet_wrap(~toupper(anc),nrow=3)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  geom_text(mapping = aes(y = Coefficient + 1.96 * `Standard.Error` + 0.1),size=3)+
  geom_hline(yintercept=0,linetype="dashed",color="pink")+
  theme(legend.position="none")+
  labs(x="P value threshold for PRS",y="Log (OR)\n(per-SD effect of PRS on MS risk)")+
	scale_y_continuous(limits = c(-0.2,1))


png("./outputs/prs_risk_plots_summ.png",res=900,units="in",width=10,height=8)
p
dev.off()


# read in MHC-less scores
# read summary of all prs at r2 of 0.1
all_prs = data.frame()
for(r2 in c(0.001,0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","afr","eur")){
	prs = read_table(paste0("./outputs/prs_risk_",r2,"_",ancestry,".prsice")) %>%
	mutate(anc=ancestry,clump_r2 = r2)
	all_prs <<- bind_rows(all_prs,prs)
}
}

all_prs_nomhc = data.frame()
for(r2 in c(0.001,0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","afr","eur")){
	prs = read_table(paste0("./outputs/prs_risk_nomhc_",r2,"_",ancestry,".prsice")) %>%
	mutate(anc=ancestry,clump_r2 = r2)
	all_prs_nomhc <<- bind_rows(all_prs_nomhc,prs)
}
}


all_prs = all_prs %>%
	mutate(MHC = "MHC") %>%
	bind_rows(all_prs_nomhc %>%
	mutate(MHC = "No MHC"))

# plot r2
p = ggplot(all_prs,
  aes(factor(Threshold),R2,fill=factor(clump_r2)))+
  geom_col(color="black",position=position_dodge())+
  facet_grid(toupper(anc)~MHC)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(x="P value threshold for PRS",
  y=bquote(Nagelkerke~R^2),
  fill=bquote(Clumping~R^2))+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))
png("./outputs/prs_risk_plots_all_clumps_summ_nagelkerke.png",res=900,units="in",width=8,height=6)
p
dev.off()

write_csv(all_prs,"./outputs/prs_risk_all_summ.csv")



# read in best scores
best_prs = data.frame()
for(r2 in c(0.001,0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","afr","eur")){
	prs = read_table(paste0("./outputs/prs_risk_",r2,"_",ancestry,".summary")) %>%
	mutate(anc=ancestry,clump_r2 = r2,MHC = "MHC")
	best_prs <<- bind_rows(best_prs,prs)
}
}

for(r2 in c(0.001,0.01,0.1,0.2,0.4,0.6,0.8)){
for(ancestry in c("sas","afr","eur")){
	prs = read_table(paste0("./outputs/prs_risk_nomhc_",r2,"_",ancestry,".summary")) %>%
	mutate(anc=ancestry,clump_r2 = r2,MHC = "No MHC")
	best_prs <<- bind_rows(best_prs,prs)
}
}

best_prs = best_prs %>%
group_by(anc) %>%
slice_min(`Empirical-P`) %>%
slice_max(PRS.R2)


# get best PRS
write_csv(best_prs,"./outputs/best_risk_prs.csv")


# read in best PRS individual scores
## manual read in
best_prs_individual_level = bind_rows(
	read_table(paste0("./outputs/prs_risk_nomhc_0.8_afr.best")) %>% mutate(anc="AFR"),
	read_table(paste0("./outputs/prs_risk_nomhc_0.8_sas.best")) %>% mutate(anc="SAS"),
	read_table(paste0("./outputs/prs_risk_0.2_eur.best")) %>% mutate(anc="EUR")
)



# read in pheno
phenos = data.frame()
for(ancestry in c("sas","afr","eur")){
	pheno = read_table(paste0("./pheno/reimputed_",ancestry,"_pheno.tsv"))
	phenos <<- bind_rows(phenos,pheno)
}

# combine
best_prs_individual_level = best_prs_individual_level %>%
	left_join(phenos %>% distinct(IID,.keep_all=T),
	by="IID")

# read in covars
covars = data.frame()
for(ancestry in c("sas","afr","eur")){
	covar =  read_table(paste0("./pheno/reimputed_",ancestry,"_covars_with_pcs.tsv"))
	covars <<- bind_rows(covars,covar)
}

# combine
best_prs = best_prs_individual_level %>%
	left_join(covars %>% distinct(IID,.keep_all=T),
	by="IID")


prevalences = 	best_prs %>%
		group_by(anc) %>%
		mutate(prs_quart = Hmisc::cut2(PRS,g=4)) %>%
		 group_by(anc,prs_quart) %>%
		 dplyr::count(MS_status) %>%
		 mutate(prop_cases = n /sum(n)) %>%
		 filter(MS_status==2) %>%
		 group_by(anc) %>%
		 mutate(prs_quart_num = row_number())

p=ggplot(prevalences,aes(prs_quart_num,prop_cases*100,fill=toupper(anc),col=toupper(anc)))+
	geom_line(show.legend=F)+
	geom_point(size=3,shape=21,color="black")+
	theme_bw()+
	scale_fill_brewer(palette="Set1")+
	scale_color_brewer(palette="Set1")+
	labs(y="% of cases in each PRS quartile",fill="Ancestry",x="PRS quartile (best PRS per-ancestry)")

png("./outputs/calibration_quartile_plot.png",res=900,units="in",width=6,height=4)
p
dev.off()

# decile plot
prevalences = 	best_prs %>%
		group_by(anc) %>%
		mutate(prs_dec = Hmisc::cut2(PRS,g=10)) %>%
		 group_by(anc,prs_dec) %>%
		 dplyr::count(MS_status) %>%
		 mutate(prop_cases = n /sum(n)) %>%
		 filter(MS_status==2) %>%
		 group_by(anc) %>%
		 mutate(prs_dec_num = row_number())

p=ggplot(prevalences,aes(prs_dec_num,prop_cases*100,fill=toupper(anc),col=toupper(anc)))+
	geom_line(show.legend=F)+
	geom_point(size=3,shape=21,color="black")+
	theme_bw()+
	scale_fill_brewer(palette="Set1")+
	scale_color_brewer(palette="Set1")+
	labs(y="% of cases in each PRS decile",fill="Ancestry",x="PRS quartile (best PRS per-ancestry)")

png("./outputs/calibration_decile_plot.png",res=900,units="in",width=6,height=4)
p
dev.off()




grouped_prs = 	best_prs %>%
		group_by(anc) %>%
		mutate(prs_quart = Hmisc::cut2(PRS,g=4))


get_model_coefs = function(this_anc){
	model = glm(data = grouped_prs %>%
		filter(anc==this_anc),
		 factor(MS_status) ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + factor(prs_quart),
			family = binomial(link="logit"))
	dat = data.frame(summary(model)$coefficients[c(13:15),]) %>%
		dplyr::select(1,2,4) %>%
		mutate(prs_quartile = c(2:4))
	colnames(dat)[c(1:2)] = c("beta","se")
	dat$or = exp(dat$beta)
	dat$lower_ci = exp(dat$beta - dat$se*1.96)
	dat$upper_ci = exp(dat$beta + dat$se*1.96)
	dat$pval = dat$`Pr...z..`
	rownames(dat) = NULL
	dat %>%
			dplyr::select(or,lower_ci,upper_ci,prs_quartile,pval) %>%
			bind_rows(
				data.frame(or=1,lower_ci = 1,upper_ci=1,prs_quartile=1)
				) %>%
				mutate(anc = this_anc)

}

all_dat = bind_rows(
	get_model_coefs("SAS"),
	get_model_coefs("AFR"),
		get_model_coefs("EUR")
	)

	p=ggplot(all_dat,aes(prs_quartile,or,fill=toupper(anc),col=toupper(anc)))+
		geom_errorbar(mapping = aes(x = prs_quartile,ymin = lower_ci,ymax=upper_ci),width=0.1,color="black")+
		geom_point(size=3,shape=21,color="black")+
		theme_bw()+
		scale_fill_brewer(palette="Set1")+
		scale_color_brewer(palette="Set1")+
		labs(y="OR for MS in each PRS quartile\n(vs lowest quartile)",fill="Ancestry",x="PRS quartile (best PRS per-ancestry)")+
		facet_wrap(~toupper(anc))

	png("./outputs/or_plot.png",res=900,units="in",width=8,height=4)
	p
	dev.off()

p = ggplot(best_prs %>% mutate(MS_status = ifelse(MS_status==2,"MS","Control")),aes(PRS,fill=factor(MS_status)))+
geom_density(alpha=0.8,color="black")+
theme_bw()+
scale_fill_brewer(palette="Set1")+
labs(x="PRS (best PRS for each ancestry)",y="Density",fill="Case/Control status")+
facet_wrap(~toupper(anc))

png("./outputs/histograms_prs.png",res=900,units="in",width=8,height=4)
p
dev.off()
