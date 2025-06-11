library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")
# set maf
maf = 0.17

sample_size = 100

beta = 0.5
n_iter=100
p_thresh = 5e-8

cv = 1

# one-tailed P for positive beta
find_power = function(beta=0.5,maf=0.4,sample_size=100,type="replication"){
beta_outs = list()
pval_outs = list()
for(i in c(1:n_iter)){

  # simulate genotypes
  genotypes = data.frame(geno = rbinom(prob = maf, size = 2, n= sample_size))
	table(genotypes$geno)
  # simulate normally-distributed outcome with mean = 1 and cv 1 in homs
  homs = genotypes %>% filter(geno == 0)
  homs$pheno = rnorm(mean = 1, sd = 1, n = nrow(homs))
  hets = genotypes %>% filter(geno == 1)
  hets$pheno = rnorm(mean = 1 + beta, sd = 1 + beta , n = nrow(hets))
  rare_homs = genotypes %>% filter(geno == 2)
  rare_homs$pheno = rnorm(mean = 1 + 2*beta, sd = 1 + 2*beta, n = nrow(rare_homs))

  # combine
  dat = bind_rows(homs,hets,rare_homs)

  # normalise outcome
  dat$pheno = RNOmni::RankNorm(dat$pheno)

  # model
  model_summ = summary(lm(data = dat, pheno ~ geno))$coefficients

	beta_out = model_summ[2,1]
  beta_outs[[i]] = beta_out

	pval_out = model_summ[2,4]
  pval_outs[[i]] = pval_out

  }

	empirical_betas = unlist(beta_outs)
	empiric_p = 1 - pnorm(mean(empirical_betas) /sd(empirical_betas) )
  crude_one_tailed_power = sum(empirical_betas>0) / n_iter
	power = sum(unlist(pval_outs)<p_thresh) / n_iter

  if(type=="discovery"){
		power
	} else {
		crude_one_tailed_power
	}
}


betas = seq(0.1,0.5,by=0.1)
mafs = seq(0.1,0.5,by=0.1)
sample_size = c(100,200,300)
params = expand.grid(b = betas, maf = mafs, n = sample_size)
n_iter=1000
p_thresh = 5e-8
powers = list()
for(i in c(1:nrow(params))){
  message(i," of ",nrow(params))
  powers[[i]] = find_power(params$b[i], params$maf[i], params$n[i],type="discovery")
}

params$power = unlist(powers)

p = ggplot(params,aes(maf,power,fill=factor(n)))+
geom_line()+
geom_point(color="black",shape=21,size=3)+
facet_wrap(~paste0("Beta=",b),nrow=1)+
theme_bw()+
scale_fill_brewer(palette="Set1")+
labs(x="MAF",y="Power",fill="N")

png("./outputs/power_curves_discovery_severity.png",res=900,units="in",width=10,height=3)
print(p)
dev.off()

# repeat for replication
p_thresh = 0.05
powers = list()
for(i in c(1:nrow(params))){
  message(i," of ",nrow(params))
  powers[[i]] = find_power(params$b[i], params$maf[i], params$n[i])
}

params$power = unlist(powers)

p = ggplot(params,aes(maf,power,fill=factor(n)))+
geom_line()+
geom_point(color="black",shape=21,size=3)+
facet_wrap(~paste0("Beta=",b),nrow=1)+
theme_bw()+
scale_fill_brewer(palette="Set1")+
labs(x="MAF",y="Power",fill="N")

png("./outputs/power_curves_replication_severity.png",res=900,units="in",width=10,height=3)
print(p)
dev.off()


