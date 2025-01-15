library(tidyverse)
setwd("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/")


# function for power
find_power = function(
    disease_af = 0.3,
    control_af = 0.15,
    n_case = 110,
    n_control = 4994,
    n_iter = 10,
    p_thresh = 5e-8){


    pvals = list()  
    betas = list()
    for(i in c(1:n_iter)){

        # simulate genotypes
        case_genotypes = data.frame(geno = rbinom(prob = disease_af, size = 2, n= n_case)) %>% mutate(dis_status = 1)
        ctrl_genotypes = data.frame(geno = rbinom(prob = control_af, size = 2, n= n_control)) %>% mutate(dis_status = 0)

        # combine
        dat = bind_rows(case_genotypes,ctrl_genotypes)

        # get exact rr 

        # model
        model_summ = broom::tidy(glm(data = dat, dis_status ~ geno,family = binomial(link="logit"))) %>% 
            filter(term == "geno")
        pvals[[i]] = model_summ$p.value
        betas[[i]] = model_summ$estimate
    }

    
    pvals = unlist(pvals)
    power = sum(pvals<5e-8)/n_iter
    or = exp(median(unlist(betas),na.rm=T))

    c(or,power)
  }

# loop through parameters 
disease_af = seq(0.05,0.5,by=0.05)
control_af = seq(0.05,0.5,by=0.05)
n_case = c(110,177)
n_control = c(4994,6785)



params = expand.grid(disease_af,control_af,n_case,n_control) %>% 
filter(
    (Var3 == 110 & Var4 == 4994 ) |
    (Var3 == 177 & Var4 == 6785 ) 
) 
powers = list()
ors = list()
for(i in c(1:nrow(params))){
  message(i," of ",nrow(params))
  res = find_power(params$Var1[i], params$Var2[i], params$Var3[i],params$Var4[i],n_iter = 10, p_thresh = 5e-8)
  powers[[i]] = res[2]
  ors[[i]] = res[1]
}
params$power = unlist(powers)
params$ancestry = ifelse(params$Var3==110,"AFR","SAS")
params$OR = unlist(ors)
params$lab = paste0("Power=",round(params$power,2)*100,"%\nOR=",round(params$OR,2))
p = ggplot(params,aes(Var1,Var2,fill=power,label=lab))+
    geom_tile(color="black")+
    facet_wrap(~ancestry)+
    geom_text(size=2)+
    theme_bw()+
    labs(x="Disease AF",y="Control AF")+
    scale_fill_gradient(low = "blue",high="orange")+
    theme(legend.position="none")

png("./outputs/power_curves_risk.png",res=900,units="in",width=16,height=8)
print(p)
dev.off()



