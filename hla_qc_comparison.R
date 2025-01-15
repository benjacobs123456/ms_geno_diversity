library(tidyverse)

# make gene list
hla_genes = c("A","B","C","DRB1","DQB1","DPB1")

# do QC comparison for each gene
for(i in c(1:6)){
    for(anc in c("sas","afr")){
        gene = hla_genes[i]
        nantes = paste0("/data/scratch/hmy117/hla_allele_nantes_",gene,"_",anc,".tsv")
        kg = paste0("/data/scratch/hmy117/hla_allele_kg_",gene,"_",anc,".tsv")
        snp2hla = paste0("/data/scratch/hmy117/filtered_genotypes_for_susceptibility_gwas_merged_all_chroms_for_gwas_",anc,"_for_hla_imp_omnibus_just_hla_genotypes.raw")

        # read in
        nantes_calls = read_tsv(nantes,col_types = "cccdd")
        kg_calls = read_tsv(kg,col_types = "cccdd")
        snp2hla_calls = read_table(snp2hla) %>% 
            dplyr::select(IID,contains("HLA")) %>% 
            pivot_longer(cols = -IID) %>% 
            separate(name,sep = "_", into = c("part1","part2","part3")) %>% 
            separate(part2,sep = "\\*", into = c("gene","allele")) %>% 
            separate(allele,sep = ":", into = c("field1","field2","field3")) %>% 
            dplyr::rename("sample.id"=IID, "n" = value) %>% 
            mutate(value = paste0(gene,"*",field1,":",field2)) %>% filter(!is.na(field2)) %>% 
            filter(gene == hla_genes[i]) %>%
            dplyr::select(sample.id,value,n)




        # spread long
        nantes_long =  nantes_calls %>%
            pivot_longer(c(allele1,allele2)) %>%
            group_by(sample.id) %>%
            mutate(value = paste0(gene,"*",as.character(value))) %>%
            dplyr::count(value) %>% ungroup()
        kg_long = kg_calls %>%
            pivot_longer(c(allele1,allele2)) %>%
            group_by(sample.id) %>%
            mutate(value = paste0(gene,"*",as.character(value))) %>%
            dplyr::count(value)%>% ungroup()

        # compare frequencies
        nantes_freq = nantes_long %>% group_by(value) %>% summarise(n = sum(n)) %>% mutate(af = n / sum(n))
        kg_freq = kg_long %>% group_by(value) %>% summarise(n = sum(n)) %>% mutate(af = n / sum(n))
        all_freqs = nantes_freq %>% inner_join(kg_freq,by="value")

        p1 = ggplot(all_freqs,aes(af.x,af.y,label = value))+
            geom_point()+
            geom_abline(intercept=0,slope=1,lwd=0.5,linetype="dashed",alpha=0.2)+
            ggrepel::geom_text_repel()+
            theme_bw()+
            labs(x="AF using SHLARC panel",y="AF using 1kg multi-ethnic panel")

        # repeat with snp2hla 
        snp2hla_freqs = snp2hla_calls %>% 
            group_by(value) %>% 
            summarise(ac = sum(n)) %>% 
            ungroup() %>% 
            mutate(total = sum(ac)) %>% 
            mutate(af = ac / total)

        all_freqs = nantes_freq %>% inner_join(snp2hla_freqs,by="value")

        p4 = ggplot(all_freqs,aes(af.x,af.y,label = value))+
            geom_point()+
            geom_abline(intercept=0,slope=1,lwd=0.5,linetype="dashed",alpha=0.2)+
            ggrepel::geom_text_repel()+
            theme_bw()+
            labs(x="AF using SHLARC panel",y="AF using SNP2HLA (with 1kg multi-ethnic panel)")


        # spread long
        combo = nantes_calls %>%
            left_join(kg_calls,by=c("sample.id")) %>%
            mutate(concordant = ifelse(allele1.x == allele1.y & allele2.x == allele2.y,"concordant","discordant"))
        pct = combo %>%
            dplyr::count(concordant) %>%
            mutate(pct = round(n/sum(n)*100,1))
        lab = paste0(pct[pct$concordant=="concordant",]$pct[1],"% concordance")

        p2 = ggplot(combo,aes(prob.x,prob.y,col=concordant))+
            geom_point()+
            geom_abline(intercept=0,slope=1,lwd=0.5,linetype="dashed",alpha=0.2)+
            theme_bw()+
            labs(x="Probability using SHLARC panel",y="Probability using 1kg multi-ethnic panel",col="Concordance")+
            scale_color_brewer(palette="Set1")+
            ggtitle(lab)

        # concordance vs SNP2HLA 
        nantes_snp2hla_concordance = nantes_long %>% 
            mutate(sample.id = as.numeric(sample.id)) %>% 
            inner_join(snp2hla_calls,by=c("sample.id","value")) %>% mutate(concordant = ifelse(n.x==n.y,"yes","no")) %>% 
            group_by(value) %>% dplyr::count(concordant)  %>% mutate(prop = n/sum(n)) %>% filter(concordant == "yes")
        p3 = ggplot(nantes_snp2hla_concordance,aes(value,prop))+geom_col(color="black")+theme_bw()+labs(x="Allele",y="Proportion concordant (HIBAG vs SNP2HLA)")+coord_flip()


        outfile = paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/hla_calls_qc_",gene,"_",anc,"plots.png")
        png(outfile,res=900,units="in",width=10,height=10)
        print(cowplot::plot_grid(p1,p2,p3,p4,ncol=2))
        dev.off()

    }
}
