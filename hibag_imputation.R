library(HIBAG)
library(tidyverse)

# make gene list 
hla_genes = c("A","B","C","DRB1","DQB1","DPB1")

# define this gene 
args = commandArgs(trailingOnly=T)
i = as.numeric(args[1])
gene = hla_genes[i]
ancestry = args[2]

# Load the model
kg_model.list = get(load("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/multiethnic_IKMB_1KG.RData"))
kg_model = kg_model.list[[gene]]

# load nantes shlarc model 
nantes_model = readRDS(paste0("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/hibag_models_nantes/model_HLA_",gene,"_bagging_100_SHLARCdata_allSNPs_anonymized.RData"))

# impute hla alleles

# impute with 1kg (hg19)
# Import plink file - hg19
plink_file = paste0("/data/scratch/hmy117/risk_gwas_genotypes_hg19_hibag",ancestry)

genotypes <- hlaBED2Geno(bed.fn=paste0(plink_file,".bed"),
  fam.fn=paste0(plink_file,".fam"),
  bim.fn=paste0(plink_file,".bim"))
summary(genotypes)

# HLA imputation at HLA-x
model_kg = hlaModelFromObj(kg_model)


# best-guess genotypes and all posterior probabilities
predicted_genos = hlaPredict(model_kg, genotypes, type="response+prob", verbose=T)

allele_calls = data.frame(predicted_genos$value)
out_file = paste0("/data/scratch/hmy117/hla_allele_kg_",gene,"_",ancestry,".tsv")
write_tsv(allele_calls,out_file)
out_file = paste0("/data/scratch/hmy117/hla_allele_kg_",gene,"_",ancestry,".rds")
saveRDS(predicted_genos,out_file)


# repeat with nantes model (hg38)
plink_file_hg38 = paste0("/data/scratch/hmy117/risk_gwas_genotypes_hg38_hibag",ancestry)

genotypes_hg38 <- hlaBED2Geno(bed.fn=paste0(plink_file_hg38,".bed"),
  fam.fn=paste0(plink_file_hg38,".fam"),
  bim.fn=paste0(plink_file_hg38,".bim"),assembly="hg38")
summary(genotypes_hg38)

# HLA imputation at HLA-x
model_nantes = hlaModelFromObj(nantes_model)

# best-guess genotypes and all posterior probabilities
predicted_genos = hlaPredict(model_nantes, genotypes_hg38, type="response+prob", verbose=T)

allele_calls = data.frame(predicted_genos$value)
out_file = paste0("/data/scratch/hmy117/hla_allele_nantes_",gene,"_",ancestry,".tsv")
write_tsv(allele_calls,out_file)
out_file = paste0("/data/scratch/hmy117/hla_allele_nantes_",gene,"_",ancestry,".rds")
saveRDS(predicted_genos,out_file)


