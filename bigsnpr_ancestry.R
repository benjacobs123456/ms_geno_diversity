library(bigsnpr)
library(tidyverse)

# download reference files
DIR = "/data/scratch/hmy117"
all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620968",
                         dir = DIR, fname = "ref_freqs.csv.gz"))
projection <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620953",
                         dir = DIR, fname = "projection.csv.gz"))

# match adams to reference

path_to_geno = "/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/combined_adams_imputed_qc_hg19_chrpos_nodups"
adams_snps = bigreadr::fread2(paste0(path_to_geno,".bim"),
  select = c(1, 4:6),
                   col.names = c("chr", "pos", "a1", "a0")) %>%
  mutate(beta = 1) %>%
  snp_match(all_freq[1:5]) %>%
  print()

# read matched SNPs
rds <- snp_readBed2(paste0(path_to_geno,".bed"),
 ind.col = adams_snps$`_NUM_ID_.ss`)
obj.bigsnp <- snp_attach(rds)
G <- obj.bigsnp$genotypes

# PCA projection
# project individuals (divided by 2) onto the PC space
PROJ <- as.matrix(projection[adams_snps$`_NUM_ID_`, -(1:5)])

correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
all_proj <- big_prodMat(G, sweep(PROJ, 2, correction / 2, '*'),
                        # scaling to get G if beta = 1 and (2 - G) if beta = -1
                        center = 1 - adams_snps$beta,
                        scale = adams_snps$beta)


X <- crossprod(PROJ,
               as.matrix(all_freq[adams_snps$`_NUM_ID_`, -(1:5)]))
               cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
Amat <- cbind(1, diag(ncol(X)))
bvec <- c(1, rep(0, ncol(X)))

# define groups
group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"
grp_fct <- factor(group, unique(group))


# assign to one group
all_centers <- t(X)
all_sq_dist <- apply(all_centers, 1, function(one_center) {
  rowSums(sweep(all_proj, 2, one_center, '-')^2)
})

THR <- 0.005  # you can adjust this threshold
thr_sq_dist <- max(dist(all_centers)^2) * THR / 0.16

cluster <- group[
  apply(all_sq_dist, 1, function(x) {
    ind <- which.min(x)
    if (isTRUE(x[ind] < thr_sq_dist)) ind else NA
  })
]

table(cluster, exclude = NULL)

# ancestry proportions
cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
Amat <- cbind(1, diag(ncol(X)))
bvec <- c(1, rep(0, ncol(X)))

# solve a QP for each projected individual
all_res <- apply(all_proj, 1, function(y) {
  quadprog::solve.QP(
    Dmat = cp_X_pd$mat,
    dvec = crossprod(y, X),
    Amat = Amat,
    bvec = bvec,
    meq  = 1
  )$sol %>%
    tapply(grp_fct, sum) %>%
    round(7)
})

# plot
dat = all_res %>% t() %>%
  data.frame() %>%
  mutate(iid = row_number()) %>%
  pivot_longer(cols = c(1:18))
ggplot(dat,aes(iid,value,fill=name))+
geom_col(position="fill")
