# Preprocessing GTEx breast tissue including cell type estimation and evaluation ---------------------------
rm(list=ls())
library(data.table)
library(tidyverse)
library(DBI)
library(RSQLite)
library(janitor)
library(ctOWAS)

#source('Github/R/ctOWAS_train_k.R')
setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')

paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
gtex_dir <- file.path(data_dir, "GTEx_Breast")
results_dir <- file.path(paper_dir, "Results")

load(file.path(results_dir, "markers.RData"))
load(file.path(results_dir, "markers_ensembl.RData"))
load(file.path(results_dir, "xCellScore.RData"))

# 1. Sample information - female white participants only

gtex_white <- read_csv(file.path(gtex_dir, "gtex_v8_race.csv")) %>% filter(RACE == "White") %>% pull(SUBJID)
subject_pheno <- fread(file.path(gtex_dir, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
female_white <- intersect(subject_pheno[SEX == 2, SUBJID],gtex_white)

# Breast expression
ensembl38 <- read_csv(file.path(data_dir, "ensembl38.txt"),show_col_types = FALSE) %>%clean_names() %>%distinct(gene_stable_id, gene_name)

breast <- fread(file.path(gtex_dir, "Breast_Mammary_Tissue.v8.normalized_expression.bed")) %>%
  mutate(gene_stable_id = sub("\\..*$", "", gene_id)) %>% inner_join(ensembl38, by = "gene_stable_id") %>%
  group_by(gene_name) %>% filter(n() == 1L) %>% ungroup()

expr <- breast %>% dplyr:: select(all_of(intersect(female_white, names(breast)))) %>% as.data.frame()
rownames(expr) <- breast$gene_id

dim(expr)# check the number
# [1] 24919   125

expr2 <- expr %>%
  tibble::rownames_to_column("gene_stable_id") %>%
  dplyr::mutate(gene_stable_id = sub("\\..*$", "", gene_stable_id)) %>%
  dplyr::inner_join(ensembl38, by = "gene_stable_id") %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(dplyr::n() == 1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(-gene_stable_id) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_name")

dim(expr2)
#[1] 24919   125

# Covariates ---------------------------------------------------------------

cov1=data.frame(fread(file.path(gtex_dir, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")))
cov0=cov1[,c("SUBJID", "AGE")]
cov2=fread(file.path(gtex_dir, "Breast_Mammary_Tissue.v8.covariates.txt"))
cov3=t(cov2[,-1])
colnames(cov3)=data.frame(cov2)[,1]
cov3=data.frame(colnames(cov2)[-1], cov3);colnames(cov3)[1]="SUBJID"
cov3=cov3[,c(1:21,67:69)] # top 15 PEER factors
cov4=data.frame(merge(cov0, cov3, by="SUBJID"))
cov <- cov4%>%filter(SUBJID %in% gtex_white) %>%filter(sex==2) %>%dplyr::select(!sex)
dim(cov)
# [1] 125  24

# Genotype ---------------------------------------------------------------
geno1=fread(file.path(gtex_dir, "shapeit_data_for_predictdb_variants-r2")) # 178698    847
geno=geno1[,c(1:9,match(cov[,1], colnames(geno1))), with=F]
dim(geno)
# [1] 178698    134

# PredictDB elastic-net reference
con <- dbConnect(SQLite(),file.path(gtex_dir, "en_Breast_Mammary_Tissue.db"))
ENextra <- dbReadTable(con, "extra") %>%mutate(gene_id = sub("\\..*$", "", gene))
ENweights <- dbReadTable(con, "weights") %>%mutate(gene_id = sub("\\..*$", "", gene))
dbDisconnect(con)

genID <- intersect(ENextra$gene, breast$gene_id)
G <- length(genID)
G
#G  6443 correct

# Save all data

save(female_white, geno, expr, expr2, genID, cov, ENweights, ensembl38, file = file.path(results_dir, "GTEx_clean.RData"))
