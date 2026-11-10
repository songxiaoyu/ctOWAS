# Train the 2-cell-type models on GTEx breast tissue ---------------------------
rm(list=ls())
# Package dependencies
library(data.table)
# library(xCell)
library(tidyverse)
library(janitor)
library(MiXcan)
library(readr)
library(dplyr)
library(glmnet)
library(janitor)
library(tibble)
library(foreach)
library(doParallel)
library(dplyr)
library(ctOWAS)


paper_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
results_dir <- file.path(paper_dir, "Results")


# ------------------------------------------------------------------------------
# 2. Gene-by-gene model training

# load  data
load(file.path(results_dir, "GTEx_clean_pi.RData"))
G=length(genID)
#G=2
# parallel version
n_cores <- max(1L, parallel::detectCores() - 2L)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# SLC4A7:   yName="ENSG00000033867.16"; j=4616
# AC244021.1: ENSG00000227082.1; j=3031
# SRGAP2C: ENSG00000171943.11; j=3030
res <- foreach(j = seq_len(G), .packages = c("glmnet", "data.table"),
  .export = c("ctOWAS_train_K", "MiXcan"),.errorhandling = "pass") %dopar% {

  # Current gene and its SNP annotation in the elastic-net reference model.
  yName=genID[j]
  xName=ENweights[which(ENweights$gene==yName), "varID"]
  xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
  nName=cov$"SUBJID" # women
  n=length(nName)

  # Align expression, genotype, covariates, and cell-type proportions by sample ID.
  yData=t(expr[which(rownames(expr)==yName), match(nName, colnames(expr))])
  xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with = FALSE]); storage.mode(xData) <- "double"
  zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
  pi1 <- as.numeric(pi2_new[,1][match(nName, rownames(pi2_new))])
  pi2 <- as.numeric(pi2_old[["mean_trim_0.05"]][match(nName, pi2_old$sample)])
  pi1Data <- cbind(Cell1 = pi1,Cell2 = 1 - pi1)
  pi2Data <- cbind(Cell1 = pi2,Cell2 = 1 - pi2)

  # Keep samples with complete genotype and expression data.
  cp.idx=complete.cases(xData) & complete.cases(yData)

  # clean x.complete: Keep SNPs with MAF genotype > 0.05 among complete samples.
  xvar0 <- which(apply(xData[cp.idx, , drop = FALSE], 2, function(f) {
    af <- mean(f) / 2
    min(af, 1 - af) > 0.05
  }))
  x.complete <- xData[cp.idx, xvar0, drop = FALSE]
  if (ncol(x.complete) == 0) return(NULL)

  z.complete=zData[cp.idx,]
  xz.complete=as.matrix(cbind(x.complete, z.complete))
  y.complete=yData[cp.idx]
  pi1.complete=pi1Data[cp.idx, ]
  pi2.complete=pi2Data[cp.idx, ]

  # Create a reproducible 10-fold CV split for this gene.
  set.seed(1335 + j*123)
  foldid= sample(rep(1:10, length.out =  length(y.complete)))

  #  y = y.complete;x = x.complete;pi_k = pi1.complete;cov = z.complete; xNameMatrix = xName.all[xvar0,]; alpha=0.5;
  # ctOWAS
  ft_s1 <- tryCatch({
    ft_s  <- ctOWAS_train_K(y = y.complete,x = x.complete,
                            pi_k = pi1.complete, cov = z.complete,foldid=foldid,
                            xNameMatrix = xName.all[xvar0,], yName=yName)

    w <- cbind(ft_s$xNameMatrix, ft_s$W)
    if (!is.null(w) && nrow(w)>0) {weights_s <- w; intercept_s = ft_s$intercept}
  }, error = function(e) {warning("ctOWAS_train_K failed: ", conditionMessage(e))
    NULL
  })

  # MiXcan - new pi
  ft_m <- MiXcan(y = y.complete, x = x.complete, pi = pi1.complete[,1], cov = z.complete,
                 xNameMatrix = xName.all[xvar0,], yName=yName, foldid = foldid)
  weights_m_new <- data.frame(ft_m$xNameMatrix, Cell1=ft_m$beta.SNP.cell1$weight, Cell2=ft_m$beta.SNP.cell2$weight)
  intercept_m_new= ft_m$beta.all.models[1,-1]

  # MiXcan - old pi
  ft_m <- MiXcan(y = y.complete, x = x.complete, pi = pi2.complete[,1], cov = z.complete,
                 xNameMatrix = xName.all[xvar0,], yName=yName, foldid = foldid)
  weights_m_old <- data.frame(ft_m$xNameMatrix, Cell1=ft_m$beta.SNP.cell1$weight, Cell2=ft_m$beta.SNP.cell2$weight)
  intercept_m_old= ft_m$beta.all.models[1,-1]
  # PrediXcan
  weights_p = data.frame(ft_m$xNameMatrix, tissue=ft_m$beta.all.models[2:(ncol(x.complete)+1), 1])

  list(gene = yName, weights_s = weights_s, intercept_s = intercept_s,
    weights_m_new = weights_m_new, weights_m_old = weights_m_old,
    intercept_m_new = intercept_m_new, intercept_m_old = intercept_m_old,
    weights_p = weights_p)
}


parallel::stopCluster(cl)
foreach::registerDoSEQ()


# ------------------------------------------------------------------------------
# 4. Combine and save gene-level weights
# ------------------------------------------------------------------------------
# Convert any foreach-level errors into standard result records
# ------------------------------------------------------------------------------
# 4. Combine and save gene-level weights
# ------------------------------------------------------------------------------
stopifnot(length(res) == G, length(genID) >= G)

# Standardize errors/NULL results and restore missing gene IDs
res <- Map(function(x, gene) {
  gene <- as.character(gene)

  if (inherits(x, "error")) {
    x <- list(error_parallel = conditionMessage(x))
  } else if (is.null(x)) {
    x <- list()
  }

  if (length(x$gene) != 1L || is.na(x$gene) || !nzchar(x$gene)) {
    x$gene <- gene
  }

  x
}, res, genID[seq_len(G)])

# Pull each component, tagged by gene
gene_names <- vapply(res, function(x) as.character(x$gene), character(1))
pull <- function(field) setNames(lapply(res, `[[`, field), gene_names)

# Empty templates preserve columns when no weights are available
empty_2ct <- tibble::tibble(
  gene = character(), Cell1 = double(), Cell2 = double()
)
empty_p <- tibble::tibble(
  gene = character(), tissue = double()
)

combine_2ct <- function(field) {
  dplyr::bind_rows(empty_2ct, dplyr::bind_rows(pull(field))) |>
    dplyr::filter(
      dplyr::coalesce(Cell1, 0) != 0 |
        dplyr::coalesce(Cell2, 0) != 0
    )
}

# Weights
filtered_weights_s     <- combine_2ct("weights_s")
filtered_weights_m_new <- combine_2ct("weights_m_new")
filtered_weights_m_old <- combine_2ct("weights_m_old")

# Intercepts
filtered_intercept_s <- dplyr::bind_rows(
  pull("intercept_s"), .id = "gene"
)
filtered_intercept_m_new <- dplyr::bind_rows(
  pull("intercept_m_new"), .id = "gene"
)
filtered_intercept_m_old <- dplyr::bind_rows(
  pull("intercept_m_old"), .id = "gene"
)

# PrediXcan
filtered_weights_p <- dplyr::bind_rows(
  empty_p, dplyr::bind_rows(pull("weights_p"))
) |>
  dplyr::filter(dplyr::coalesce(tissue, 0) != 0)

save(
  filtered_weights_s, filtered_intercept_s,
  filtered_weights_m_new, filtered_intercept_m_new,
  filtered_weights_m_old, filtered_intercept_m_old,
  filtered_weights_p,
  file = file.path(results_dir, "GTEx_Breast_pi2_ThreeModel_Summary.RData")
)

sapply(
  list(
    ctOWAS = filtered_weights_s,
    MiXcan_new = filtered_weights_m_new,
    MiXcan_old = filtered_weights_m_old,
    PrediXcan = filtered_weights_p
  ),
  function(d) length(unique(d$gene))
)
