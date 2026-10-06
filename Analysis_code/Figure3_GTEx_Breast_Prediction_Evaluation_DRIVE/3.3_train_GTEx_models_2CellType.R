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

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results")

source('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/R/ctOWAS_train_k.R')

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


res <- foreach(j = seq_len(G), .packages = c("glmnet", "data.table"),
  .export = c("ctOWAS_train_K", "MiXcan"),.errorhandling = "pass") %dopar% {


  # Current gene and its SNP annotation in the elastic-net reference model.
  yName=genID[j]
  xName=ENweights[which(ENweights$gene==yName), "varID"]
  xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
  nName=cov$"SUBJID" # women
  n=length(nName)

  # Align expression, genotype, covariates, and cell-type proportions by sample ID.
  yData=t(expr2[which(rownames(expr)==yName), match(nName, colnames(expr))])
  xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with = FALSE]);
  storage.mode(xData) <- "double"
  zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
  pi1 <- as.numeric(pi2_new[,1][match(nName, rownames(pi2_new))])
  pi2 <- as.numeric(pi2_old[["mean_trim_0.05"]][match(nName, pi2_old$sample)])
  pi1Data <- cbind(Cell1 = pi1,Cell2 = 1 - pi1)
  pi2Data <- cbind(Cell1 = pi2,Cell2 = 1 - pi2)

  # Keep samples with complete genotype and expression data.
  cp.idx=complete.cases(xData) & complete.cases(yData)

  # Keep SNPs with mean genotype > 0.05 among complete samples.
  # Handle the single-SNP case separately so matrix dimensions stay valid.
  px=ncol(xData)
  if (px>1) {
    xvar0=which(apply(xData[cp.idx,], 2, function(f) mean(f)>0.05))
    x.complete=xData[cp.idx,xvar0]
  }
  if (px==1) {
    if (mean(xData[cp.idx, 1]) > 0.05) {
      xvar0 <- 1L
      x.complete <- matrix(xData[cp.idx, 1], ncol = 1)
    } else {
      xvar0 <- integer(0)
      x.complete <- matrix(numeric(0), nrow = sum(cp.idx), ncol = 0)
    }
  }
  if (ncol(x.complete) == 0 ||is.null(nrow(x.complete))) {next}

  z.complete=zData[cp.idx,]
  xz.complete=as.matrix(cbind(x.complete, z.complete))
  y.complete=yData[cp.idx]
  pi1.complete=pi1Data[cp.idx, ]
  pi2.complete=pi2Data[cp.idx, ]
  # Create a reproducible 10-fold CV split for this gene.
  set.seed(1334 + j*14905)
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
res <- Map(function(x, gene) if (inherits(x, "error")) list(gene = gene, error_parallel = conditionMessage(x)) else x, res, genID[seq_len(G)])

res_weights_s <- lapply(res, `[[`, "weights_s")
res_intercept_s <- lapply(res, `[[`, "intercept_s")
res_weights_m_new <- lapply(res, `[[`, "weights_m_new")
res_intercept_m_new <- lapply(res, `[[`, "intercept_m_new")
res_weights_m_old <- lapply(res, `[[`, "weights_m_old")
res_intercept_m_old <- lapply(res, `[[`, "intercept_m_old")
res_weights_p <- lapply(res, `[[`, "weights_p")

gene_names <- vapply(res, `[[`, character(1), "gene")
names(res_weights_s) <- names(res_weights_m_new) <-  names(res_weights_m_new) <-names(res_weights_p) <- names(res_intercept_s) <-
  names(res_intercept_m_new) <-   names(res_intercept_m_new) <- gene_names

keep_2ct <- \(x) {
  if (is.null(x)) return(NULL)

  keep <- dplyr::coalesce(x$Cell1, 0) != 0 |
    dplyr::coalesce(x$Cell2, 0) != 0

  x[keep, , drop = FALSE]
}
filtered_weights_s <- res_weights_s |>
  lapply(keep_2ct) |>
  dplyr::bind_rows()

filtered_intercept_s <- res_intercept_s |>
  dplyr::bind_rows(.id = "gene")

filtered_weights_m_new <- res_weights_m_new |>
  lapply(keep_2ct) |>
  dplyr::bind_rows()

filtered_intercept_m_new <- res_intercept_m_new |>
  dplyr::bind_rows(.id = "gene")

filtered_weights_m_old <- res_weights_m_old |>
  lapply(keep_2ct) |>
  dplyr::bind_rows()

filtered_intercept_m_old <- res_intercept_m_old |>
  dplyr::bind_rows(.id = "gene")


filtered_weights_p <- res_weights_p |>
  dplyr::bind_rows() |>
  dplyr::filter(dplyr::coalesce(tissue, 0) != 0)


save(filtered_weights_s, filtered_intercept_s,
     filtered_weights_m_new, filtered_intercept_m_new,
     filtered_weights_m_old, filtered_intercept_m_old,
     filtered_weights_p, file = file.path(results_dir, "GTEx_Breast_pi2_ThreeModel_Summary.RData"))

length(unique(filtered_weights_s$gene))
length(unique(filtered_weights_m_new$gene))
length(unique(filtered_weights_m_old$gene))
length(unique(filtered_weights_p$gene))




