# Train the 3 vs 2 cell-type ctOWAS model on GTEx breast tissue for comparison ---------------------------

# Package dependencies
library(data.table)
library(tidyverse)
library(janitor)
library(ctOWAS)
library(readr)
library(dplyr)
library(glmnet)
library(janitor)
library(tibble)
library(doParallel)
library(dplyr)
setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results")

# load  data
load(file.path(results_dir, "GTEx_clean_pi.RData"))
G=length(genID)
# ------------------------------------------------------------------------------
# 3. Gene-by-gene model training
# ------------------------------------------------------------------------------

# Run the 3-cell-type ctOWAS training loop with a fixed CV seed per gene.

# parallel version
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# G=5

res <- foreach(j = seq_len(G), .packages = c("glmnet", "data.table"),
               .export = c("ctOWAS_train_K"),.errorhandling = "pass") %dopar% {

                 # Current gene and its SNP annotation in the elastic-net reference model.
                 yName=genID[j]
                 xName=ENweights[which(ENweights$gene==yName), "varID"]
                 xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
                 nName=cov$"SUBJID" # women
                 n=length(nName)

                 # Align expression, genotype, covariates, and cell-type proportions by sample ID.
                 yData=t(expr[which(rownames(expr)==yName), match(nName, colnames(expr))])
                 xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with = FALSE]);
                 storage.mode(xData) <- "double"
                 zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
                 # pi2Data=cbind(pi2[match(nName, pi2$SampleID),2], 1-pi2[match(nName, pi2$SampleID),2])
                 pi3Data=pi3_new[match(nName, rownames(pi3_new)),]
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
                 # pi2.complete=pi2Data[cp.idx, ]
                 pi3.complete=pi3Data[cp.idx, ]

                 # Create a reproducible 10-fold CV split for this gene.
                 set.seed(1334 + j*14905)
                 foldid= sample(rep(1:10, length.out =  length(y.complete)))

                 # #  y = y.complete;x = x.complete;pi_k = pi3.complete;cov = z.complete;xNameMatrix = xName.all[xvar0,]; fdrcut=0.25;alpha=0.5; seed=NULL
                 # # ctOWAS K=2
                 # ft_s_2ct <- tryCatch({
                 #   ft_s_2ct  <- ctOWAS_train_K(y = y.complete,x = x.complete,pi_k = pi2.complete,
                 #                            cov = z.complete, xNameMatrix = xName.all[xvar0,], yName=yName,
                 #                            L=10, pip_cutoff=0.1,shared_prior_weight=0.2,
                 #                            cov_prior_weight=1, hard_threshold=T, coverage=0.99)
                 #   w <- cbind(ft_s_2ct$xNameMatrix, ft_s_2ct$W)
                 #   if (!is.null(w) && nrow(w)>0) {weights_s_2ct <- w; intercept_s_2ct = ft_s_2ct$beta.all.models[1,]}
                 # }, error = function(e) {
                 #   cat("ctOWAS_train_K has no predictor.")})

                 # ctOWAS K=3
                 ft_s_3ct1 <- tryCatch({
                   ft_s_3ct  <- ctOWAS_train_K(y = y.complete,x = x.complete,pi_k = pi3.complete,
                                            cov = z.complete, xNameMatrix = xName.all[xvar0,], yName=yName, foldid = foldid)
                   w <- cbind(ft_s_3ct$xNameMatrix, ft_s_3ct$W)
                   if (!is.null(w) && nrow(w)>0) {weights_s_3ct <- w; intercept_s_3ct = ft_s_3ct$beta.all.models[1,]}
                 }, error = function(e) {
                   cat("ctOWAS_train_K has no predictor.")})

                 list(gene = yName, weights_s_3ct = weights_s_3ct, intercept_s_3ct = intercept_s_3ct)
               }


parallel::stopCluster(cl)
foreach::registerDoSEQ()
# ------------------------------------------------------------------------------
# 4. Combine and save gene-level weights
# ------------------------------------------------------------------------------

# Convert any foreach-level errors into standard result records
res <- Map(function(x, gene) if (inherits(x, "error")) list(gene = gene, error_parallel = conditionMessage(x)) else x, res, genID[seq_len(G)])

res_weights_s_3ct <- lapply(res, `[[`, "weights_s_3ct")
res_intercept_s_3ct <- lapply(res, `[[`, "intercept_s_3ct")

gene_names <- vapply(res, `[[`, character(1), "gene")
names(res_weights_s_3ct) <- names(res_intercept_s_3ct)<- gene_names

keep_3ct <- \(x) {x |>
    dplyr::filter(
      dplyr::coalesce(Cell1, 0) != 0 |
        dplyr::coalesce(Cell2, 0) != 0 |
        dplyr::coalesce(Cell3, 0) != 0
    )}


filtered_weights_s_3ct <- res_weights_s_3ct |> dplyr::bind_rows(.id = "gene") |> keep_3ct()
filtered_intercept_s_3ct <- res_intercept_s_3ct |> dplyr::bind_rows(.id = "gene") |> keep_3ct()



save(filtered_weights_s_3ct, filtered_intercept_s_3ct,
     file = file.path(results_dir, "GTEx_Breast_pi3_Comparison.RData"))

