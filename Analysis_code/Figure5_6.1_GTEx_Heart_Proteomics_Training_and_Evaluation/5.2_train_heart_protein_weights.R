# Train 2-component PWAS MiXcan models on GTEx heart proteomics -----------------
#
# Goal:
#   Build protein prediction weights for a downstream PWAS/ctOWAS analysis.
#   Each protein is modeled separately using cis genotype predictors, covariates,
#   and two cell-fraction components:
#     1. Cardiomyocytes
#     2. Others = Fibroblasts + Endothelial + Smooth muscle/pericyte + Immune
#
# Workflow steps:
#   Step 1. Set paths and load protein, cell-fraction, covariate, and EA inputs.
#   Step 2. Align all training inputs by GTEx donor ID.
#   Step 3. Train one PWAS MiXcan model per protein:
#           - load matching GTEx heart genotype dosage
#           - keep cis SNPs within +/- 1 Mb when using chromosome-wide dosage
#           - remove incomplete samples and low-dosage SNPs
#           - fit MiXcan_train_K with K = 2 cell components
#   Step 4. Save nonzero SNP weights and a skipped-protein log.
#
# Helper functions (sample IDs, protein/cell-fraction/annotation loading,
# genotype dosage parsing, output formatting) are defined elsewhere and must be
# on the search path before this script runs.
#
# This script does not read GWAS summary statistics and does not run liftover.

# ------------------------------------------------------------------------------
# Step 1. Set paths and load inputs
# ------------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ctOWAS)
library(MiXcan)
library(foreach)
library(doParallel)



setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
source('Github/Analysis_code/Figure6_GTEx_Heart_Proteomics_Training_and_Evaluation/helper_function.R')

# set directory

paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
heart_dir <- file.path(data_dir, "GTEx_Heart")
protein_dir <- file.path(heart_dir, "GTEx_Pi_Estimate")
results_dir <- file.path(paper_dir, "Results", "heart_protein_training")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# find data

protein <- load_protein_data(file.path(protein_dir,"Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData"))
protein_sample_cols <- names(protein)[-(1:5)]


cov_data <- read_covariates(file.path(paper_dir, "IntermediateResults", "covariate_EA_with_age.txt")) %>%
  distinct(DonorID, .keep_all = TRUE)
# GWS: Paper_ctOWAS/New generated files/codes/by_chr_nomiss -- unpruned dosage files w/o missingness.
geno_raw_dir <- file.path(paper_dir, "IntermediateResults", "WGS", "WGS_pruned_by_chr_500kb_1_r2_0.8")

sample_map <- data.frame(
  RawSampleID = protein_sample_cols,
  SampleID = normalize_sample_id(protein_sample_cols),
  DonorID = donor_id(protein_sample_cols),
  stringsAsFactors = FALSE
) %>% distinct(DonorID, .keep_all = TRUE)

# estimate pi
library(xCell2)

markers <- list(
  Cardiomyocytes = c("ENSG00000198626", "ENSG00000092054"),
  Others = c("ENSG00000011465", "ENSG00000110799", "ENSG00000261371", "ENSG00000133392", "ENSG00000065534",
              "ENSG00000069431", "ENSG00000170458", "ENSG00000173372", "ENSG00000163751")
)


# Cardiomyocytes
# [1] "RYR2" "MYH7"
#
# Others
# [1] "DCN"    "VWF"    "PECAM1" "MYH11"
# [5] "MYLK"   "ABCC9"  "CD14"   "C1QA"
# [9] "CPA3"

# Prepare numeric matrix: genes × samples
protein2 <- as.data.frame(protein)
stopifnot("gene_id" %in% names(protein2))
gene_ids <- sub("\\..*$", "", trimws(as.character(protein2$gene_id)))
sample_cols <- grep("^GTEX[.-]", names(protein2), value = TRUE)
if (!length(sample_cols)) stop("No GTEX sample columns found.")
# This must line up 1:1 with protein_sample_cols (used above for the
# RawSampleID -> DonorID sample map) or the pi_heart rownames assigned below
# will be silently mismatched against the wrong samples.
stopifnot("sample_cols must match protein_sample_cols for correct sample alignment" =
            identical(sample_cols, protein_sample_cols))

protein_mat <- as.matrix(protein2[, sample_cols, drop = FALSE])
storage.mode(protein_mat) <- "numeric"
keep <- !is.na(gene_ids) & nzchar(gene_ids)
protein_mat <- protein_mat[keep, , drop = FALSE]
rownames(protein_mat) <- gene_ids[keep]

if (!nrow(protein_mat)) stop("No rows with valid gene IDs.")
if (any(!is.finite(protein_mat))) stop("Protein matrix contains NA/NaN/Inf; handle missing values before scoring.")

# Average duplicate gene IDs on the supplied abundance scale
if (anyDuplicated(rownames(protein_mat))) {
  gene_n <- table(rownames(protein_mat))
  protein_mat <- rowsum(protein_mat, group = rownames(protein_mat), reorder = FALSE)
  protein_mat <- sweep(protein_mat, 1, as.numeric(gene_n[rownames(protein_mat)]), "/")
}

# Check marker availability
markers_present <- lapply(markers, intersect, y = rownames(protein_mat))
print(data.frame(Cell_type = names(markers), Requested = lengths(markers), Available = lengths(markers_present)))
missing_markers <- lapply(markers, setdiff, y = rownames(protein_mat))
if (any(lengths(missing_markers) > 0L)) print(missing_markers)

# Generate three signatures per cell type; retain identical signatures when necessary
make_xcell2_signatures <- function(marker_list, available_genes, n_signatures = 3L, fraction = 0.8, min_size = 2L) {
  stopifnot(n_signatures >= 3L, n_signatures == as.integer(n_signatures), fraction > 0, fraction <= 1, min_size >= 1L)
  marker_list <- lapply(marker_list, intersect, y = available_genes)
  bad <- lengths(marker_list) < min_size
  if (any(bad)) stop("Fewer than ", min_size, " markers available for: ", paste(names(marker_list)[bad], collapse = ", "))

  signatures <- list()
  for (cell_type in names(marker_list)) {
    genes <- marker_list[[cell_type]]
    n_genes <- length(genes)
    signature_size <- min(n_genes, max(min_size, floor(fraction * n_genes)))
    candidate <- lapply(seq_len(n_signatures), function(j) {
      start <- floor((j - 1L) * n_genes / n_signatures)
      index <- (start + seq_len(signature_size) - 1L) %% n_genes + 1L
      genes[index]
    })
    names(candidate) <- paste0(cell_type, "#custom", seq_len(n_signatures))
    signatures <- c(signatures, candidate)
  }
  signatures
}

custom_signatures <- make_xcell2_signatures(markers, rownames(protein_mat), n_signatures = 3L, fraction = 0.8, min_size = 2L)
cell_names <- names(markers)
spill_identity <- diag(length(cell_names))
dimnames(spill_identity) <- list(cell_names, cell_names)

# Placeholder parameters are appropriate only for the raw-score call below
custom_params <- data.frame(celltype = cell_names, a = 1, b = 1, m = 1, n = 0, stringsAsFactors = FALSE)
custom_xcell2 <- methods::new("xCell2Object", signatures = custom_signatures, dependencies = list(),
                              params = custom_params, spill_mat = spill_identity, genes_used = rownames(protein_mat))
methods::validObject(custom_xcell2)

# Confirm three usable signatures per cell type
signature_counts <- table(sub("#.*$", "", names(custom_signatures)))
print(signature_counts)
stopifnot(all(signature_counts >= 3L), all(lengths(custom_signatures) >= 2L))

# Calculate scores using the full measured gene background
xCellScore <- xCell2::xCell2Analysis(mix = protein_mat, xcell2object = custom_xcell2, minSharedGenes = 0.8,
                                     rawScores = TRUE, spillover = FALSE, BPPARAM = BiocParallel::SerialParam())
xCellScore <- xCellScore[cell_names, , drop = FALSE]
print(xCellScore[, seq_len(min(6L, ncol(xCellScore))), drop = FALSE])

# Save to your existing results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
save(xCellScore, file = file.path(results_dir, "xCellScore2_GTEx_Heart.RData"))


source('Github/R/deNet_composition_K.R')

cell_type <- list(Cardiomyocytes = "Cardiomyocytes", Others = "Others")
pi_prior <- estimate_prior(xCellScore = xCellScore, cell_types = cell_type, mu = c(0.5, 0.5))
# cardiomyocyte: 50%:
# cite: https://www.nature.com/articles/s41586-020-2797-4?utm_source=chatgpt.com


pi_once=deNet_composition_K_once(expr =t(protein_mat) , composition = pi_prior,
                                 marker_list=markers)

pi_heart <- pi_estimation_K(expr = protein_mat, prior = pi_prior, marker_list=markers,
                           n_iteration = 100, seed = 1, verbose = F)


save(pi_heart, file = file.path(results_dir, "pi_heart_GTEx.RData"))



# ------------------------------------------------------------------------------
# Step 2. Align all training inputs by GTEx donor ID
# ------------------------------------------------------------------------------

pi_heart = pi_heart %>% as.data.frame() %>%rownames_to_column("RawSampleID")
# Final training sample set: donors with protein abundance, cell fractions, and
# covariates available. The covariate file is already EA-filtered; the separate
# EA keep file is applied only when it has enough overlap with these samples.
sample_map <- sample_map %>%
  inner_join(pi_heart, by = "RawSampleID", suffix = c("", "_pi")) %>%
  inner_join(cov_data, by = "DonorID", suffix = c("", "_cov"))

cov_cols <- setdiff(names(sample_map),c("RawSampleID", "SampleID", "DonorID", "SampleID_pi", "SampleID_cov","Cardiomyocytes", "Others"))
cov_cols <- cov_cols[vapply(sample_map[cov_cols], is.numeric, logical(1))]
cell_type_cols <- c("Cardiomyocytes","Others")
protein_ann <- annotation_from_protein(protein)
protein_ann <- fill_missing_annotation(protein_ann, file.path(data_dir, "ensembl38.txt"))
protein_ann <- protein_ann %>% distinct(gene_id, .keep_all = TRUE)

# protein_ann=protein_ann[1:20,]
# ------------------------------------------------------------------------------
# Step 3. Train one PWAS MiXcan model per protein
# ------------------------------------------------------------------------------



process_protein <- function(j) {
  print(j)
  target <- protein_ann[j, ]
  gene_id <- target$gene_id
  errors <- character()

  empty_result <- function(reason) list(gene_id = gene_id, weights_s = NULL, intercept_s = NULL,
                                        weights_m = NULL, intercept_m = NULL, weights_p = NULL, skip = reason)

  # cat("Processing", j, "of", nrow(protein_ann), gene_id, "\n")

  if (is.na(target$chr) || is.na(target$start) || is.na(target$end))
    return(empty_result("missing_annotation"))

  geno_target <- read_genotype_for_target(target_id = gene_id, chr = target$chr,
                                          start = target$start, end = target$end)
  if (is.null(geno_target)) return(empty_result("missing_or_empty_genotype_raw"))

  y_row <- which(protein$gene_id == gene_id)
  if (length(y_row) != 1L) return(empty_result("protein_row_not_unique"))

  donors <- intersect(sample_map$DonorID, geno_target$x$DonorID)
  if (length(donors) < 20L) return(empty_result("too_few_genotyped_samples"))

  aligned <- sample_map[match(donors, sample_map$DonorID), ]
  geno_rows <- geno_target$x[match(donors, geno_target$x$DonorID), ]

  y <- as.numeric(protein[y_row, aligned$RawSampleID, drop = TRUE])
  x <- as.matrix(geno_rows[, geno_target$snp_annot$raw_col, drop = FALSE])
  z <- as.matrix(aligned[, cov_cols, drop = FALSE])
  pi_k <- as.matrix(aligned[, cell_type_cols, drop = FALSE])

  storage.mode(x) <- "numeric"
  rownames(x) <- donors
  colnames(x) <- geno_target$snp_annot$varID

  complete_idx <- complete.cases(y) & complete.cases(x) & complete.cases(z) & complete.cases(pi_k)
  if (sum(complete_idx) < 20L) return(empty_result("too_few_complete_samples"))

  x <- x[complete_idx, , drop = FALSE]
  y <- y[complete_idx]
  z <- z[complete_idx, , drop = FALSE]
  pi_k <- pi_k[complete_idx, , drop = FALSE]

  keep_snps <- colMeans(x, na.rm = TRUE) > 0.05
  if (!any(keep_snps)) return(empty_result("no_snps_after_maf_filter"))

  x <- x[, keep_snps, drop = FALSE]
  snp_annot <- geno_target$snp_annot[match(colnames(x), geno_target$snp_annot$varID), ]

  x_name_matrix <- snp_annot %>%
    transmute(varID, gene_id = target$gene_id, gene_name = target$gene_name, chr, pos,
              ref_allele, eff_allele, dosed_allele)

  set.seed(1334 + j * 149053)
  foldid <- sample(rep(seq_len(10L), length.out = length(y)))


  weights_s <- intercept_s <- weights_m <- intercept_m <- weights_p <- NULL
  # ctOWAS
  fit_s <- tryCatch(
    ctOWAS_train_K(y = y, x = x, pi_k = pi_k, cov = z, xNameMatrix = x_name_matrix,
                    yName = gene_id, foldid = foldid),
    error = function(e) {
      errors <<- c(errors, paste0("ctOWAS_train_K_failed: ", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(fit_s)) {
    weights_s <- combine_weights(fit = fit_s, snp_annot=snp_annot, target=target,
                                 cell_type_cols=cell_type_cols, model = "ctOWAS")
    weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
    weights_s <- weights_s[rowSums(abs(weights_s[, weight_cols, drop = FALSE]) > 2e-16) > 0, ]
    if (nrow(weights_s)==0) {weights_s <- NULL} else {
      intercept_s <- extract_intercepts(fit = fit_s, target, cell_type_cols, model = "ctOWAS")}
  }

  # MiXcan
  fit_m <- tryCatch(
    MiXcan::MiXcan(y = y, x = x, pi = pi_k[, 1], cov = z, xNameMatrix = x_name_matrix,
           yName = gene_id, foldid = foldid),
    error = function(e) {
      errors <<- c(errors, paste0("MiXcan_train_failed: ", conditionMessage(e)))
      NULL
    }
  )


  if (!is.null(fit_m) && !identical(fit_m$type, "NoPredictor")) {
    weights_m <- combine_weights(fit = fit_m, snp_annot=snp_annot, target=target,
                                 cell_type_cols=cell_type_cols,model = "MiXcan")
    weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
    weights_m <- weights_m[rowSums(abs(weights_m[, weight_cols, drop = FALSE]) > 2e-16) > 0, ]

    if (!nrow(weights_m)) {weights_m <- NULL} else {
      intercept_m <- extract_intercepts(fit = fit_m, target, cell_type_cols, model = "MiXcan")
    }
  }

  # PrediXcan weights from the MiXcan fit
  if (!is.null(fit_m)) {
    weights_p <- combine_weights(fit = fit_m, snp_annot, target, cell_type_cols, model = "PrediXcan")
    weights_p <- weights_p[abs(weights_p$tissue) > 2e-16, ]
    if (!nrow(weights_p)) weights_p <- NULL
  }

  list(gene_id = gene_id, weights_s = weights_s, intercept_s = intercept_s,
       weights_m = weights_m, intercept_m = intercept_m, weights_p = weights_p)
}


# a sequential run for debug
t0=Sys.time()
res <- foreach(j = seq_len(nrow(protein_ann)), .errorhandling = "pass",
               .packages = c("dplyr", "data.table")) %do% process_protein(j)
t1=Sys.time()
t1-t0


#
# # parallel run
# library(foreach)
# library(doSNOW)
# Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
# chr_key <- ifelse(is.na(protein_ann$chr), "missing", as.character(protein_ann$chr)) # Group all proteins by chromosome - for fast parallel computing
# chr_groups <- split(seq_len(nrow(protein_ann)), chr_key)
# n_cores <- max(1L, parallel::detectCores() - 2L)
# cl <- parallel::makeCluster(n_cores)
# parallel::clusterSetRNGStream(cl, iseed = 1334L)
# doSNOW::registerDoSNOW(cl)
#
# pb <- txtProgressBar(min = 0, max = length(chr_groups), style = 3)
# progress_fun <- function(n) setTxtProgressBar(pb, n)
# snow_opts <- list(progress = progress_fun, preschedule = FALSE)
#
# t0 <- Sys.time()
# res <- foreach(
#   chr_idx = unname(chr_groups),
#   .combine = c,
#   .multicombine = TRUE,
#   .errorhandling = "pass",
#   .packages = c("dplyr", "data.table", "ctOWAS", "MiXcan"),
#   .options.snow = snow_opts
# ) %dopar% {
#
#   ans <- lapply(chr_idx, function(j) {
#     process_protein(j)
#   })
#
#   chr_value <- protein_ann$chr[chr_idx[1]]
#
#   if (!is.na(chr_value)) {
#     cache_key <- as.character(chr_value)
#
#     if (exists(cache_key, envir = genotype_chr_cache, inherits = FALSE)) {
#       rm(list = cache_key, envir = genotype_chr_cache)
#     }
#   }
#
#   gc(verbose = FALSE)
#   ans
# }
#
# close(pb)
# parallel::stopCluster(cl)
# t1=Sys.time()
# t1-t0
#
#
#







# ------------------------------------------------------------------------------
# Step 4. Save heart protein weights and skipped-protein log
# ------------------------------------------------------------------------------

# Bind the per-protein result list, write it under the shared naming scheme, and
# report where it went.

weights_ctOWAS_GTEx_heart    <- lapply(res, `[[`, "weights_s")
intercept_ctOWAS_GTEx_heart  <- lapply(res, `[[`, "intercept_s")
weights_MiXcan_GTEx_heart     <- lapply(res, `[[`, "weights_m")
intercept_MiXcan_GTEx_heart   <- lapply(res, `[[`, "intercept_m")
weights_PrediXcan_GTEx_heart  <- lapply(res, `[[`, "weights_p")

gene_names <- vapply(res, `[[`, character(1), "gene_id")
names(weights_ctOWAS_GTEx_heart) <- names(intercept_ctOWAS_GTEx_heart) <-
  names(weights_MiXcan_GTEx_heart) <- names(intercept_MiXcan_GTEx_heart) <-
  names(weights_PrediXcan_GTEx_heart) <- gene_names

weights_ctOWAS_GTEx_heart   <- data.table::rbindlist(lapply(res, `[[`, "weights_s"), fill = TRUE)
intercept_ctOWAS_GTEx_heart <- data.table::rbindlist(lapply(res, `[[`, "intercept_s"), fill = TRUE)
weights_MiXcan_GTEx_heart    <- data.table::rbindlist(lapply(res, `[[`, "weights_m"), fill = TRUE)
intercept_MiXcan_GTEx_heart  <- data.table::rbindlist(lapply(res, `[[`, "intercept_m"), fill = TRUE)
weights_PrediXcan_GTEx_heart <- data.table::rbindlist(lapply(res, `[[`, "weights_p"), fill = TRUE)

save(weights_ctOWAS_GTEx_heart, intercept_ctOWAS_GTEx_heart,
     weights_MiXcan_GTEx_heart, intercept_MiXcan_GTEx_heart,
     weights_PrediXcan_GTEx_heart, file = file.path(results_dir, "GTEx_heart_pi2_ThreeModel_Summary.RData"))

