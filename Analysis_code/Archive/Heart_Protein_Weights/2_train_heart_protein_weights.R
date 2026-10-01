# Train 2-component PWAS MiXcan models on GTEx heart proteomics -----------------
#
# Goal:
#   Build protein prediction weights for a downstream PWAS/ctOWAS analysis.
#   Each protein is modeled separately using cis genotype predictors, covariates,
#   and two cell-fraction components:
#     1. Cardiomyocytes
#     2. Other = Fibroblasts + Endothelial + Smooth muscle/pericyte + Immune
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

library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ctOWAS)

source('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/R/ctOWAS_train_k.R')
source('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/Analysis_code/1_Prediction_Evaluation_in_Heart/2_training_helper_function.R')

# set directory
paper_dir <- Sys.getenv( "PAPER_ctOWAS_DIR", unset = '/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
data_dir <- file.path(paper_dir, "Data")
heart_dir <- file.path(data_dir, "Heart")
protein_dir <- file.path(heart_dir, "GTEx_Pi_Estimate")
results_dir <- file.path(paper_dir, "Results", "heart_protein_weights", "training_model_weights")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# find file

protein <- load_protein_data(file.path(protein_dir,"Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData"))
protein_sample_cols <- names(protein)[-(1:5)]
pi_data <- read_pi_data(file.path(protein_dir, "BayesDeBulk_pi.tsv"), protein_sample_cols) %>%
  distinct(DonorID, .keep_all = TRUE)
cov_data <- read_covariates(file.path(paper_dir, "IntermediateResults", "covariate_EA_with_age.txt")) %>%
  distinct(DonorID, .keep_all = TRUE) 
# GWS: Paper_ctOWAS/New generated files/codes/by_chr_nomiss -- unpruned dosage files w/o missingness. 
geno_raw_dir <- Sys.getenv("PWAS_GENO_RAW_DIR", unset = file.path(paper_dir, "IntermediateResults", "WGS", "by_chr_nomiss"))

sample_map <- data.frame(
  RawSampleID = protein_sample_cols,
  SampleID = normalize_sample_id(protein_sample_cols),
  DonorID = donor_id(protein_sample_cols),
  stringsAsFactors = FALSE
) %>% distinct(DonorID, .keep_all = TRUE)



# ------------------------------------------------------------------------------
# Step 2. Align all training inputs by GTEx donor ID
# ------------------------------------------------------------------------------

# Final training sample set: donors with protein abundance, cell fractions, and
# covariates available. The covariate file is already EA-filtered; the separate
# EA keep file is applied only when it has enough overlap with these samples.
sample_map <- sample_map %>%
  inner_join(pi_data, by = "DonorID", suffix = c("", "_pi")) %>%
  inner_join(cov_data, by = "DonorID", suffix = c("", "_cov"))


protein_ann <- annotation_from_protein(protein)
protein_ann <- fill_missing_annotation(protein_ann, file.path(data_dir, "ensembl38.txt"))
protein_ann <- protein_ann %>% distinct(gene_id, .keep_all = TRUE)


# ------------------------------------------------------------------------------
# Step 3. Train one PWAS MiXcan model per protein
# ------------------------------------------------------------------------------

# res_weights_all <- vector("list", nrow(protein_ann))
# skipped <- list()
# diagnostics <- list()

# Log a protein we cannot train and move on. Uses <<- so it also works from
# inside the MiXcan_train_K error handler below.
add_skip <- function(gene_id, reason) {
  skipped[[length(skipped) + 1L]] <<- data.frame(gene_id = gene_id, reason = reason)
}

for (j in seq_len(nrow(protein_ann))) {
  target <- protein_ann[j, ]
  cat("Processing", j, "of", nrow(protein_ann), target$gene_id, "\n")

  # Without gene coordinates we cannot define the cis genotype window.
  if (is.na(target$chr) || is.na(target$start) || is.na(target$end)) {
    add_skip(target$gene_id, "missing_annotation")
    next
  }

  # Load per-protein or chromosome-wide genotype dosages. Proteins without a
  # matching raw file are logged rather than stopping the full run.
  geno_target <- read_genotype_for_target(
    target_id = target$gene_id,
    chr = target$chr,
    start = target$start,
    end = target$end
  )
  if (is.null(geno_target)) {
    add_skip(target$gene_id, "missing_or_empty_genotype_raw")
    next
  }

  y_row <- which(protein$gene_id == target$gene_id)
  if (length(y_row) != 1L) {
    add_skip(target$gene_id, "protein_row_not_unique")
    next
  }

  # Match genotype donors to the aligned protein/pi/covariate donor set.
  donors <- intersect(sample_map$DonorID, geno_target$x$DonorID)
  if (length(donors) < 20) {
    add_skip(target$gene_id, "too_few_genotyped_samples")
    next
  }

  aligned <- sample_map[match(donors, sample_map$DonorID), ]
  geno_rows <- geno_target$x[match(donors, geno_target$x$DonorID), ]

  # Final model matrices:
  #   y    = protein abundance for this protein
  #   x    = cis genotype dosage matrix
  #   z    = covariate matrix
  #   pi_k = two cell-fraction columns: Cardiomyocytes and Other
  y <- as.numeric(protein[y_row, aligned$RawSampleID, drop = TRUE])
  x <- as.matrix(geno_rows[, geno_target$snp_annot$raw_col, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- donors
  colnames(x) <- geno_target$snp_annot$varID
  candidate_snp_num <- ncol(x)

  z <- as.matrix(aligned[, cov_cols, drop = FALSE])
  pi_k <- as.matrix(aligned[, cell_type_cols, drop = FALSE])

  # Remove samples with any missing value in model inputs.
  complete_idx <- complete.cases(y) &
    complete.cases(x) &
    complete.cases(z) &
    complete.cases(pi_k)
  if (sum(complete_idx) < 20) {
    add_skip(target$gene_id, "too_few_complete_samples")
    next
  }

  x <- x[complete_idx, , drop = FALSE]
  y <- y[complete_idx]
  z <- z[complete_idx, , drop = FALSE]
  pi_k <- pi_k[complete_idx, , drop = FALSE]
  complete_sample_num <- length(y)

  # Drop extremely rare/unused dosage columns. This mirrors the Breast_BCAC_2CellTypes
  # training script's mean-dosage filter.
  keep_snps <- colMeans(x, na.rm = TRUE) > 0.05
  if (!any(keep_snps)) {
    add_skip(target$gene_id, "no_snps_after_maf_filter")
    next
  }

  x <- x[, keep_snps, drop = FALSE]
  snp_num_after_maf <- ncol(x)
  snp_annot <- geno_target$snp_annot[match(colnames(x), geno_target$snp_annot$varID), ]
  x_name_matrix <- snp_annot %>%
    transmute(
      varID = varID,
      gene_id = target$gene_id,
      gene_name = target$gene_name,
      chr = chr,
      pos = pos,
      ref_allele = ref_allele,
      eff_allele = eff_allele,
      dosed_allele = dosed_allele
    )

  # Use a deterministic but protein-specific CV split so reruns are reproducible.
  set.seed(1334 + j * 149053)
  foldid <- sample(1:10, length(y), replace = TRUE)

  # Fit the two-component MiXcan model. Errors are logged per protein so one
  # failed fit does not stop training for the remaining proteins.
  fit <- tryCatch(
    ctOWAS::MiXcan_train_K(
      y = y,
      x = x,
      pi_k = pi_k,
      cov = z,
      xNameMatrix = x_name_matrix,
      yName = target$gene_id,
      foldid = foldid,
      alpha = mixcan_alpha,
      lambda_choice = lambda_choice
    ),
    error = function(e) {
      add_skip(target$gene_id, paste0("MiXcan_train_K_failed: ", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    next
  }

  # Save only SNPs with at least one meaningfully nonzero cell-component weight.
  # glmnet can return tiny numerical values around 1e-16, which should not be
  # counted as real selected SNP weights.
  weights <- combine_weights(fit, snp_annot, target, cell_type_cols)
  weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
  weights <- weights[rowSums(abs(weights[, weight_cols, drop = FALSE]) > weight_eps) > 0, ]
  if (nrow(weights)) {
    res_weights_all[[j]] <- weights
  }
  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    gene_id = target$gene_id,
    gene_name = target$gene_name,
    chr = target$chr,
    sample_num = complete_sample_num,
    candidate_snp_num = candidate_snp_num,
    snp_num_after_maf = snp_num_after_maf,
    selected_snp_num = nrow(weights),
    type = fit$type,
    alpha = mixcan_alpha,
    lambda_choice = lambda_choice,
    lambda_cell = fit$lambda_cell,
    lambda_tissue = fit$lambda_tissue,
    cv_mse_cell = fit$cv_mse_cell,
    cv_mse_tissue = fit$cv_mse_tissue,
    cv_r2_cell = fit$cv_r2_cell,
    cv_r2_tissue = fit$cv_r2_tissue
  )

  if (j %% 100 == 0) {
    cat("Processed", j, "proteins\n")
  }
}

# ------------------------------------------------------------------------------
# Step 4. Save heart protein weights and skipped-protein log
# ------------------------------------------------------------------------------

# Bind the per-protein result list, write it under the shared naming scheme, and
# report where it went.
save_output <- function(rows, label, tag) {
  file <- file.path(
    results_dir,
    paste0("weights_heart_protein_cardiomyocytes_other", tag, output_suffix, ".csv")
  )
  write_csv(bind_rows(rows), file)
  cat("Saved ", label, ": ", file, "\n", sep = "")
}

save_output(res_weights_all, "weights", "")
save_output(skipped, "skipped log", "_skipped")
save_output(diagnostics, "diagnostics", "_diagnostics")
