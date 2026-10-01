library(data.table)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
source("Github/Analysis_code/Figure6_GTEx_Heart_Proteomics_Training_and_Evaluation/helper_function.R")

# Load data -------------------------------------------------------------------

paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results", "heart_protein_training")
load(model_file <- file.path(results_dir, "GTEx_heart_pi2_ThreeModel_Summary.RData"))
# "weights_ctOWAS_GTEx_heart", "intercept_ctOWAS_GTEx_heart",
#  "weights_MiXcan_GTEx_heart", "intercept_MiXcan_GTEx_heart",
#  "weights_PrediXcan_GTEx_heart")
dosage_dir <- file.path(paper_dir, "IntermediateResults/WGS/by_chr_heart_left_ventricle_351_samples_variant_id")
out_dir <- file.path(paper_dir, "Results", "predicted_protein_HLV351_with_intercepts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Standardize the three prediction models ------------------------------------

weights <- rbindlist(list(
  copy(as.data.table(weights_ctOWAS_GTEx_heart))[, `:=`(model = "ctOWAS", weight_tissue = NA_real_)],
  copy(as.data.table(weights_MiXcan_GTEx_heart))[, `:=`(model = "MiXcan", weight_tissue = NA_real_)],
  copy(as.data.table(weights_PrediXcan_GTEx_heart))[, `:=`(
    model = "PrediXcan", weight_cardiomyocytes = NA_real_,
    weight_others = NA_real_, weight_tissue = tissue
  )][, tissue := NULL]
), use.names = TRUE, fill = TRUE)

intercepts <- rbindlist(list(
  copy(as.data.table(intercept_ctOWAS_GTEx_heart))[
    , `:=`(model = "ctOWAS", intercept_tissue = NA_real_)],
  copy(as.data.table(intercept_MiXcan_GTEx_heart))[
    , `:=`(model = "MiXcan", intercept_tissue = NA_real_)],
  unique(as.data.table(weights_PrediXcan_GTEx_heart)[, .(gene_id, gene_name, chr)])[
    , `:=`(model = "PrediXcan", intercept_cardiomyocytes = NA_real_,
           intercept_others = NA_real_, intercept_tissue = 0)]), use.names = TRUE, fill = TRUE)
weights[, chr := as.character(chr)]
intercepts[, chr := as.character(chr)]


# Predict by chromosome -------------------------------------------------------

prediction_files <- character()
match_reports <- list()

for (chr_number in 1:22) {
  message("===== Predicting chr", chr_number, " =====")

  dosage_file <- file.path(dosage_dir, paste0("chr", chr_number, "_HLV351_dosage_nomiss_variant_id.raw"))
  raw_header <- names(fread(file = dosage_file,nrows = 0L))
  id_cols <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), raw_header)
  snp_header <- setdiff(raw_header, id_cols)
  # Create a varID-to-dosage-column map directly from the PLINK header
  dosage_map <- data.table(dosage_col = snp_header)
  dosage_map[, `:=`(
    varID = sub("_[^_]+$", "", dosage_col),
    raw_dosed_allele = toupper(sub("^.*_", "", dosage_col))
  )]

  intercept_chr <- intercepts[as.integer(chr) == chr_number]
  weights_chr <- copy(weights[as.integer(chr) == chr_number])

  weights_chr[, `:=`(varID = as.character(varID),dosed_allele = toupper(dosed_allele))]

  weights_chr <- merge(weights_chr, dosage_map,
    by = "varID",
    all.x = TRUE,
    sort = FALSE)

  # Compute the match report BEFORE dropping unmatched rows below, so it
  # reflects the true match rate rather than trivially showing 100%.
  match_report_chr <- weights_chr[
    ,
    .(
      chr = chr_number,
      n_weight_rows = .N,
      n_with_dosage_col = sum(!is.na(dosage_col)),
      n_missing_dosage_col = sum(is.na(dosage_col))
    ),
    by = model
  ]

  match_reports[[as.character(chr_number)]] <- match_report_chr

  # SNPs in the trained model that have no matching variant in this
  # chromosome's HLV351 dosage file (not genotyped / filtered out in the
  # prediction cohort) get dosage_col = NA from the merge above. They cannot
  # contribute to the prediction, so drop them here (rather than at the
  # matrix-multiply step, where an NA column name would error out) and log
  # them per chromosome so missingness is auditable.
  missing_weights_chr <- weights_chr[is.na(dosage_col)]
  if (nrow(missing_weights_chr)) {
    fwrite(missing_weights_chr, file.path(out_dir, paste0("missing_weight_snps_chr", chr_number, ".csv")))
  }
  weights_chr <- weights_chr[!is.na(dosage_col)]

  use_weights <- weights_chr
  dosage_cols <- unique(use_weights$dosage_col)

  dosage <- fread(file = dosage_file,select = unique(c(id_cols, dosage_cols)))

  sample_dt <- dosage[, .(FID, IID)]
  sample_dt[, sample_id := IID]
  n_samples <- nrow(sample_dt)

  if (length(dosage_cols)) {
    dosage_mat <- as.matrix(dosage[, ..dosage_cols])
    storage.mode(dosage_mat) <- "numeric"
    colnames(dosage_mat) <- dosage_cols
  } else {
    dosage_mat <- matrix(numeric(), nrow = n_samples, ncol = 0)
  }

  pred_list <- vector("list", nrow(intercept_chr))

  for (i in seq_len(nrow(intercept_chr))) {
    model_name <- intercept_chr$model[i]
    gene <- intercept_chr$gene_id[i]
    gene_weights <- use_weights[model == model_name & gene_id == gene]

    if (model_name %in% c("ctOWAS", "MiXcan")) {
      pred_cardio <- rep(intercept_chr$intercept_cardiomyocytes[i], n_samples)
      pred_other <- rep(intercept_chr$intercept_others[i], n_samples)

      if (nrow(gene_weights)) {
        x <- dosage_mat[, gene_weights$dosage_col, drop = FALSE]
        pred_cardio <- pred_cardio + as.numeric(x %*% gene_weights$weight_cardiomyocytes)
        pred_other <- pred_other + as.numeric(x %*% gene_weights$weight_others)
      }
      pred_tissue = pred_cardio +pred_other

      pred_list[[i]] <- data.table(
        FID = sample_dt$FID, IID = sample_dt$IID, sample_id = sample_dt$sample_id,
        model = model_name, gene_id = gene, gene_name = intercept_chr$gene_name[i], chr = chr_number,
        pred_cardiomyocytes = pred_cardio, pred_other = pred_other, pred_tissue = pred_tissue
      )

    # predixcan
    } else if (model_name == "PrediXcan") {
      pred_tissue <- rep(intercept_chr$intercept_tissue[i], n_samples)

      if (nrow(gene_weights)) {
        x <- dosage_mat[, gene_weights$dosage_col, drop = FALSE]
        pred_tissue <- pred_tissue + as.numeric(x %*% gene_weights$weight_tissue)
      }

      pred_list[[i]] <- data.table(
        FID = sample_dt$FID, IID = sample_dt$IID, sample_id = sample_dt$sample_id,
        model = model_name, gene_id = gene, gene_name = intercept_chr$gene_name[i], chr = chr_number,
        pred_cardiomyocytes = NA_real_, pred_other = NA_real_, pred_tissue = pred_tissue
      )
    }
  }

  pred_chr <- rbindlist(pred_list, use.names = TRUE, fill = TRUE)
  pred_file <- file.path(out_dir, paste0("predicted_heart_protein_HLV351_chr", chr_number, "_three_models.csv"))
  fwrite(pred_chr, pred_file)
  prediction_files <- c(prediction_files, pred_file)

  message("Rows written chr", chr_number, ": ", nrow(pred_chr))
  for (j in seq_len(nrow(match_report_chr))) {
    message("Matched ", match_report_chr$model[j], " weight rows chr", chr_number, ": ",
            match_report_chr$n_with_dosage_col[j], " / ", match_report_chr$n_weight_rows[j])
  }

  rm(dosage, dosage_mat, pred_chr, pred_list)
  gc()
}

# Combine and summarize -------------------------------------------------------

if (!length(prediction_files)) stop("No chromosome prediction files were generated.")

match_report <- rbindlist(match_reports, use.names = TRUE, fill = TRUE)
match_report_file <- file.path(out_dir, "snp_matching_report_by_chr_and_model.csv")
fwrite(match_report, match_report_file)

message("Combining chromosome predictions...")
pred_all <- rbindlist(lapply(prediction_files, fread), use.names = TRUE, fill = TRUE)
combined_file <- file.path(out_dir, "predicted_heart_protein_HLV351_all_chr_three_models.csv")
fwrite(pred_all, combined_file)

model_summary <- pred_all[, .(
  samples = uniqueN(sample_id), genes = uniqueN(gene_id), rows = .N
), by = model]
model_match_summary <- match_report[, .(
  total_weight_rows = sum(n_weight_rows), matched_weight_rows = sum(n_with_dosage_col),
  missing_dosage_column_rows = sum(n_missing_dosage_col)
), by = model]
model_summary <- merge(model_summary, model_match_summary, by = "model", all = TRUE)

summary_file <- file.path(out_dir, "prediction_summary.txt")
summary_lines <- c(
  paste0("Prediction date: ", Sys.time()), paste0("Model file: ", model_file),
  paste0("Dosage directory: ", dosage_dir),
  paste0("Combined prediction file: ", combined_file), "", "Summary by model:",
  paste(capture.output(print(model_summary)), collapse = "\n"), "",
  "Note: PrediXcan predictions use a zero intercept because no separate PrediXcan intercept object was supplied."
)
writeLines(summary_lines, summary_file)

message("Saved combined predictions: ", combined_file)
message("Saved matching report: ", match_report_file)
message("Saved summary: ", summary_file)
