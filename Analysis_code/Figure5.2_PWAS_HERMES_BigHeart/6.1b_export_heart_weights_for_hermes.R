library(data.table)

paper_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
training_dir <- file.path(paper_dir, "Results", "heart_protein_training")
model_file <- file.path(training_dir, "GTEx_heart_pi2_ThreeModel_Summary.RData")


loaded <- load(model_file)
cat("Loaded objects from", model_file, ":\n  ", paste(loaded, collapse = ", "), "\n")

required_objs <- c("weights_ctOWAS_GTEx_heart", "weights_MiXcan_GTEx_heart")
missing_objs <- setdiff(required_objs, loaded)
if (length(missing_objs)) {
  stop("Model file is missing expected objects: ", paste(missing_objs, collapse = ", "),
       call. = FALSE)
}

export_weights <- function(weights_obj, type_label, out_file) {
  dt <- as.data.table(weights_obj)
  setnames(dt, old = "weight_others", new = "weight_other", skip_absent = TRUE)
  required_cols <- c("gene_id", "gene_name", "varID", "chr", "pos",
                      "ref_allele", "eff_allele", "weight_cardiomyocytes", "weight_other")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols)) {
    stop("Exported ", type_label, " weights are missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  dt[, type := type_label]
  fwrite(dt, out_file)
  cat("Wrote", nrow(dt), "rows (", uniqueN(dt$gene_id), "genes ) to", out_file, "\n")
}

export_weights(weights_ctOWAS_GTEx_heart, "ctOWAS",
                file.path(training_dir, "weights_ctOWAS_for_hermes.csv"))
export_weights(weights_MiXcan_GTEx_heart, "MiXcan",
                file.path(training_dir, "weights_MiXcan_for_hermes.csv"))

cat("\nDone. Next: run 7.2_HERMES_prepare_data_pwas.R twice, once with\n",
    "  HERMES_MODEL_TAG=ctOWAS Rscript 7.2_HERMES_prepare_data_pwas.R\n",
    "and once with\n",
    "  HERMES_MODEL_TAG=MiXcan Rscript 7.2_HERMES_prepare_data_pwas.R\n")
