# Bridges Figure6's freshly-trained heart protein weights
# (Results/heart_protein_training/GTEx_heart_pi2_ThreeModel_Summary.RData)
# into the flat per-model weights CSV that 7.2_HERMES_prepare_data_pwas.R
# expects as input.
#
# Why this script exists: 6.2_train_heart_protein_weights.R's cell type
# columns are named c("Cardiomyocytes", "Others") (plural "Others"), so the
# RData objects it saves carry columns weight_others / intercept_others
# (plural). 7.2 (and the HERMES association code in 7.4) were written
# against an older convention with a SINGULAR "weight_other" column. Rather
# than touch either of those working scripts, this script renames at
# export time.
#
# This also replaces the old weights_ctOWAS_all.csv / weights_MiXcan_all.csv
# files in this same folder (dated before the 6.2 object-name bugfix): those
# were written by a since-changed version of 6.2 and can no longer be
# reproduced by any script currently in this repo, so they should be treated
# as stale. This script always regenerates fresh CSVs directly from whatever
# RData 6.2 most recently produced.
#
# Run this after 6.2_train_heart_protein_weights.R and before
# 7.2_HERMES_prepare_data_pwas.R.

library(data.table)

paper_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
training_dir <- file.path(paper_dir, "Results", "heart_protein_training")
model_file <- file.path(training_dir, "GTEx_heart_pi2_ThreeModel_Summary.RData")

if (!file.exists(model_file)) {
  stop("Missing trained model file: ", model_file,
       " -- run 6.2_train_heart_protein_weights.R first.", call. = FALSE)
}

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
