#!/usr/bin/env Rscript

library(data.table)
setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')

# load data
paper_dir <- getwd()
out_dir <-file.path(paper_dir,"Results/heart_protein_prediction_correlation_HLV351")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("Github/Analysis_code/Figure6_GTEx_Heart_Proteomics_Training_and_Evaluation/helper_function.R")

# load real expr
expr <- fread(file.path(paper_dir, "Data/GTEx_Heart/gene_tpm_2022-06-06_v10_heart_left_ventricle.gct"), skip = 2, check.names = FALSE)

# Load predictions
predictions <- data.table::fread(file.path(paper_dir,
                                           "Results/predicted_protein_HLV351_with_intercepts",
                                           "predicted_heart_protein_HLV351_all_chr_three_models.csv"))

ctOWAS <- predictions[model == "ctOWAS"]
MiXcan <- predictions[model == "MiXcan"]

predixcan <- predictions[model == "PrediXcan",
                         .(sample_id, gene_id, gene_name, pred_tissue_predixcan = pred_tissue)]
# clean merge and compare data

extract_donor_id <- function(x) {
  vapply(
    strsplit(as.character(x), "-", fixed = TRUE),
    function(z) paste(utils::head(z, 2L), collapse = "-"),
    character(1)
  )
}


pick_dosage_column <- function(rsid, dosed, raw_names) {
  candidate <- paste0(rsid, "_", dosed)

  if (!is.na(candidate) && candidate %in% raw_names) {
    candidate
  } else {
    NA_character_
  }
}

predictions[, donor_id := extract_donor_id(sample_id)]

ctOWAS <- predictions[
  model == "ctOWAS",
  .(
    donor_id,
    gene_id,
    gene_name,
    method = "ctOWAS",
    predicted = pred_tissue
  )
]

MiXcan <- predictions[
  model == "MiXcan",
  .(
    donor_id,
    gene_id,
    gene_name,
    method = "MiXcan",
    predicted = pred_tissue
  )
]

predixcan <- predictions[
  model == "PrediXcan",
  .(
    donor_id,
    gene_id,
    gene_name,
    method = "PrediXcan",
    predicted = pred_tissue
  )
]

pred <- rbindlist(
  list(ctOWAS, MiXcan, predixcan),
  use.names = TRUE
)

pred <- pred[
  !is.na(donor_id) &
    nzchar(donor_id) &
    !is.na(gene_id) &
    is.finite(predicted)
]

# GTEx TPM gene IDs and protein-annotation gene IDs both carry a GENCODE/
# Ensembl version suffix (e.g. "ENSG00000141510.16"), but the two annotation
# sources can pin different releases, so matching on the raw versioned ID can
# silently drop genes whose version string differs even though it's the same
# gene. Match on the version-stripped ID instead, and carry it through to the
# expression side so the merge below joins on it rather than the versioned ID.
pred[, gene_id_novers := sub("\\..*$", "", as.character(gene_id))]
target_genes <- unique(pred$gene_id_novers)
target_donors <- unique(pred$donor_id)
expr[, gene_id_novers := sub("\\..*$", "", as.character(Name))]
expr <- expr[gene_id_novers %in% target_genes]
sample_cols <- setdiff(names(expr), c("Name", "Description"))
sample_donors <- extract_donor_id(sample_cols)
keep_sample_cols <- sample_cols[sample_donors %in% target_donors]

if (!length(keep_sample_cols)) {
  stop("No expression samples overlap prediction donors.", call. = FALSE)
}

expr_long <- melt(
  expr[, c("Name", "Description", "gene_id_novers", keep_sample_cols), with = FALSE],
  id.vars = c("Name", "Description", "gene_id_novers"),
  variable.name = "expression_sample_id",
  value.name = "tpm"
)
setnames(expr_long, c("Name", "Description"), c("expression_gene_id", "expression_gene_name"))
expr_long[, donor_id := extract_donor_id(expression_sample_id)]
expr_long[, tpm := as.numeric(tpm)]

# If a donor has more than one RNA sample, average expression at donor level.
expr_donor <- expr_long[, .(
  expression_gene_name = expression_gene_name[1],
  tpm = mean(tpm, na.rm = TRUE)
), by = .(gene_id_novers, donor_id)]
expr_donor[, log2_tpm_plus1 := log2(tpm + 1)]

dat <- merge(
  pred,
  expr_donor,
  by = c("gene_id_novers", "donor_id"),
  all = FALSE,
  sort = FALSE
)

cor_by_gene <- dat[, .(
  n = .N,
  pearson_tpm = safe_cor(predicted, tpm, "pearson"),
  spearman_tpm = safe_cor(predicted, tpm, "spearman"),
  pearson_log2_tpm_plus1 = safe_cor(predicted, log2_tpm_plus1, "pearson"),
  spearman_log2_tpm_plus1 = safe_cor(predicted, log2_tpm_plus1, "spearman")
), by = .(method, gene_id, gene_name)]

wide <- dcast(
  cor_by_gene,
  gene_id + gene_name ~ method,
  value.var = c("n", "pearson_tpm", "spearman_tpm",
                "pearson_log2_tpm_plus1", "spearman_log2_tpm_plus1")
)
if (all(c("pearson_log2_tpm_plus1_ctOWAS",
          "pearson_log2_tpm_plus1_PrediXcan") %in% names(wide))) {
  wide[, delta_pearson_log2_tpm_plus1 :=
         pearson_log2_tpm_plus1_ctOWAS -
         pearson_log2_tpm_plus1_PrediXcan]
}

cor_file <- file.path(out_dir, "prediction_expression_correlation_by_gene.csv")
wide_file <- file.path(out_dir, "prediction_expression_correlation_by_gene_wide.csv")
plot_file <- file.path(out_dir, "prediction_expression_correlation_boxplot.png")
summary_file <- file.path(out_dir, "prediction_expression_correlation_summary.txt")

fwrite(cor_by_gene, cor_file)
fwrite(wide, wide_file)

plot_dt <- cor_by_gene[!is.na(pearson_log2_tpm_plus1)]
png(plot_file, width = 1800, height = 1400, res = 180)
boxplot(
  pearson_log2_tpm_plus1 ~ method,
  data = plot_dt,
  col = c("#4C78A8", "#F58518"),
  ylab = "Pearson correlation with log2(TPM + 1)",
  xlab = "",
  main = "Predicted Protein vs Real Gene Expression"
)
stripchart(
  pearson_log2_tpm_plus1 ~ method,
  data = plot_dt,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  col = rgb(0, 0, 0, 0.18),
  add = TRUE
)
abline(h = 0, lty = 2, col = "gray40")
dev.off()

summary_dt <- cor_by_gene[, .(
  genes = .N,
  genes_with_nonmissing_pearson_log = sum(!is.na(pearson_log2_tpm_plus1)),
  median_pearson_log2_tpm_plus1 = median(pearson_log2_tpm_plus1, na.rm = TRUE),
  mean_pearson_log2_tpm_plus1 = mean(pearson_log2_tpm_plus1, na.rm = TRUE),
  median_spearman_log2_tpm_plus1 = median(spearman_log2_tpm_plus1, na.rm = TRUE),
  mean_spearman_log2_tpm_plus1 = mean(spearman_log2_tpm_plus1, na.rm = TRUE),
  median_pearson_tpm = median(pearson_tpm, na.rm = TRUE),
  mean_pearson_tpm = mean(pearson_tpm, na.rm = TRUE)
), by = method]

if ("delta_pearson_log2_tpm_plus1" %in% names(wide)) {
  paired <- wide[!is.na(delta_pearson_log2_tpm_plus1)]
  delta_lines <- c(
    paste0("Paired genes: ", nrow(paired)),
    paste0("Median ctOWAS - PrediXcan Pearson log2(TPM+1): ",
           median(paired$delta_pearson_log2_tpm_plus1, na.rm = TRUE)),
    paste0("Mean ctOWAS - PrediXcan Pearson log2(TPM+1): ",
           mean(paired$delta_pearson_log2_tpm_plus1, na.rm = TRUE)),
    paste0("Genes where ctOWAS > PrediXcan: ",
           sum(paired$delta_pearson_log2_tpm_plus1 > 0, na.rm = TRUE)),
    paste0("Genes where PrediXcan > ctOWAS: ",
           sum(paired$delta_pearson_log2_tpm_plus1 < 0, na.rm = TRUE))
  )
} else {
  delta_lines <- "Paired comparison unavailable."
}

summary_lines <- c(
  paste0("Correlation analysis date: ", Sys.time()),
  paste0("Overlapping donors: ", uniqueN(dat$donor_id)),
  paste0("Overlapping genes: ", uniqueN(dat$gene_id)),
  paste0("Rows used after merge: ", nrow(dat)),
  "",
  "Method summaries:",
  capture.output(print(summary_dt)),
  "",
  "Paired comparison:",
  delta_lines,
  "",
  paste0("Saved per-gene correlations: ", cor_file),
  paste0("Saved wide comparison table: ", wide_file),
  paste0("Saved boxplot: ", plot_file)
)
writeLines(summary_lines, summary_file)
