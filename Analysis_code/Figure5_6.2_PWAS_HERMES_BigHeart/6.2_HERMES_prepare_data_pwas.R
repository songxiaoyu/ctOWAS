# Prepare HERMES DCM GWAS for PWAS/ctOWAS analysis.
#
# Inputs:
#   1. Heart protein weights from Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R
#   2. HERMES GWAS lifted from hg19 to hg38 by 1_liftover_hermes_gwas.py
#
# Output:
#   1. Per-chromosome SNP ID lists for 1000Genome LD extraction
#   2. Per-chromosome merged weight + GWAS RDS files for ctOWAS

library(data.table)
library(dplyr)
library(readr)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()

# Which trained model to prepare: "ctOWAS" or "MiXcan" (produced by
# 7.1b_export_heart_weights_for_hermes.R). Run this script once per tag.
model_tag <- Sys.getenv("HERMES_MODEL_TAG", "ctOWAS")
if (!model_tag %in% c("ctOWAS", "MiXcan")) {
  stop("HERMES_MODEL_TAG must be 'ctOWAS' or 'MiXcan', got: ", model_tag, call. = FALSE)
}
cat("Preparing HERMES input for model:", model_tag, "")

workspace_dir <- file.path(paper_dir, "Results", "hermes_pwas",paste0("hermes_workspace_", model_tag))

# load HERMES
gwas_file <- file.path(paper_dir, "Data",  "HERMES_GWAS_DCM_EUR","FORMAT-METAL_Pheno5_EUR_hg38_rsid.tsv.gz")

# load weights (from 7.1b_export_heart_weights_for_hermes.R)
weights_dir <- file.path(paper_dir, "Results", "heart_protein_training")
weights_file <- file.path(weights_dir, paste0("weights_", model_tag, "_for_hermes.csv"))
if (!file.exists(weights_file)) {
  stop("Missing weights file: ", weights_file,
       " -- run 7.1b_export_heart_weights_for_hermes.R first.", call. = FALSE)
}
mw_input <- fread(weights_file)

# create new folders
hermes_dir <- file.path(paper_dir,"Results", "hermes_pwas")
input_dir <- file.path(workspace_dir, "hermes_input")
filtered_id_dir <- file.path(workspace_dir, "hermes_filtered_id")
dir.create(hermes_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtered_id_dir, recursive = TRUE, showWarnings = FALSE)


normalize_chr <- function(x) {sub("^chr", "", as.character(x))}
normalize_varid <- function(x) {x <- as.character(x);x <- sub("^chr", "", x);
  x <- sub("_b38$", "", x); gsub("_", ":", x)}

# Step 1. Read trained PWAS model weights.

required_weight_cols <- c("gene_id", "gene_name", "varID", "chr", "pos",
  "ref_allele", "eff_allele", "weight_cardiomyocytes", "weight_other", "type")
missing_weight_cols <- setdiff(required_weight_cols, names(mw_input))
if (length(missing_weight_cols)) {
  stop("Weights file is missing columns: ", paste(missing_weight_cols, collapse = ", "))
}

mw_input <- mw_input %>%
  mutate(
    CHR = paste0("chr", normalize_chr(chr)),
    varID = normalize_varid(varID)
  )
setDT(mw_input)
setkey(mw_input, CHR, varID)

# Step 2. Read lifted HERMES GWAS once and standardize its columns.
gwas <- fread(gwas_file)
required_gwas_cols <- c("CHR38", "POS38", "A1", "A2", "A1_beta", "se", "pval")
missing_gwas_cols <- setdiff(required_gwas_cols, names(gwas))
if (length(missing_gwas_cols)) {
  stop("HERMES GWAS is missing columns: ", paste(missing_gwas_cols, collapse = ", "))
}
if (!"rsid" %in% names(gwas)) {
  gwas[, rsid := NA_character_]
}

# HERMES A1_beta is the beta for A1. Use A1 as Effect.Gwas and A2 as baseline.
gwas <- gwas %>%
  mutate(
    CHR = paste0("chr", normalize_chr(CHR38)),
    POS = as.integer(POS38),
    Baseline.Gwas = as.character(A2),
    Effect.Gwas = as.character(A1),
    beta.Gwas = as.numeric(A1_beta),
    SE.Gwas = as.numeric(se),
    pval.Gwas = as.numeric(pval),
    rsid = as.character(rsid)
  ) %>%
  filter(!is.na(POS), !is.na(beta.Gwas), !is.na(SE.Gwas))
setDT(gwas)

# Step 3. For each chromosome, match forward and reverse allele orientation.
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")

  chr_label <- paste0("chr", chr)
  gwas_chr <- gwas[CHR == chr_label]
  if (nrow(gwas_chr) == 0) {
    cat("  --> No GWAS rows, skipping\n")
    next
  }

  gwas_fwd <- copy(gwas_chr)
  gwas_fwd[, varID := paste0(chr, ":", POS, ":", Baseline.Gwas, ":", Effect.Gwas)]
  gwas_fwd[, flip := FALSE]

  gwas_rev <- copy(gwas_chr)
  gwas_rev[, varID := paste0(chr, ":", POS, ":", Effect.Gwas, ":", Baseline.Gwas)]
  gwas_rev[, flip := TRUE]

  gwas_long <- rbindlist(list(gwas_fwd, gwas_rev), use.names = TRUE, fill = TRUE)
  mw_gwas_input <- mw_input[CHR == chr_label][gwas_long, on = .(varID), nomatch = 0L]

  if (nrow(mw_gwas_input) == 0) {
    cat("  --> No overlapping weighted GWAS SNPs\n")
  } else {
    mw_gwas_input[flip == TRUE, beta.Gwas := -beta.Gwas]
    mw_gwas_input[flip == TRUE, c("Baseline.Gwas", "Effect.Gwas") := list(Effect.Gwas, Baseline.Gwas)]
  }

  write.table(
    data.frame(mw_gwas_input$varID),
    file = file.path(filtered_id_dir, sprintf("hermes_filtered_chr%d_gwas_id_pwas.txt", chr)),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(
    mw_gwas_input,
    file = file.path(input_dir, sprintf("chr%d_mw_gwas_input_hermes_pwas.rds", chr))
  )
}
