rm(list = ls())
library(dplyr)
library(ctOWAS)
library(MiXcan)
library(purrr)
library(data.table)

paper_dir <- '/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS'
results_dir <- file.path(paper_dir, 'Results')
drive_pheno_dir <- file.path(paper_dir, 'Data', 'DRIVE')
drive_geno_dir <- file.path(paper_dir, '../Paper_MiXcan/Data/OncoArray')
cache_dir <- file.path(results_dir, 'DRIVE_cached_inputs')
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

run_preparation <- FALSE       # FALSE for subsequent ctOWAS-only runs
run_ctowas <- TRUE
chromosomes <- 1:22
mc_cores <- 1L                # macOS/Linux; use 1 on Windows
ctowas_tag <- 'default'       # Change to preserve results from different runs
ctowas_reg_scale <- 0.01      # NULL uses the package default
mixcan_reg_scale <- 0.01      # ctOWAS with MiXcan weights

# Preparation-only setup ---------------------------------------------
if (run_preparation) {
  load(file.path(results_dir, 'GTEx_Breast_pi2_ThreeModel_Summary.RData'))
  merge_keys <- c('gene', 'rsid', 'varID', 'ref_allele', 'eff_allele')
  mw_input <- list(
    filtered_weights_s %>% rename(ctOWAS_Cell1 = Cell1, ctOWAS_Cell2 = Cell2),
    filtered_weights_m_new %>% rename(MiXcan_Cell1 = Cell1, MiXcan_Cell2 = Cell2),
    filtered_weights_p %>% rename(PrediXcan_tissue = tissue)
  ) %>%
    purrr::reduce(dplyr::full_join, by = merge_keys) %>%
    mutate(
      across(c(ctOWAS_Cell1, ctOWAS_Cell2, MiXcan_Cell1, MiXcan_Cell2,
               PrediXcan_tissue), ~ coalesce(.x, 0)),
      CHR = as.integer(sub('^chr([0-9]+)_.*$', '\\1', varID)),
      POS = as.integer(sub('^chr[0-9]+_([0-9]+)_.*$', '\\1', varID)),
      ref_allele = toupper(ref_allele),
      eff_allele = toupper(eff_allele)
    ) %>%
    select(gene, rsid, varID, CHR, POS, ref_allele, eff_allele,
           ctOWAS_Cell1, ctOWAS_Cell2, MiXcan_Cell1, MiXcan_Cell2,
           PrediXcan_tissue) %>% arrange(CHR, POS, gene)
  setDT(mw_input)

  drive_pheno <- read.csv(file.path(drive_pheno_dir, 'oncoarray-drive.pheno.csv'))
  drive_pheno$subject_index <- as.character(drive_pheno$subject_index)
  pca_result <- fread(file.path(
    drive_geno_dir, 'oncoarray_principal_components_fastmode_2020-10-17.txt'
  ))
  names(pca_result) <- c('FID', paste0('PC', 1:10))
  covariates <- drive_pheno %>%
    select(subject_index, FID, genotyping_center, age_int) %>%
    inner_join(pca_result, by = 'FID') %>% select(-FID) %>%
    mutate(age_int = if_else(age_int == 888, NA_real_, as.numeric(age_int))) %>%
    filter(PC1 <= 0.006) %>% tidyr::drop_na()
  stopifnot(!anyDuplicated(covariates$subject_index),
            !anyDuplicated(drive_pheno$subject_index))
}

run_chr <- function(chr) {
  message('\nPreparing chromosome ', chr)
  model_chr <- mw_input[CHR == chr]
  if (!nrow(model_chr)) stop('No model variants for chromosome ', chr)

  # STEP 1: Read only variants required by the prediction models.
  genotype_path <- file.path(
    drive_geno_dir, sprintf('oncoarray_dosages_chr%02d.txt.gz', chr)
  )
  key_path <- tempfile(fileext = '.keys')
  subset_file <- tempfile(fileext = '.txt')
  on.exit(unlink(c(key_path, subset_file)), add = TRUE)
  variant_keys <- unique(model_chr[, .(
    CHR_key = paste0('chr', CHR), POS,
    REF = ref_allele, ALT = eff_allele
  )])
  fwrite(variant_keys, key_path, sep = '\t', col.names = FALSE)
  awk_script <- 'BEGIN { FS = "\t" }
    FILENAME != "-" { k[$1 FS $2 FS $3 FS $4] = 1; next }
    FNR == 1 { print; next }
    match($0, /^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*/) {
      split(substr($0, 1, RLENGTH), a, FS)
      c = a[1]; sub(/^chr/, "", c); c = "chr" c
      if ((c FS a[2] FS toupper(a[4]) FS toupper(a[5])) in k) print
    }'
  cmd <- sprintf('gzip -dc %s | awk %s %s - > %s',
                 shQuote(genotype_path), shQuote(awk_script),
                 shQuote(key_path), shQuote(subset_file))
  status <- system2('bash', c('-o', 'pipefail', '-c', shQuote(cmd)))
  if (status != 0L) stop('Genotype filtering failed for chromosome ', chr)
  drive_genome <- fread(subset_file, header = TRUE, check.names = FALSE)
  setnames(drive_genome, 1L, 'CHR')
  drive_genome[, `:=`(
    CHR = as.integer(sub('^chr', '', CHR)), POS = as.integer(POS),
    REF = toupper(REF), ALT = toupper(ALT)
  )]
  if (!nrow(drive_genome)) stop('No matching genotype variants: chromosome ', chr)

  # Assumes the original genotype file has nine metadata columns.
  SampleID <- intersect(names(drive_genome)[-(1:9)], covariates$subject_index)
  D <- drive_pheno$affection_status[match(SampleID, drive_pheno$subject_index)]
  # Use the same observed outcomes for GWAS and all association tests.
  keep <- !is.na(D)
  SampleID <- SampleID[keep]
  D <- D[keep]
  if (!length(D) || !all(D %in% c(0, 1)) || length(unique(D)) != 2L) {
    stop('affection_status must contain both 0 and 1 among retained samples.')
  }
  cov <- as.data.frame(covariates)[
    match(SampleID, covariates$subject_index),
    setdiff(names(covariates), 'subject_index'), drop = FALSE
  ]
  cov_mat <- model.matrix(~ ., data = cov)[, -1, drop = FALSE]
  stopifnot(nrow(cov_mat) == length(D))
  geno_cols <- c(names(drive_genome)[1:9], SampleID)
  drive_genome2 <- drive_genome[, ..geno_cols]
  mw_drive_input <- merge(
    model_chr, drive_genome2,
    by.x = c('CHR', 'POS', 'ref_allele', 'eff_allele'),
    by.y = c('CHR', 'POS', 'REF', 'ALT'), all = FALSE
  )
  if (!nrow(mw_drive_input)) stop('No merged variants: chromosome ', chr)
  X_genome <- t(as.matrix(mw_drive_input[, ..SampleID]))

  # STEP 2: GWAS (same method as the original script).
  drive_gwas_result <- as.data.frame(run_gwas(
    X = X_genome, D = D, family0 = stats::binomial(),
    covar = cov_mat, method_binomial = 'score'
  ))
  stopifnot(nrow(drive_gwas_result) == nrow(mw_drive_input),
            'Z' %in% names(drive_gwas_result))
  fwrite(
    cbind(mw_drive_input[, .(CHR, POS, ref_allele, eff_allele)],
          drive_gwas_result),
    file.path(results_dir, sprintf('drive_gwas_chr%02d.tsv', chr)), sep = '\t'
  )

  # STEP 3: Cache inputs before running the other association models.
  # Cache contains individual-level genotypes; compress=FALSE favors speed.
  drive_input_total <- cbind(drive_gwas_result, as.data.frame(mw_drive_input))
  saveRDS(list(
    input = drive_input_total, SampleID = SampleID,
    n0 = sum(D == 0), n1 = sum(D == 1)
  ), file.path(cache_dir, sprintf('DRIVE_input_chr%02d.rds', chr)),
  compress = FALSE)

  # STEP 4: Original MiXcan and PrediXcan analyses.
  drive_split_df <- split(drive_input_total, drive_input_total$gene)
  drive_gene <- names(drive_split_df)
  drive_result <- data.frame(
    Gene_ID = drive_gene, p_m_1 = NA_real_, p_m_2 = NA_real_,
    p_m = NA_real_, p_predixcan = NA_real_
  )
  for (g in seq_along(drive_gene)) {
    gene <- drive_gene[g]
    gene_df <- drive_split_df[[gene]]
    message('Processing gene: ', gene)
    new_x <- t(as.matrix(gene_df[, SampleID, drop = FALSE]))

    # Run MiXcan
    W <- as.matrix(gene_df[, c('MiXcan_Cell1', 'MiXcan_Cell2'), drop = FALSE])
    if (any(W != 0)) {
      pred <- MiXcan::MiXcan_prediction(weight = W, new_x = new_x)
      mixcan_res <- MiXcan_association(
        new_y = pred, new_outcome = D, new_cov = cov_mat, family = 'binomial'
      )
      drive_result[g, c('p_m_1', 'p_m_2', 'p_m')] <- c(
        mixcan_res$cell1_p, mixcan_res$cell2_p, mixcan_res$p_combined
      )
    }
    # Run PrediXcan
    W <- as.matrix(gene_df[, 'PrediXcan_tissue', drop = FALSE])
    if (any(W != 0)) {
      pred <- new_x %*% W
      PrediXcan_res <- glm(D ~ pred + cov_mat, family = 'binomial')
      drive_result[g, 'p_predixcan'] <- summary(PrediXcan_res)$coefficients[2, 4]
    }
  }
  fwrite(drive_result, file.path(
    results_dir, sprintf('DRIVE_other_models_chr%02d.csv', chr)
  ))
  drive_result
}

# Separate ctOWAS analysis -------------------------------------------
run_ctowas_chr <- function(chr, tag = 'default', reg_scale = NULL,
                           mixcan_reg_scale = 0.01) {
  message('\nRunning ctOWAS: chromosome ', chr)
  cached <- readRDS(file.path(cache_dir, sprintf('DRIVE_input_chr%02d.rds', chr)))
  input <- as.data.frame(cached$input)
  SampleID <- cached$SampleID
  drive_split_df <- split(input, input$gene)
  drive_gene <- names(drive_split_df)
  drive_result <- data.frame(
    Gene_ID = drive_gene,
    z_o_1 = NA_real_, z_o_2 = NA_real_,
    p_o_1 = NA_real_, p_o_2 = NA_real_, p_o = NA_real_,
    p_s_1 = NA_real_, p_s_2 = NA_real_, p_s = NA_real_
  )
  for (g in seq_along(drive_gene)) {
    gene_df <- drive_split_df[[drive_gene[g]]]
    message('Processing gene: ', drive_gene[g])
    new_x <- t(as.matrix(gene_df[, SampleID, drop = FALSE]))
    args <- list(gwas_z_score = gene_df$Z, x_g = new_x,
                 n0 = cached$n0, n1 = cached$n1, family = 'binomial')

    # ctOWAS with its own prediction weights
    W <- as.matrix(gene_df[, c('ctOWAS_Cell1', 'ctOWAS_Cell2'), drop = FALSE])
    if (any(W != 0)) {
      own_args <- c(list(W = W), args)
      if (!is.null(reg_scale)) own_args$reg_scale <- reg_scale
      fit <- do.call(ctOWAS_assoc_test_K, own_args)
      drive_result[g, c('z_o_1', 'z_o_2', 'p_o_1', 'p_o_2', 'p_o')] <-
        c(fit$Z_join, fit$p_join_vec, fit$p_join)
    }
    # ctOWAS with MiXcan prediction weights
    W <- as.matrix(gene_df[, c('MiXcan_Cell1', 'MiXcan_Cell2'), drop = FALSE])
    if (any(W != 0)) {
      mix_args <- c(list(W = W), args)
      mix_args$reg_scale <- mixcan_reg_scale
      fit <- do.call(ctOWAS_assoc_test_K, mix_args)
      drive_result[g, c('p_s_1', 'p_s_2', 'p_s')] <-
        c(fit$p_join_vec, fit$p_join)
    }
  }
  fwrite(drive_result, file.path(
    results_dir, sprintf('DRIVE_ctOWAS_%s_chr%02d.csv', tag, chr)
  ))

  # Also produce the original combined column layout, when available.
  other_path <- file.path(results_dir, sprintf('DRIVE_other_models_chr%02d.csv', chr))
  if (file.exists(other_path)) {
    combined <- merge(drive_result, fread(other_path), by = 'Gene_ID', all = TRUE)
    fwrite(combined, file.path(
      results_dir, sprintf('DRIVE_four_models_%s_chr%02d.csv', tag, chr)
    ))
  }
  drive_result
}

# Parallel execution; propagate worker errors to the main session.
run_chromosomes <- function(FUN, ...) {
  if (.Platform$OS.type == 'windows' || mc_cores == 1L) {
    out <- lapply(chromosomes, FUN, ...)
  } else {
    out <- parallel::mclapply(chromosomes, FUN, ..., mc.cores = mc_cores)
  }
  failed <- vapply(out, inherits, logical(1), what = 'try-error')
  if (any(failed)) stop('Failed chromosomes: ', paste(chromosomes[failed], collapse = ', '))
  out
}

if (run_preparation) res_other <- run_chromosomes(run_chr)
if (run_ctowas) {
  res_ctowas <- run_chromosomes(
    run_ctowas_chr, tag = ctowas_tag,
    reg_scale = ctowas_reg_scale, mixcan_reg_scale = mixcan_reg_scale
  )
}

# Subsequent runs: set run_preparation <- FALSE near the top, and change
# ctowas_tag, ctowas_reg_scale, and/or mixcan_reg_scale as desired.
# Alternatively, rerun just this call in the current R session:
# res_ctowas_005 <- run_chromosomes(
#   run_ctowas_chr, tag = 'reg005', reg_scale = 0.05, mixcan_reg_scale = 0.05
# )


# ------ Harvest results -------------------

library(data.table)

# Match the tag used in run_ctowas_chr()
ctowas_tag <- "default"
chromosomes <- 1:22

# Read and combine chromosome-specific results
combine_chr <- function(pattern) {
  files <- file.path(results_dir, sprintf(pattern, chromosomes))

  if (any(!file.exists(files))) {
    stop("Missing files:\n", paste(files[!file.exists(files)], collapse = "\n"))
  }

  rbindlist(
    lapply(seq_along(files), function(i) {
      dat <- fread(files[i])
      dat[, CHR := chromosomes[i]]
      dat
    }),
    use.names = TRUE
  )
}

# MiXcan and PrediXcan
results_other <- combine_chr("DRIVE_other_models_chr%02d.csv")

# ctOWAS and ctOWAS with MiXcan weights
results_ctowas <- combine_chr(
  paste0("DRIVE_ctOWAS_", ctowas_tag, "_chr%02d.csv")
)

# Check uniqueness before merging
keys <- c("CHR", "Gene_ID")

# Combine all models
results_all <- merge(results_ctowas,results_other,by = keys, all = TRUE)

setorder(results_all, CHR, Gene_ID)

# Save combined results
fwrite(results_all,  file.path(results_dir, paste0("DRIVE_all_models_", ctowas_tag, "_all_chr.csv")))

dim(results_all)
head(results_all)
