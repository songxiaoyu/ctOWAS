rm(list = ls())
library(dplyr)
library(ctOWAS)
library(data.table)

paper_dir <- '/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS'
results_dir <- file.path(paper_dir, 'Results')
drive_pheno_dir <- file.path(paper_dir, 'Data', 'DRIVE')
drive_geno_dir <- file.path(paper_dir, '../Paper_MiXcan/Data/OncoArray')
cache_dir <- file.path(results_dir, 'DRIVE_3ct_cached_inputs')
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

run_preparation <- F       # FALSE for subsequent ctOWAS-only runs
run_ctowas <- TRUE
chromosomes <- 22:1
mc_cores <- 1L                # macOS/Linux; use 1 on Windows
ctowas_tag <- 'default'       # Change to preserve results from different runs
ctowas_reg_scale <- 0.01      # NULL uses the package default

# Preparation-only setup ---------------------------------------------
if (run_preparation) {
  weights_env <- new.env(parent = emptyenv())
  load(file.path(results_dir, 'GTEx_Breast_pi3_Comparison.RData'),
       envir = weights_env)
  if (!exists('filtered_weights_s_3ct', envir = weights_env, inherits = FALSE)) {
    stop('The weights file must contain filtered_weights_s_3ct.')
  }
  mw_input <- copy(as.data.table(weights_env$filtered_weights_s_3ct))
  required <- c('gene', 'varID', 'ref_allele', 'eff_allele',
                'Cell1', 'Cell2', 'Cell3')
  stopifnot(all(required %in% names(mw_input)))
  mw_input[, `:=`(
    CHR = suppressWarnings(as.integer(sub(
      '^(?:chr)?([^_]+)_.*$', '\\1', varID, perl = TRUE))),
    POS = suppressWarnings(as.integer(sub(
      '^(?:chr)?[^_]+_([0-9]+)_.*$', '\\1', varID, perl = TRUE))),
    ref_allele = toupper(ref_allele),
    eff_allele = toupper(eff_allele)
  )]
  if (anyNA(mw_input[, .(CHR, POS, Cell1, Cell2, Cell3)])) {
    stop('Missing/invalid variant coordinates or prediction weights.')
  }
  setorder(mw_input, CHR, POS, gene)

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
    file.path(results_dir, sprintf('drive_gwas_3ct_chr%02d.tsv', chr)), sep = '\t'
  )

  # STEP 3: Cache inputs before running the other association models.
  # Cache contains individual-level genotypes; compress=FALSE favors speed.
  drive_input_total <- cbind(drive_gwas_result, as.data.frame(mw_drive_input))
  saveRDS(list(
    input = drive_input_total, SampleID = SampleID,
    n0 = sum(D == 0), n1 = sum(D == 1)
  ), file.path(cache_dir, sprintf('DRIVE_ctOWAS_3ct_input_chr%02d.rds', chr)),
  compress = FALSE)

  invisible(chr)
}

# Separate ctOWAS analysis -------------------------------------------
run_ctowas_chr <- function(chr, tag = 'default', reg_scale = NULL) {
  message('\nRunning ctOWAS: chromosome ', chr)
  cached <- readRDS(file.path(cache_dir, sprintf('DRIVE_ctOWAS_3ct_input_chr%02d.rds', chr)))
  input <- as.data.frame(cached$input)
  SampleID <- cached$SampleID
  drive_split_df <- split(input, input$gene)
  drive_gene <- names(drive_split_df)
  drive_result <- data.frame(
    Gene_ID = drive_gene,
    z_1 = NA_real_, z_2 = NA_real_, z_3 = NA_real_,
    p_1 = NA_real_, p_2 = NA_real_, p_3 = NA_real_, p_joint = NA_real_
  )
  for (g in seq_along(drive_gene)) {
    gene_df <- drive_split_df[[drive_gene[g]]]
    message('Processing gene: ', drive_gene[g])
    new_x <- t(as.matrix(gene_df[, SampleID, drop = FALSE]))
    args <- list(gwas_z_score = gene_df$Z, x_g = new_x,
                 n0 = cached$n0, n1 = cached$n1, family = 'binomial')

    # ctOWAS with its own prediction weights
    W <- as.matrix(gene_df[, c('Cell1', 'Cell2', 'Cell3'), drop = FALSE])
    if (any(W != 0)) {
      own_args <- c(list(W = W), args)
      if (!is.null(reg_scale)) own_args$reg_scale <- reg_scale
      fit <- do.call(ctOWAS_assoc_test_K, own_args)
      drive_result[g, c('z_1', 'z_2', 'z_3', 'p_1', 'p_2', 'p_3', 'p_joint')] <-
        c(fit$Z_join, fit$p_join_vec, fit$p_join)
    }

  }
  fwrite(drive_result, file.path(
    results_dir, sprintf('DRIVE_ctOWAS_3ct_%s_chr%02d.csv', tag, chr)
  ))

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

if (run_preparation) res_input <- run_chromosomes(run_chr)
if (run_ctowas) {
  res_ctowas <- run_chromosomes(
    run_ctowas_chr, tag = ctowas_tag,
    reg_scale = ctowas_reg_scale
  )
}

# Merge all chromosome results; require every requested chromosome.
combine_ctowas_results <- function(tag = 'default', chr = chromosomes) {
  files <- file.path(results_dir,
    sprintf('DRIVE_ctOWAS_3ct_%s_chr%02d.csv', tag, chr))
  if (any(!file.exists(files))) {
    stop('Missing chromosome files:\n',
         paste(files[!file.exists(files)], collapse = '\n'))
  }
  result <- rbindlist(lapply(seq_along(files), function(i) {
    dat <- fread(files[i])
    dat[, CHR := chr[i]]
    dat
  }), use.names = TRUE)
  setcolorder(result, c('Gene_ID', 'CHR', 'z_1', 'z_2', 'z_3',
                       'p_1', 'p_2', 'p_3', 'p_joint'))
  setorder(result, CHR, Gene_ID)
  fwrite(result, file.path(results_dir,
    sprintf('DRIVE_ctOWAS_3ct_%s_all_chr.csv', tag)))
  result
}

if (run_ctowas) {
  drive_result_all <- combine_ctowas_results(ctowas_tag)
  if (interactive()) View(drive_result_all)
}



# ------ Harvest three-cell-type ctOWAS results -----------------------

library(data.table)

ctowas_tag <- "default"
chromosomes <- 1:22

files <- file.path(
  results_dir,
  sprintf("DRIVE_ctOWAS_3ct_%s_chr%02d.csv", ctowas_tag, chromosomes)
)

if (any(!file.exists(files))) {
  stop("Missing files:\n", paste(files[!file.exists(files)], collapse = "\n"))
}

results_all_3ct <- rbindlist(
  lapply(seq_along(files), function(i) {
    dat <- fread(files[i])
    dat[, CHR := chromosomes[i]]
    dat
  }),
  use.names = TRUE
)

# Check uniqueness
stopifnot(
  !anyDuplicated(results_all_3ct[, .(CHR, Gene_ID)])
)

setcolorder(
  results_all_3ct,
  c("Gene_ID", "CHR", "z_1", "z_2", "z_3",
    "p_1", "p_2", "p_3", "p_joint")
)

setorder(results_all_3ct, CHR, Gene_ID)

fwrite(
  results_all_3ct,
  file.path(
    results_dir,
    sprintf("DRIVE_ctOWAS_3ct_%s_all_chr.csv", ctowas_tag)
  )
)

dim(results_all_3ct)
head(results_all_3ct)
