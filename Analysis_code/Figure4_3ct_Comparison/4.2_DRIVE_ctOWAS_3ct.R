# DRIVE association using ONLY the three-cell-type ctOWAS prediction model.
# Cell1/Cell2/Cell3 follow the column order of pi3_new in the training script.
# Requires genotype coordinates on the same genome build as the model variants.
# Matches model ref/eff alleles to dosage REF/ALT exactly; does not flip alleles.
library(data.table)
library(ctOWAS)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results")
drive_pheno_dir <- file.path(paper_dir, "Data", "DRIVE")
drive_geno_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_MiXcan/Data/OncoArray"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
weight_cols <- c("Cell1", "Cell2", "Cell3")
p_cols <- c("p_o_1", "p_o_2", "p_o_3", "p_o")
variant_key <- function(chr, pos, ref, alt) paste(chr, pos, toupper(ref), toupper(alt), sep = ":")

# 1. Load three-cell-type ctOWAS weights into an isolated environment.
models <- new.env()
load(file.path(results_dir, "GTEx_Breast_pi3_Comparison.RData"), envir = models)
if (!exists("filtered_weights_s_3ct", models, inherits = FALSE)) stop("Missing filtered_weights_s_3ct.")
w <- as.data.table(models$filtered_weights_s_3ct)
stopifnot(all(c("gene", "varID", "ref_allele", "eff_allele", weight_cols) %in% names(w)))
w <- w[, c("gene", "varID", "ref_allele", "eff_allele", weight_cols), with = FALSE]
if (!all(vapply(w[, ..weight_cols], is.numeric, logical(1)))) stop("Weight columns must be numeric.")
if (any(!is.finite(as.matrix(w[, ..weight_cols])))) stop("Weights contain NA/NaN/Inf; check model output.")
w <- w[rowSums(abs(as.matrix(w[, ..weight_cols]))) > 0]
w[, `:=`(gene = as.character(gene), CHR = suppressWarnings(as.integer(sub("^(?:chr)?([^_]+)_.*$", "\\1", varID, perl = TRUE))),
          POS = suppressWarnings(as.integer(sub("^(?:chr)?[^_]+_([0-9]+)_.*$", "\\1", varID, perl = TRUE))),
          ref_allele = toupper(ref_allele), eff_allele = toupper(eff_allele))]
if (anyNA(w[, .(gene, CHR, POS, ref_allele, eff_allele)])) stop("Missing gene/alleles or unrecognized varID format.")
w <- w[CHR %in% 1:22]
if (!nrow(w)) stop("No nonzero autosomal weights.")
w[, key := variant_key(CHR, POS, ref_allele, eff_allele)]
w <- unique(w)
if (anyDuplicated(w[, .(gene, key)])) stop("Conflicting duplicate gene-variant weights.")
setorder(w, CHR, POS, gene)

# 2. Match phenotype and PCs explicitly, preserving the original PC1 ancestry filter.
pheno <- fread(file.path(drive_pheno_dir, "oncoarray-drive.pheno.csv"), colClasses = c(subject_index = "character", FID = "character"))
pc <- fread(file.path(drive_geno_dir, "oncoarray_principal_components_fastmode_2020-10-17.txt"), colClasses = c(subject_ID = "character"))
stopifnot(all(c("subject_index", "FID", "affection_status", "genotyping_center", "age_int") %in% names(pheno)), ncol(pc) >= 11L)
setnames(pc, "subject_ID", "FID")
setnames(pc, names(pc)[2:11], paste0("PC", 1:10))
if (anyNA(pheno$subject_index) || anyDuplicated(pheno$subject_index) || anyNA(pc$FID) || anyDuplicated(pc$FID)) stop("Sample IDs must be unique and nonmissing.")
pheno <- pheno[, .(subject_index, FID, affection_status, genotyping_center, age_int)]
pheno[, `:=`(age_int = as.numeric(age_int), affection_status = as.numeric(affection_status))]
pheno[age_int == 888, age_int := NA_real_]
pheno <- merge(pheno, pc[, c("FID", paste0("PC", 1:10)), with = FALSE], by = "FID", sort = FALSE)
pheno <- pheno[is.finite(PC1) & PC1 <= 0.006]
if (!all(na.omit(pheno$affection_status) %in% c(0, 1))) stop("affection_status must be 0=control, 1=case; recode explicitly if needed.")
pheno <- pheno[complete.cases(pheno[, c("affection_status", "genotyping_center", "age_int", paste0("PC", 1:10)), with = FALSE])]
pc_cols <- paste0("PC", 1:10)
pheno <- pheno[is.finite(pheno$age_int) & rowSums(!is.finite(as.matrix(pheno[, ..pc_cols]))) == 0]
if (!nrow(pheno)) stop("No eligible subjects after phenotype/covariate filtering.")

# 3. Stream only required variants; original dosage format has CHR/POS/.../REF/ALT in columns 1/2/4/5.
read_chr <- function(chr, keys) {
  infile <- file.path(drive_geno_dir, sprintf("oncoarray_dosages_chr%02d.txt.gz", chr))
  if (!file.exists(infile)) stop("Missing genotype file: ", infile)
  keyfile <- tempfile(fileext = ".keys"); subsetfile <- tempfile(fileext = ".txt")
  on.exit(unlink(c(keyfile, subsetfile)), add = TRUE)
  fwrite(unique(keys[, .(CHR, POS, ref_allele, eff_allele)]), keyfile, sep = "\t", col.names = FALSE)
  awk <- 'NR==FNR {k[$1 FS $2 FS $3 FS $4]=1; next} FNR==1 {print; next} {c=$1; sub(/^chr/, "", c); if ((c FS $2 FS toupper($4) FS toupper($5)) in k) print}'
  cmd <- sprintf("gzip -dc %s | awk %s %s - > %s", shQuote(infile), shQuote(awk), shQuote(keyfile), shQuote(subsetfile))
  if (system2("bash", c("-o", "pipefail", "-c", shQuote(cmd))) != 0L) stop("Genotype extraction failed.")
  con <- file(subsetfile, "r"); first <- readLines(con, n = 2L); close(con)
  if (length(first) < 2L) return(NULL)
  header <- strsplit(trimws(first[1]), "[[:space:]]+")[[1]]
  header[1] <- "CHR"
  g <- fread(subsetfile, header = FALSE, skip = 1L, check.names = FALSE)
  if (ncol(g) != length(header) || anyDuplicated(header)) stop("Invalid genotype header.")
  setnames(g, header)
  if (!all(c("CHR", "POS", "REF", "ALT") %in% names(g))) stop("Expected CHR, POS, REF, ALT dosage columns.")
  g[, `:=`(CHR = as.integer(sub("^chr", "", CHR)), POS = as.integer(POS))]
  g[, key := variant_key(CHR, POS, REF, ALT)]
  if (anyDuplicated(g$key)) stop("Duplicate allele-specific genotype variants.")
  g
}

# 4. Run each SNP GWAS once per chromosome, then test each gene using aligned W/Z/LD inputs.
run_chr <- function(chr) {
  wc <- w[CHR == chr]
  out <- unique(wc[, .(Gene_ID = gene, CHR)])
  out[, `:=`(p_o_1 = NA_real_, p_o_2 = NA_real_, p_o_3 = NA_real_, p_o = NA_real_,
             n0 = NA_integer_, n1 = NA_integer_, n_snps = 0L, status = "pending", error = "")]
  out[, n_model_snps := as.integer(table(wc$gene)[Gene_ID])]
  tryCatch({
    g <- read_chr(chr, wc)
    if (is.null(g)) { out[, status := "no_matching_variants"]; return(out) }
    sample_ids <- intersect(names(g), pheno$subject_index)
    if (length(sample_ids) < 20L) stop("Fewer than 20 matched eligible samples.")
    a <- pheno[match(sample_ids, subject_index)]
    D <- a$affection_status
    n0 <- sum(D == 0); n1 <- sum(D == 1)
    set(out, j = "n0", value = as.integer(n0)); set(out, j = "n1", value = as.integer(n1))
    if (!n0 || !n1) stop("Both cases and controls are required.")
    cov_df <- as.data.frame(a[, c("genotyping_center", "age_int", paste0("PC", 1:10)), with = FALSE])
    cov_df$genotyping_center <- factor(cov_df$genotyping_center)
    cov_df <- cov_df[, vapply(cov_df, function(v) length(unique(v)) > 1L, logical(1)), drop = FALSE]
    cov <- if (ncol(cov_df)) model.matrix(~ ., data = cov_df)[, -1L, drop = FALSE] else matrix(numeric(0), nrow(a), 0L)
    # Remove linear dependencies including dependencies on the intercept.
    design <- cbind(Intercept = 1, cov); q <- qr(design)
    independent <- sort(q$pivot[seq_len(q$rank)]); independent <- independent[independent != 1L]
    cov <- design[, independent, drop = FALSE]
    X <- t(as.matrix(g[, ..sample_ids])); storage.mode(X) <- "double"
    dimnames(X) <- list(sample_ids, g$key)
    # Use a common sample set for GWAS and LD: drop incomplete/constant variants.
    keep <- colSums(!is.finite(X)) == 0L
    keep[keep] <- apply(X[, keep, drop = FALSE], 2, var) > 0
    X <- X[, keep, drop = FALSE]
    if (!ncol(X)) stop("No complete, variable SNPs.")
    if (any(X < 0 | X > 2)) stop("Expected ALT-allele dosages in [0, 2].")
    gw <- as.data.frame(run_gwas(X = X, D = D, family0 = stats::binomial(), covar = cov, method_binomial = "score"))
    if (!"Z" %in% names(gw) || nrow(gw) != ncol(X)) stop("run_gwas must return one row per X column, in input order, with column Z.")
    z <- setNames(gw$Z, colnames(X))
    for (i in seq_len(nrow(out))) {
      gene_id <- out$Gene_ID[i]
      wg <- wc[gene == gene_id & key %in% names(z)[is.finite(z)]]
      out$n_snps[i] <- nrow(wg)
      if (!nrow(wg)) { out$status[i] <- "no_usable_snps"; next }
      W <- as.matrix(wg[, ..weight_cols]); rownames(W) <- wg$key
      x_g <- X[, wg$key, drop = FALSE]
      tryCatch({
        fit <- ctOWAS_assoc_test_K(W = W, gwas_z_score = unname(z[wg$key]), x_g = x_g, n0 = n0, n1 = n1, family = "binomial")
        p <- c(fit$p_join_vec, fit$p_join)
        for (k in seq_along(p_cols)) set(out, i, p_cols[k], p[k])
        out$status[i] <- if (nrow(wg) < out$n_model_snps[i]) "ok_partial_snps" else "ok"
      }, error = function(e) { out$status[i] <<- "association_failed"; out$error[i] <<- conditionMessage(e) })
    }
    out
  }, error = function(e) { out[status == "pending", `:=`(status = "chromosome_failed", error = conditionMessage(e))]; out })
}

# 5. Save per-chromosome results and global BH FDR across genes for each test.
# NA marks failed/unavailable tests; it is never replaced with an observed p-value.
res <- lapply(sort(unique(w$CHR)), function(chr) {
  message("Processing chromosome ", chr)
  out <- run_chr(chr)
  fwrite(out, file.path(results_dir, sprintf("DRIVE_ctOWAS_3ct_chr%02d.csv", chr)))
  gc(verbose = FALSE)
  out
})
drive_result_all <- rbindlist(res, use.names = TRUE)
for (k in seq_along(p_cols)) {
  p <- drive_result_all[[p_cols[k]]]; fdr <- rep(NA_real_, length(p)); ok <- is.finite(p)
  fdr[ok] <- p.adjust(p[ok], method = "BH")
  set(drive_result_all, j = sub("^p_", "fdr_", p_cols[k]), value = fdr)
}
fwrite(drive_result_all, file.path(results_dir, "DRIVE_ctOWAS_3ct_all_chr.csv"))
print(drive_result_all[, .N, by = status])
print(drive_result_all[, lapply(.SD, function(p) sum(p < 0.1, na.rm = TRUE)), .SDcols = sub("^p_", "fdr_", p_cols)])
