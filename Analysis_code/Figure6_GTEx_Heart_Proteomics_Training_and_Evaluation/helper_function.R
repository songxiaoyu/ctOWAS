# ------------------------------------------------------------------------------
# Step 2. Define helper functions
# ------------------------------------------------------------------------------
get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}


# Convert GTEx sample IDs from R-safe column names back to dash-separated IDs.
# Example: GTEX.QESD.0526.SM.DEJHB -> GTEX-QESD-0526-SM-DEJHB.
normalize_sample_id <- function(x) {gsub("\\.", "-", x)}

# Reduce a GTEx sample ID to the donor ID used for matching genotype, covariate, protein, and cell-fraction records.
# Example: GTEX-QESD-0526-SM-DEJHB -> GTEX-QESD.
donor_id <- function(x) {
  x <- normalize_sample_id(x)
  vapply(strsplit(x, "-"), function(parts) {paste(parts[seq_len(min(2L, length(parts)))], collapse = "-")}, character(1))
}


safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  x <- x[keep]
  y <- y[keep]
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  stats::cor(x, y, method = method)
}

# Verify that a required input file exists before the program proceeds
required_file <- function(path, label) {
  if (!file.exists(path)) {stop(label, " not found: ", path, call. = FALSE)}
  path
}

# Covariates are stored as one row per covariate and one column per sample.
# MiXcan needs one row per sample, so transpose the table and add donor IDs.
read_covariates <- function(path) {
  cov_wide <- fread(required_file(path, "Covariate file"), check.names = FALSE)
  cov_name_col <- names(cov_wide)[1]

  cov <- as.data.frame(t(as.data.frame(cov_wide[, -1, with = FALSE])))
  colnames(cov) <- cov_wide[[cov_name_col]]
  cov$SampleID <- rownames(cov)
  cov$DonorID <- donor_id(cov$SampleID)
  rownames(cov) <- NULL

  value_cols <- setdiff(names(cov), c("SampleID", "DonorID"))
  cov[value_cols] <- lapply(cov[value_cols], function(x) as.numeric(as.character(x)))
  cov
}

# The imputed protein file contains the training phenotypes. The original
# non-imputed protein file keeps more proteins but has missing values; this
# script uses Imputed_protein because model fitting expects complete y values.
load_protein_data <- function(path) {
  e <- new.env()
  load(required_file(path, "Protein RData file"), envir = e)
  if (!exists("Imputed_protein", envir = e)) {
    stop("Expected object Imputed_protein in ", path, call. = FALSE)
  }
  get("Imputed_protein", envir = e)
}

# Load five BayesDeBulk fractions, then collapse them to the requested 2-cell
# model: Cardiomyocytes versus all other heart cell types combined.
# The pi table has no sample ID column, so rows are matched to protein sample
# columns by order, following the pi-estimation workflow.
read_pi_data <- function(path, protein_sample_cols) {
  pi <- fread(required_file(path, "BayesDeBulk pi file"))
  source_cell_type_cols <- c(
    "Cardiomyocytes",
    "Fibroblasts",
    "Endothelial",
    "Smooth_muscle_cells_Pericyte",
    "Immune_cells"
  )
  missing_cols <- setdiff(source_cell_type_cols, names(pi))
  if (length(missing_cols)) {
    stop(
      "BayesDeBulk pi file is missing cell-type columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(pi) != length(protein_sample_cols)) {
    stop(
      "BayesDeBulk pi rows (", nrow(pi),
      ") do not match protein sample columns (", length(protein_sample_cols), ").",
      call. = FALSE
    )
  }

  pi$SampleID <- normalize_sample_id(protein_sample_cols)
  pi$DonorID <- donor_id(pi$SampleID)
  pi[, (source_cell_type_cols) := lapply(.SD, as.numeric), .SDcols = source_cell_type_cols]
  pi[, Other := Fibroblasts + Endothelial + Smooth_muscle_cells_Pericyte + Immune_cells]
  pi[, c("SampleID", "DonorID", "Cardiomyocytes", "Other"), with = FALSE]
}

# Keep only EA samples
read_ea_donors <- function(path) {
  keep <- fread(required_file(path, "EA donor keep file"), header = FALSE)
  unique(as.character(keep[[1]]))
}

# Keep the protein annotation needed for cis-window SNP selection and output.
annotation_from_protein <- function(protein) {
  ann <- protein[, c("gene_name", "gene_id", "chr", "start", "end")]
  names(ann) <- c("gene_name", "gene_id", "chr", "start", "end")
  ann$gene_id_no_version <- sub("\\..*$", "", ann$gene_id)
  ann
}

# Fill missing chromosome/start/end fields from ensembl38.txt when possible.
# Most records should already have annotation in the protein file; this is a
# fallback to reduce skipped proteins.
fill_missing_annotation <- function(ann, ensembl_file) {
  if (!file.exists(ensembl_file)) {
    return(ann)
  }

  ens <- fread(ensembl_file)
  setnames(
    ens,
    old = names(ens),
    new = make_clean_names(names(ens))
  )
  needed_cols <- c("gene_stable_id", "chromosome_scaffold_name",
                   "transcript_start_bp", "transcript_end_bp")
  if (!all(needed_cols %in% names(ens))) {
    return(ann)
  }

  ens <- ens[, ..needed_cols] %>%
    group_by(gene_stable_id) %>%
    summarise(
      chromosome_scaffold_name = first(chromosome_scaffold_name),
      transcript_start_bp = min(transcript_start_bp, na.rm = TRUE),
      transcript_end_bp = max(transcript_end_bp, na.rm = TRUE),
      .groups = "drop"
    )
  ann <- ann %>%
    left_join(
      ens,
      by = c("gene_id_no_version" = "gene_stable_id")
    ) %>%
    mutate(
      chr = ifelse(is.na(chr), chromosome_scaffold_name, chr),
      start = ifelse(is.na(start), transcript_start_bp, start),
      end = ifelse(is.na(end), transcript_end_bp, end)
    ) %>%
    dplyr::select(gene_name, gene_id, gene_id_no_version, chr, start, end)
  ann
}

make_clean_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# Genotype input is chromosome-wide PLINK raw dosage:
#   chr<chr>_dosage_nomiss.raw
# The script filters each chromosome file to cis SNPs using the protein gene
# coordinates below.
candidate_raw_files <- function(target_id, chr) {
  parent_raw_dir <- dirname(geno_raw_dir)
  primary_matches <- Sys.glob(file.path(geno_raw_dir, sprintf("chr%s_dosage*.raw", chr)))
  search_dirs <- unique(c(
    geno_raw_dir,
    file.path(geno_raw_dir, "pruned_by_chr"),
    file.path(parent_raw_dir, "pruned_by_chr")
  ))
  raw_names <- c(
    sprintf("chr%s_dosage_nomiss.raw", chr),
    sprintf("chr%s_dosage.raw", chr),
    sprintf("chr%s_dosage_nomiss_LDpruned_500kb_1_r2_0.8.raw", chr)
  )
  unique(c(primary_matches, unlist(lapply(search_dirs, file.path, raw_names), use.names = FALSE)))
}

# PLINK --export A may produce columns like:
#   1:11260293:G:A_G
# or, when variant IDs are rsIDs:
#   rs10904045_C
# The first format carries its own coordinates. The second format needs the
# matching chromosome pvar file to recover chr/pos/ref/alt annotation.
candidate_pvar_files <- function(chr) {
  parent_raw_dir <- dirname(geno_raw_dir)
  unique(c(
    file.path(geno_raw_dir, sprintf("GTEx_EA_chr%s_nomiss.pvar", chr)),
    file.path(geno_raw_dir, sprintf("chr%s_hg38.pvar", chr)),
    file.path(parent_raw_dir, "pruned_by_chr", sprintf("GTEx_EA_chr%s_nomiss.pvar", chr)),
    file.path(parent_raw_dir, "pruned_by_chr", sprintf("chr%s_hg38.pvar", chr))
  ))
}

read_pvar_annotation <- function(chr) {
  pvar_file <- candidate_pvar_files(chr)
  pvar_file <- pvar_file[file.exists(pvar_file)][1]
  if (is.na(pvar_file)) {
    return(NULL)
  }

  pvar <- fread(pvar_file, skip = "#CHROM", check.names = FALSE)
  names(pvar) <- sub("^#", "", names(pvar))
  needed <- c("CHROM", "POS", "ID", "REF", "ALT")
  if (!all(needed %in% names(pvar))) {
    return(NULL)
  }

  pvar[, c("CHROM", "POS", "ID", "REF", "ALT"), with = FALSE] %>%
    rename(
      chr = CHROM,
      pos = POS,
      rsid = ID,
      ref_allele = REF,
      alt_allele = ALT
    ) %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      rsid = as.character(rsid),
      ref_allele = as.character(ref_allele),
      alt_allele = as.character(alt_allele)
    )
}

parse_position_variant_cols <- function(raw_cols) {
  var_id <- sub("_[^_]+$", "", raw_cols)
  dosed_allele <- sub("^.*_", "", raw_cols)
  parts <- strsplit(var_id, ":", fixed = TRUE)

  data.frame(
    varID = var_id,
    raw_col = raw_cols,
    chr = vapply(parts, function(x) x[1], character(1)),
    pos = as.integer(vapply(parts, function(x) x[2], character(1))),
    ref_allele = vapply(parts, function(x) x[3], character(1)),
    eff_allele = vapply(parts, function(x) x[4], character(1)),
    dosed_allele = dosed_allele,
    stringsAsFactors = FALSE
  )
}

read_rsid_position_annotation <- function(chr) {
  if (!nzchar(rsid_annotation_file) || !file.exists(rsid_annotation_file)) {
    return(NULL)
  }

  rsid_pos <- fread(rsid_annotation_file, check.names = FALSE)
  names(rsid_pos) <- make_clean_names(names(rsid_pos))
  rsid_col <- intersect(names(rsid_pos), c("rsid", "rs_id", "id"))[1]
  chr_col <- intersect(names(rsid_pos), c("chr38", "chrom38", "chromosome38", "chr", "chromosome"))[1]
  pos_col <- intersect(names(rsid_pos), c("pos38", "position38", "bp38", "pos", "position", "bp"))[1]
  ref_col <- intersect(names(rsid_pos), c("ref", "ref_allele", "reference_allele"))[1]
  alt_col <- intersect(names(rsid_pos), c("alt", "alt_allele", "effect_allele"))[1]
  dosage_col <- intersect(names(rsid_pos), c("dosage_col", "raw_col", "dosage_column"))[1]
  if (is.na(rsid_col) || is.na(chr_col) || is.na(pos_col)) {
    stop(
      "PWAS_RSID_ANNOT_FILE must contain rsid plus chromosome and position columns.",
      call. = FALSE
    )
  }

  out <- rsid_pos %>%
    transmute(
      rsid = as.character(.data[[rsid_col]]),
      chr = sub("^chr", "", as.character(.data[[chr_col]])),
      pos = as.integer(.data[[pos_col]]),
      ref_allele = if (!is.na(ref_col)) as.character(.data[[ref_col]]) else NA_character_,
      alt_allele = if (!is.na(alt_col)) as.character(.data[[alt_col]]) else NA_character_,
      annotation_raw_col = if (!is.na(dosage_col)) as.character(.data[[dosage_col]]) else NA_character_
    ) %>%
    filter(chr == sub("^chr", "", as.character(chr))) %>%
    distinct(rsid, .keep_all = TRUE)
  out
}

parse_rsid_variant_cols <- function(raw_cols, chr) {
  pvar <- read_pvar_annotation(chr)

  raw_map <- data.frame(
    rsid = sub("_[^_]+$", "", raw_cols),
    raw_col = raw_cols,
    dosed_allele = sub("^.*_", "", raw_cols),
    stringsAsFactors = FALSE
  )

  # First try direct matching when the pvar ID column is also rsID.
  if (!is.null(pvar)) {
    snp_annot <- raw_map %>%
      inner_join(pvar, by = "rsid")
  } else {
    snp_annot <- data.frame()
  }

  # Most current pvar files use coordinate IDs instead. In that case, use the
  # pruned annotation file in pruned_by_chr. If it has ref/alt, no pvar join is
  # needed; otherwise recover ref/alt from pvar by chr/pos.
  if (nrow(snp_annot) == 0) {
    rsid_pos <- read_rsid_position_annotation(chr)
    if (!is.null(rsid_pos)) {
      snp_annot <- raw_map %>%
        inner_join(rsid_pos, by = "rsid")
      if (
        nrow(snp_annot) > 0 &&
        !is.null(pvar) &&
        any(is.na(snp_annot$ref_allele) | is.na(snp_annot$alt_allele))
      ) {
        snp_annot <- snp_annot %>%
          dplyr::select(-ref_allele, -alt_allele) %>%
          inner_join(
            pvar %>% dplyr::select(chr, pos, ref_allele, alt_allele),
            by = c("chr", "pos")
          )
      }
    }
  }

  if (nrow(snp_annot) == 0) {
    return(data.frame())
  }

  snp_annot %>%
    mutate(
      alt_choices = strsplit(alt_allele, ",", fixed = TRUE),
      eff_allele = vapply(seq_along(alt_choices), function(i) {
        choices <- alt_choices[[i]]
        if (dosed_allele[i] %in% choices) {
          dosed_allele[i]
        } else {
          choices[1]
        }
      }, character(1)),
      allele_match = dosed_allele == ref_allele | dosed_allele == eff_allele,
      varID = paste(chr, pos, ref_allele, eff_allele, sep = ":")
    ) %>%
    arrange(raw_col, desc(allele_match)) %>%
    distinct(raw_col, .keep_all = TRUE) %>%
    dplyr::select(varID, raw_col, chr, pos, ref_allele, eff_allele, dosed_allele) %>%
    as.data.frame()
}

parse_raw_variant_cols <- function(raw_cols, chr) {
  var_id <- sub("_[^_]+$", "", raw_cols)
  if (all(grepl("^[^:]+:[0-9]+:[^:]+:[^:]+$", var_id))) {
    return(parse_position_variant_cols(raw_cols))
  }
  parse_rsid_variant_cols(raw_cols, chr)
}

genotype_chr_cache <- new.env(parent = emptyenv())

# Read one chromosome's raw dosage and annotation once, then reuse it for all
# proteins on that chromosome.
read_genotype_for_chr <- function(target_id, chr) {
  cache_key <- as.character(chr)
  if (exists(cache_key, envir = genotype_chr_cache, inherits = FALSE)) {
    return(get(cache_key, envir = genotype_chr_cache, inherits = FALSE))
  }

  files <- candidate_raw_files(target_id, chr)
  raw_file <- files[file.exists(files)][1]
  if (is.na(raw_file)) {
    assign(cache_key, NULL, envir = genotype_chr_cache)
    return(NULL)
  }

  cat("Reading genotype raw file:", raw_file, "\n")
  raw <- fread(raw_file, check.names = FALSE)
  sample_col <- if ("IID" %in% names(raw)) "IID" else names(raw)[2]
  raw$DonorID <- donor_id(as.character(raw[[sample_col]]))

  variant_cols <- setdiff(
    names(raw),
    c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "DonorID")
  )
  snp_annot <- parse_raw_variant_cols(variant_cols, chr)
  out <- list(raw_file = raw_file, raw = raw, snp_annot = snp_annot)
  assign(cache_key, out, envir = genotype_chr_cache)
  out
}

# Select genotype dosages for one protein's +/- 1 Mb cis window.
read_genotype_for_target <- function(target_id, chr, start, end, window_bp = 1e6) {
  geno_chr <- read_genotype_for_chr(target_id, chr)
  if (is.null(geno_chr)) {
    return(NULL)
  }

  snp_annot <- geno_chr$snp_annot
  cis_start <- max(0, as.integer(start) - window_bp)
  cis_end <- as.integer(end) + window_bp
  keep <- snp_annot$chr == as.character(chr) &
    snp_annot$pos >= cis_start &
    snp_annot$pos <= cis_end
  snp_annot <- snp_annot[keep, , drop = FALSE]

  if (nrow(snp_annot) == 0) {
    return(NULL)
  }

  x <- as.data.frame(geno_chr$raw[, c("DonorID", snp_annot$raw_col), with = FALSE])
  list(raw_file = geno_chr$raw_file, x = x, snp_annot = snp_annot)
}

# Turn MiXcan and ctOWAS training model outputs into the flat weight table expected by downstream
# PWAS/GWAS preparation code.
combine_weights <- function(fit, snp_annot, target, cell_type_cols,
                            model=c("ctOWAS", "MiXcan", "PrediXcan")) {

  if (model==("ctOWAS")) {
    weights <- as.data.frame(fit$W)
    weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
    colnames(weights) <- weight_cols
    weights$varID <- rownames(fit$W)
  }
  if (model==("MiXcan")) {
    weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
    weights <- data.frame(varID=fit$beta.SNP.cell1$varID, fit$beta.SNP.cell1$weight, fit$beta.SNP.cell2$weight)
    colnames(weights)[2:3] <- weight_cols

  }
  if (model=="PrediXcan") {
    weights <-fit$beta.all.models[2: (nrow(fit$beta.SNP.cell1)+1),"Tissue"]
    weights<-data.frame(tissue=weights)
    weight_cols="tissue"
    weights$varID <- fit$beta.SNP.cell1$varID
  }

  out <- snp_annot %>%
    dplyr::select(varID, chr, pos, ref_allele, eff_allele, dosed_allele) %>%
    inner_join(weights, by = "varID") %>%
    mutate(
      gene_id = target$gene_id,
      gene_name = target$gene_name
    ) %>%
    dplyr::select(
      gene_id, gene_name, varID, chr, pos, ref_allele, eff_allele,
      dosed_allele, all_of(weight_cols)
    )
  out
}



extract_intercepts <- function(fit, target, cell_type_cols, model=c("ctOWAS", "MiXcan")) {

  if (model=="ctOWAS") {
    intercept_row <- fit$beta.all.models[1, , drop = FALSE]
  }
  if (model=="MiXcan"){
    intercept_row <- fit$beta.all.models[1, -1, drop = FALSE]
  }
  intercept_names <- paste0("intercept_", make_clean_names(cell_type_cols))

  colnames(intercept_row)=intercept_names
  out <- data.frame(
    gene_id   = target$gene_id,
    gene_name = target$gene_name,
    chr       = target$chr,
    intercept_row,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out
}

intercepts_as_model_terms <- function(intercepts, cell_type_cols) {
  weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
  out <- data.frame(
    gene_id = intercepts$gene_id,
    gene_name = intercepts$gene_name,
    term_type = "intercept",
    varID = "Intercept",
    chr = intercepts$chr,
    pos = NA_integer_,
    ref_allele = NA_character_,
    eff_allele = NA_character_,
    dosed_allele = NA_character_,
    type = intercepts$type,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(cell_type_cols)) {
    out[[weight_cols[k]]] <- intercepts[[paste0("intercept_", make_clean_names(cell_type_cols[k]))]]
  }
  out[, c("gene_id", "gene_name", "term_type", "varID", "chr", "pos",
          "ref_allele", "eff_allele", "dosed_allele", weight_cols, "type")]
}


