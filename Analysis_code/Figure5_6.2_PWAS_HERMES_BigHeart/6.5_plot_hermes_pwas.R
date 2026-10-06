# Plot the output of the HERMES PWAS association script.
# QQ: joint p-values; Venn: two cell types at BH FDR < 0.1.
# Overlap indicates significance in both cell types, not a test of specificity.

library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggforce)
library(bacon)


paper_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
model_tag <- Sys.getenv("HERMES_MODEL_TAG", "ctOWAS")
workspace_dir <- file.path(paper_dir, "Results", "hermes_pwas", paste0("hermes_workspace_", model_tag))
result_table <- file.path(workspace_dir, "hermes_result", "hermes_result_pwas.csv")
figure_dir <- file.path(paper_dir, "Figure")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fdr_cutoff <- 0.1

if (!file.exists(result_table)) stop("Missing combined HERMES result table: ", result_table, call. = FALSE)
tab <- fread(result_table)
required_cols <- c("gene_id", "gene_name", "p_join", "p_cardiomyocytes", "p_other")
missing_cols <- setdiff(required_cols, names(tab))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
if (!nrow(tab)) stop("The HERMES result table is empty.", call. = FALSE)
if (anyNA(tab$gene_id) || anyDuplicated(tab$gene_id)) stop("Expected one row per nonmissing gene_id.", call. = FALSE)

valid_p <- function(p) is.finite(p) & p >= 0 & p <= 1
adjust_p <- function(p, method = "BH") {
  out <- rep(NA_real_, length(p)); keep <- valid_p(p)
  out[keep] <- p.adjust(p[keep], method = method)
  out
}
tab[, `:=`(fwer_p_joint = adjust_p(p_join, "bonferroni"), fdr_p_joint = adjust_p(p_join),
           fdr_cardiomyocytes = adjust_p(p_cardiomyocytes), fdr_other = adjust_p(p_other))]


make_qq_plot <- function(tab) {
  qq <- copy(tab[valid_p(p_join)])
  if (!nrow(qq)) stop("No valid joint p-values for QQ plot.", call. = FALSE)
  setorder(qq, p_join)
  p <- pmax(qq$p_join, .Machine$double.xmin)
  p_bacon <- pmin(pmax(qq$p_join, .Machine$double.xmin), 1 - .Machine$double.eps)
  z_values <- qnorm(p_bacon, lower.tail = FALSE)
  set.seed(20260918)
  bc <- bacon::bacon(teststatistics = z_values)
  lambda_gc <- as.numeric(bacon::inflation(bc))
  qq[, `:=`(obs = -log10(p), exp = -log10(ppoints(.N)))]
  qq[, gene := fifelse(is.na(gene_name) | gene_name == "", as.character(gene_id), as.character(gene_name))]
  qq[, sig_group := fifelse(fwer_p_joint < 0.05, "FWER < 0.05",
                            fifelse(fdr_p_joint < 0.1, "FDR < 0.1", "Not labelled"))]
  labels <- qq[sig_group != "Not labelled"]
  qq[, plot_label := fifelse(sig_group == "Not labelled", "", gene)]
  ggplot(qq, aes(exp, obs)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
    geom_point(data = labels, aes(color = sig_group), size = 2.4) +
    geom_text_repel(
      data = qq, aes(label = plot_label, color = sig_group),
      size = 3.3, fontface = "bold", seed = 1,
      max.overlaps = Inf, max.iter = 100000, max.time = 5,
      box.padding = 0.6, point.padding = 0.5, point.size = 2.4,
      force = 3, min.segment.length = 0,
      segment.color = "grey50", segment.size = 0.3
    ) +
    annotate("text", x = 0, y = Inf, label = sprintf("lambda[GC] == %.3f", lambda_gc), parse = TRUE, hjust = 0, vjust = 1.3, size = 5) +
    scale_color_manual(values = c("FWER < 0.05" = "red", "FDR < 0.1" = "black")) +
    labs(x = expression(Expected ~ -log[10](p)), y = expression(Observed ~ -log[10](p))) +
    theme_classic(base_size = 12) + theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10))
}

p1 = make_qq_plot(tab)

# Count Number  -----------------------
make_count_venn <- function(tab) {
  # Compare genes with valid tests in both cell types; BH adjustment above uses all valid tests per set.
  eligible <- tab[is.finite(fdr_cardiomyocytes) & is.finite(fdr_other)]
  cm <- eligible$gene_id[eligible$fdr_cardiomyocytes < fdr_cutoff]
  other <- eligible$gene_id[eligible$fdr_other < fdr_cutoff]
  counts <- c(length(setdiff(cm, other)), length(intersect(cm, other)), length(setdiff(other, cm)))
  cat("Cardiomyocyte only / Both / Other only:", counts, "\n")
  cat("Genes excluded from Venn due to missing/invalid cell-type tests:", nrow(tab) - nrow(eligible), "\n")
  circles <- data.frame(x0 = c(-0.8, 0.8), y0 = 0, r = 2, type = c("Cardiomyocyte", "Non-cardiomyocyte"))
  colors <- c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")
  ggplot() +
    geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = type, color = type), alpha = 0.5, linewidth = 0.5) +
    annotate("text", x = c(-1.9, 0, 1.9), y = 0, label = counts, size = 6, fontface = "bold") +
    annotate("text", x = c(-1.5, 1.5), y = 2.35, label = c("Cardiomyocyte", "Non-cardiomyocyte"), size = 3.8, fontface = "bold") +
    coord_fixed(xlim = c(-3, 3), ylim = c(-2.2, 2.7)) +
    scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
    labs(title = paste0("Cell-type associations: FDR < ", fdr_cutoff), subtitle = "BH adjustment within each cell type; circles not to scale") +
    theme_void(base_size = 12) + theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10))
}

p2=make_count_venn(tab)



# Write gene name --------------

# Write gene names: red first, then black, from top to bottom.
tab[, `:=`(fwer_cardiomyocytes = adjust_p(p_cardiomyocytes, "bonferroni"),
           fwer_other = adjust_p(p_other, "bonferroni"))]
tab[, gene_label := fifelse(
  is.na(gene_name) | trimws(gene_name) == "",
  as.character(gene_id), as.character(gene_name)
)]
sig_colors <- c("FWER < 0.05" = "red", "FDR < 0.1 only" = "black")

# Use the same eligible genes as the count Venn.
eligible <- tab[is.finite(fdr_cardiomyocytes) & is.finite(fdr_other)]
card_ids <- eligible[fdr_cardiomyocytes < fdr_cutoff, gene_id]
other_ids <- eligible[fdr_other < fdr_cutoff, gene_id]
card_only <- setdiff(card_ids, other_ids)
other_only <- setdiff(other_ids, card_ids)
shared_ids <- intersect(card_ids, other_ids)

make_gene_labels <- function(ids, region, x, size) {
  d <- copy(tab[gene_id %in% ids])
  if (!nrow(d)) return(NULL)

  d[, `:=`(
    q_label = switch(region, card = fdr_cardiomyocytes, other = fdr_other,
                     shared = pmax(fdr_cardiomyocytes, fdr_other)),
    fwer_label = switch(region, card = fwer_cardiomyocytes, other = fwer_other,
                        shared = pmax(fwer_cardiomyocytes, fwer_other))
  )]
  d[, red := fwer_label < 0.05]
  d[, sig_group := fifelse(red, "FWER < 0.05", "FDR < 0.1 only")]

  # Red first, then black; ascending FDR within each group.
  setorder(d, -red, q_label, gene_label)
  spacing <- if (nrow(d) > 1L) min(0.28, 2.8 / (nrow(d) - 1L)) else 0
  d[, `:=`(x = x, y = ((.N + 1) / 2 - seq_len(.N)) * spacing)]

  geom_text(data = d, aes(x, y, label = gene_label, color = sig_group),
            fontface = "bold", size = size, inherit.aes = FALSE)
}

circles <- data.table(x0 = c(-1.05, 1.05), y0 = 0, r = 2.15,
                      type = c("Cardiomyocyte", "Non-cardiomyocyte"))

plot_gene_venn <- ggplot() +
  geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = type),
              alpha = 0.25, color = "black", linewidth = 0.8) +
  make_gene_labels(card_only, "card", x = -2, size = 3.7) +
  make_gene_labels(other_only, "other", x = 2, size = 3.7) +
  make_gene_labels(shared_ids, "shared", x = 0, size = 3.2) +
  annotate("text", x = -1.75, y = 2.45,
           label = sprintf("Cardiomyocyte\n(p = %d)", length(card_ids)), size = 4, fontface = "bold") +
  annotate("text", x = 1.75, y = 2.45,
           label = sprintf("Non-cardiomyocyte\n(p = %d)", length(other_ids)), size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")) +
  scale_color_manual(values = sig_colors) +
  coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-2.55, 2.9)) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

p3 <- plot_gene_venn

# combine

combined_plot <- plot_grid(p1, p3, ncol = 2, labels = c("A", "B"),
                           label_size = 16, label_fontface = "bold", rel_widths = c(1.05, 0.95))



figure_path <- file.path(figure_dir, "Figure_HERMES_QQ_venn.pdf")
ggsave(figure_path, combined_plot, width = 10.2, height = 5.2)

# ---------Table S--------------------

# Input and annotation.
ensembl38 <- janitor::clean_names(readr::read_csv(file.path(data_dir, "ensembl38.txt"), show_col_types = FALSE))
cols <- c("gene_id", "gene_name", "chr", "input_snp_num", "Z_cardiomyocytes", "p_cardiomyocytes", "Z_other", "p_other", "p_join")
d <- copy(as.data.table(tab)[, ..cols])
stopifnot(!anyNA(d$gene_id), !anyDuplicated(d$gene_id))
setnames(d, c("Z_other", "p_other"), c("Z_non_cardiomyocytes", "p_non_cardiomyocytes"))
p_cols <- c(cardiomyocytes = "p_cardiomyocytes", non_cardiomyocytes = "p_non_cardiomyocytes", tissue = "p_join")
stopifnot(all(vapply(d[, unname(p_cols), with = FALSE], is.numeric, logical(1))))

d[, ensembl_id := sub("\\..*$", "", gene_id)]
a <- copy(as.data.table(ensembl38))
a[, ensembl_id := sub("\\..*$", "", gene_stable_id)]
a <- a[!is.na(karyotype_band) & nzchar(trimws(karyotype_band)) & !is.na(chromosome_scaffold_name),
       .(cytoband = paste(sort(unique(paste0(chromosome_scaffold_name, trimws(karyotype_band)))), collapse = ";")), by = ensembl_id]
d[, cytoband := NA_character_]
d[a, on = "ensembl_id", cytoband := i.cytoband]
d[, ensembl_id := NULL]

# Adjust ALL valid tests before filtering; use strict adjusted-p thresholds throughout.
valid_p <- function(p) is.finite(p) & p >= 0 & p <= 1
adjust_p <- function(p, method) {
  ok <- valid_p(p); ans <- rep(NA_real_, length(p))
  ans[ok] <- p.adjust(p[ok], method = method); ans
}
max_sig <- function(p, sig) if (any(sig)) max(p[sig]) else NA_real_

stats <- rbindlist(lapply(names(p_cols), function(level) {
  p <- d[[p_cols[[level]]]]; n <- sum(valid_p(p))
  fwer <- adjust_p(p, "bonferroni"); fdr <- adjust_p(p, "BH")
  s3 <- !is.na(fwer) & fwer < 0.05; s2 <- !is.na(fdr) & fdr < 0.1
  set(d, j = paste0("FWER_", level), value = fwer)
  set(d, j = paste0("FDR_", level), value = fdr)
  data.table(Level = level, N_tests = n, FWER_raw_p_strict_cutoff = if (n) 0.05/n else NA_real_,
             FWER_largest_observed_p = max_sig(p, s3), FDR_largest_observed_p = max_sig(p, s2),
             N_FWER = sum(s3), N_FDR = sum(s2), N_three_stars = sum(s3), N_two_stars = sum(s2 & !s3))
}))
print(stats, digits = 12)  # N_FDR includes genes marked ***.

# Filter using numeric FDR values, then format for export.
# Reset from numeric results, retaining genes significant in any test.
p_cols <- c(cardiomyocytes = "p_cardiomyocytes", non_cardiomyocytes = "p_non_cardiomyocytes", tissue = "p_join")
fdr_cols <- paste0("FDR_", names(p_cols))
out <- copy(d[rowSums(as.matrix(d[, ..fdr_cols]) < 0.1, na.rm = TRUE) > 0])

for (level in names(p_cols)) {
  p_col <- p_cols[[level]]
  fwer_col <- paste0("FWER_", level); fdr_col <- paste0("FDR_", level)
  stopifnot(all(c(p_col, fwer_col, fdr_col) %in% names(out)), is.numeric(out[[p_col]]))

  p <- out[[p_col]]
  s3 <- !is.na(out[[fwer_col]]) & out[[fwer_col]] < 0.05
  s2 <- !is.na(out[[fdr_col]]) & out[[fdr_col]] < 0.1
  stars <- fifelse(s3, "***", fifelse(s2, "**", ""))

  # Verify against the same numeric source, avoiding stale stats.
  full3 <- !is.na(d[[fwer_col]]) & d[[fwer_col]] < 0.05
  full2 <- !is.na(d[[fdr_col]]) & d[[fdr_col]] < 0.1
  stopifnot(sum(s3) == sum(full3), sum(s2 & !s3) == sum(full2 & !full3))

  txt <- rep(NA_character_, length(p)); ok <- valid_p(p)
  txt[ok] <- paste0(sprintf("%.2e", p[ok]), stars[ok])
  set(out, j = p_col, value = txt)
}

z_cols <- c("Z_cardiomyocytes", "Z_non_cardiomyocytes")
out[, (z_cols) := lapply(.SD, function(z) fifelse(is.finite(z), sprintf("%.2f", z), NA_character_)), .SDcols = z_cols]
out[, (grep("^(FWER|FDR)_", names(out), value = TRUE)) := NULL]
setcolorder(out, c("gene_id", "gene_name", "chr", "cytoband", "input_snp_num",
                   "Z_cardiomyocytes", "p_cardiomyocytes", "Z_non_cardiomyocytes", "p_non_cardiomyocytes", "p_join"))


output_file <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/result/Table_S1.csv"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(out, output_file, na = "NA")
cat("Saved", nrow(out), "genes to:", output_file, "\n")

