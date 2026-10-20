rm(list=ls())

library(xCell2)
library(BiocParallel)
library(tidyverse)
library(janitor)
library(ggplot2)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir=getwd()
data_dir <- file.path(paper_dir, "Data")
results_dir <- file.path(paper_dir, "Results")

load(file.path(results_dir, "GTEx_clean.RData"))

# known markers -----------
markers <- list( Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1", "CLDN3", "CLDN4",     ## Broad/luminal epithelial
                 "KRT5", "KRT14", "KRT15", "KRT17", "TP63"   ), ## Basal-myoepithelial
  Adipocyte = c("ADIPOQ", "PLIN1", "PLIN4", "FABP4", "LPL", "LEP", "CIDEC", "CIDEA", "LIPE", "DGAT2", "PNPLA2", "G0S2"), ## Primarily mature white adipocytes
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "DPT", "COL14A1", "FBLN1", "FBLN2", "PCOLCE", "PCOLCE2", "WISP2") )    ## Collagen/ECM fibroblast program

lengths(markers)
# Epithelial  Adipocyte Fibroblast
#         13         12         12
save(markers, file = file.path(results_dir, "markers.RData"))

#  marker ensemble data --------------
load(file.path(results_dir, "markers.RData"))
markers_ensembl <- lapply(markers, \(x) na.omit(ensembl38$gene_stable_id[match(x, ensembl38$gene_name)]) |>as.character())
lengths(markers_ensembl)
# Epithelial  Adipocyte Fibroblast
# 13         12         12
save(markers_ensembl, file = file.path(results_dir, "markers_ensembl.RData"))

# load(file.path(results_dir, "markers_ensembl.RData"))



#  xcell2 score --------------
make_xcell2_signatures <- function(marker_list, available_genes, n_signatures = 3L, fraction = 0.8, min_size = 8L) {

  marker_list <- lapply( marker_list,intersect,y = available_genes)

  if (any(lengths(marker_list) < min_size)) {
    stop("Fewer than ", min_size, " markers available for: ",
      paste(names(marker_list)[lengths(marker_list) < min_size], collapse = ", "))}

  signatures <- list()

  for (cell_type in names(marker_list)) {

    genes <- marker_list[[cell_type]]
    n_genes <- length(genes)

    signature_size <- max(min_size,floor(fraction * n_genes))
    signature_size <- min(signature_size, n_genes)

    for (j in seq_len(n_signatures)) {

      start <- floor( (j - 1L) * n_genes / n_signatures)
      index <- (start + seq_len(signature_size) - 1L) %% n_genes + 1L
      signature_name <- paste0(cell_type, "#custom",j)
      signatures[[signature_name]] <- genes[index]
    }
  }
  signatures
}

custom_signatures <- make_xcell2_signatures(marker_list = markers, available_genes = rownames(expr2), n_signatures = 3, fraction = 0.8, min_size = 8)
cell_names <- unique(sub("#.*$", "", names(custom_signatures)))
genes_used <- unique(unlist(custom_signatures, use.names = FALSE))
spill_identity <- diag(length(cell_names))
dimnames(spill_identity) <- list(cell_names, cell_names)

## Placeholder parameters are unused with rawScores = TRUE
custom_params <- data.frame(
  celltype = cell_names,
  a = rep(1, length(cell_names)),
  b = rep(1, length(cell_names)),
  m = rep(1, length(cell_names)),
  n = rep(0, length(cell_names)),
  stringsAsFactors = FALSE
)

custom_xcell2 <- methods::new(
  "xCell2Object",
  signatures = custom_signatures,
  dependencies = list(),
  params = custom_params,
  spill_mat = spill_identity,
  genes_used = genes_used
)

methods::validObject(custom_xcell2)
custom_xcell2


xCellScore<- xCell2::xCell2Analysis(
  mix = expr2,
  xcell2object = custom_xcell2,
  minSharedGenes = 0.8,
  rawScores = TRUE,
  spillover = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
save(xCellScore, file = file.path(results_dir, "xCellScore.RData"))

load(file.path(results_dir, "xCellScore.RData"))
dim(xCellScore)







# ----------------------------------------------------------------------
# -------------------- Plot the results ----------
# ----------------------------------------------------------------------


cor_by_marker <- function(marker, proportion) {
  apply(expr2[intersect(marker, rownames(expr2)), ], 1,
        \(x) cor(x, proportion, use = "pairwise.complete.obs"))}

plot_df <- bind_rows(
  # Two cell types: epithelial vs stromal (adi + fibr)
  tibble(correlation = cor_by_marker(markers$Epithelial, pi2_old$mean_trim_0.05),
         cell_type = "Epithelial", method = "MiXcan", panel = "Two Cell Types"),
  tibble(correlation = c(
    cor_by_marker(markers$Adipocyte, 1 - pi2_old$mean_trim_0.05),
    cor_by_marker(markers$Fibroblast, 1 - pi2_old$mean_trim_0.05)),
    cell_type = "Stromal", method = "MiXcan", panel = "Two Cell Types"),
  tibble(correlation = cor_by_marker(markers$Epithelial, pi2_new[, 1]),
         cell_type = "Epithelial", method = "ctOWAS", panel = "Two Cell Types"),
  tibble(correlation = c(
    cor_by_marker(markers$Adipocyte, pi2_new[, 2]),
    cor_by_marker(markers$Fibroblast, pi2_new[, 2])),
    cell_type = "Stromal", method = "ctOWAS", panel = "Two Cell Types"),

  # Three cell types: ctOWAS
  tibble(correlation = cor_by_marker(markers$Epithelial, pi_once$Ecomposition[, 1]),
         cell_type = "Epithelial", method = "ctOWAS", panel = "Three Cell Types"),
  tibble(correlation = cor_by_marker(markers$Adipocyte, pi_once$Ecomposition[, 2]),
         cell_type = "Adipocyte", method = "ctOWAS", panel = "Three Cell Types"),
  tibble(correlation = cor_by_marker(markers$Fibroblast, pi_once$Ecomposition[, 3]),
         cell_type = "Fibroblast", method = "ctOWAS", panel = "Three Cell Types")
) %>%
  filter(is.finite(correlation))

plot_df <- plot_df %>%
  mutate(
    panel = factor(panel, levels = c("Two Cell Types", "Three Cell Types")),
    cell_type = factor(
      cell_type,
      levels = c("Epithelial", "Stromal", "Adipocyte", "Fibroblast")
    ),
    method = factor(method, levels = c("MiXcan", "ctOWAS"))
  )

ggplot(plot_df, aes(cell_type, correlation, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.7) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(MiXcan = "#4C78A8", ctOWAS = "#CC79A7")) +
  labs(x = NULL, y = "Marker gene–proportion correlation", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


# two cell type plot
two_ct_df <- plot_df %>%
  filter(panel == "Two Cell Types") %>%
  mutate(
    cell_type = factor(cell_type,
                       levels = c("Epithelial", "Stromal")),
    method = factor(method, levels = c("MiXcan", "ctOWAS"))
  )

ggplot(two_ct_df, aes(cell_type, correlation, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.7) +
  scale_fill_manual(values = c(MiXcan = "#4C78A8", ctOWAS = "#CC79A7")) +
  labs(
    x = NULL,
    y = "Marker Gene–Proportion Correlation",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
