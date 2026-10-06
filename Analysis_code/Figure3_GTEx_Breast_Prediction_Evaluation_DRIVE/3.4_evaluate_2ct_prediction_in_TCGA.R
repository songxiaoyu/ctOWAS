rm(list=ls())
library(data.table)
library(ctOWAS)
library(doMC)
library(tidyverse)
library(janitor)
library(ggpubr)
library(RSQLite)
library(xCell)
library(BayesDeBulk)


setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
tcga_dir <- file.path(data_dir, "TCGA")
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")
source('Github/R/deNet_purity.R')
source('Github/R/pi_estimation.R')

# load TCGA data and clean -----------------
tcga_white_id <- read_csv(file.path(tcga_dir, "tcga_race.csv")) %>%filter(race == "white") %>% pull(case_submitter_id)
length(tcga_white_id)
# [1] 7415
# genotype data
geno=readRDS(file.path(tcga_dir, "TCGA_genofile.rds"))
dim(geno)
# [1] 171652    110

# expression data
expr_TCGA = data.table::fread(file.path(tcga_dir, "tcga_quantile_normalized_2020-08-17.txt"))
dim(expr_TCGA)
# [1] 25702   112
expr_TCGA <- expr_TCGA %>%
  rename_with(~ str_remove(.x, "-11[AB]$"))%>%
  dplyr::select(gene_id, any_of(tcga_white_id)) %>%
  as_tibble() %>%
  mutate(gene_id = str_sub(gene_id, 1, 15)) %>%
  rename(gene_stable_id = gene_id)
dim(expr_TCGA)
# [1] 25702    104

# sample ID
SampleID=intersect(colnames(expr_TCGA), colnames(geno))
length(SampleID)
# [1] 102
setdiff(colnames(expr_TCGA), colnames(geno))

# dictionary -- Need to find gene name for cell type composition estimation as cell markers are in gene name.
ensembl38 <-read_csv(file.path(data_dir, "ensembl38.txt")) %>%
  janitor::clean_names() %>%
  dplyr::select(gene_stable_id, gene_name) %>% unique


# # 2. cell type composition estimation (skip unless we want to re-estimate pi) ------------
# # TSNets pi - narrow to the marker genes
#
#
# expr_TCGA2 <- expr_TCGA %>%
#   inner_join(ensembl38 %>% distinct(gene_stable_id, .keep_all = TRUE),by = "gene_stable_id") %>%
#   filter(!is.na(gene_name), gene_name != "") %>%
#   mutate(row_variance = apply(pick(-gene_stable_id, -gene_name),1,var,na.rm = TRUE)) %>%
#   arrange(gene_name, desc(row_variance)) %>%
#   distinct(gene_name, .keep_all = TRUE) %>%
#   select(-gene_stable_id, -row_variance) %>%
#   column_to_rownames("gene_name")
#
# make_xcell2_signatures <- function(marker_list, available_genes, n_signatures = 3L, fraction = 0.8, min_size = 8L) {
#
#   marker_list <- lapply( marker_list,intersect,y = available_genes)
#
#   if (any(lengths(marker_list) < min_size)) {
#     stop("Fewer than ", min_size, " markers available for: ",
#          paste(names(marker_list)[lengths(marker_list) < min_size], collapse = ", "))}
#
#   signatures <- list()
#
#   for (cell_type in names(marker_list)) {
#
#     genes <- marker_list[[cell_type]]
#     n_genes <- length(genes)
#
#     signature_size <- max(min_size,floor(fraction * n_genes))
#     signature_size <- min(signature_size, n_genes)
#
#     for (j in seq_len(n_signatures)) {
#
#       start <- floor( (j - 1L) * n_genes / n_signatures)
#       index <- (start + seq_len(signature_size) - 1L) %% n_genes + 1L
#       signature_name <- paste0(cell_type, "#custom",j)
#       signatures[[signature_name]] <- genes[index]
#     }
#   }
#   signatures
# }
# load(file.path(results_dir, "markers.RData"))
# custom_signatures <- make_xcell2_signatures(marker_list = markers, available_genes = rownames(expr_TCGA2), n_signatures = 3, fraction = 0.8, min_size = 8)
# cell_names <- unique(sub("#.*$", "", names(custom_signatures)))
# genes_used <- unique(unlist(custom_signatures, use.names = FALSE))
# spill_identity <- diag(length(cell_names))
# dimnames(spill_identity) <- list(cell_names, cell_names)
#
# ## Placeholder parameters are unused with rawScores = TRUE
# custom_params <- data.frame(
#   celltype = cell_names,
#   a = rep(1, length(cell_names)),
#   b = rep(1, length(cell_names)),
#   m = rep(1, length(cell_names)),
#   n = rep(0, length(cell_names)),
#   stringsAsFactors = FALSE
# )
#
# custom_xcell2 <- methods::new(
#   "xCell2Object",
#   signatures = custom_signatures,
#   dependencies = list(),
#   params = custom_params,
#   spill_mat = spill_identity,
#   genes_used = genes_used
# )
#
# methods::validObject(custom_xcell2)
# custom_xcell2
#
#
# xCellScore<- xCell2::xCell2Analysis(
#   mix = expr_TCGA2,
#   xcell2object = custom_xcell2,
#   minSharedGenes = 0.8,
#   rawScores = TRUE,
#   spillover = FALSE,
#   BPPARAM = BiocParallel::SerialParam()
# )
# save(xCellScore, file = file.path(results_dir, "xCellScore_TCGA.RData"))
#
# source('Github/R/deNet_composition_K.R')
# cell_type_2 <- list(Epithelial = "Epithelial", Stromal = c("Adipocyte", "Fibroblast"))
# pi_prior <- estimate_prior(xCellScore = xCellScore, cell_types = cell_type_2, mu = c(0.4, 0.6))
#
# apply(pi_prior, 2 ,summary)
# ## work on unique markers
# all_markers <- unlist(markers, use.names = FALSE)
# expr2_all <- expr_TCGA2[all_markers,,drop = FALSE]
# dim(expr2_all)
# marker_list_2 <- list(Epithelial = markers[["Epithelial"]],
#                       Stromal = unique(c(markers[["Adipocyte"]], markers[["Fibroblast"]])))
#
# pi2_TCGA_new <- pi_estimation_K(expr = expr2_all, prior = pi_prior, marker_list=marker_list_2,
#                            n_iteration = 100, seed = 1, verbose = T)
# pi2_TCGA_new |>
#   as.data.frame() |>
#   tibble::rownames_to_column("sample") |>
#   readr::write_csv(file.path(results_dir, "pi_TCGA.csv"))



# 3. load pi2 prediction models -------------------
load('Results/GTEx_Breast_pi2_ThreeModel_Summary.Rdata')
# filtered_weights_m;filtered_weights_p;filtered_weights_s;
# filtered_intercept_m;filtered_intercept_s

pi2 <- read_csv(file.path(results_dir, "pi_TCGA.csv"))
gene=unique(c(filtered_weights_s$gene,filtered_weights_m_new$gene, filtered_weights_m_old$gene,filtered_weights_p$gene ))

cor=  NULL
for (g in  seq_along(gene)){

  yName=gene[g];
  yName_TCGA = sub("\\..*$", "", yName)

  y=expr_TCGA[match(yName_TCGA, expr_TCGA$gene_stable_id), SampleID, drop = FALSE]%>% as.numeric()
  if (mean(is.na(y))>0.5) {next}
  pi <- pi2$Epithelial[match(SampleID, pi2$sample)]
  # ctOWAS

  w_s <- filtered_weights_s[filtered_weights_s$gene == yName, ]
  i_s <- filtered_intercept_s[match(yName, filtered_intercept_s$gene), ]
  i_s[-1][is.na(i_s[-1])] <- 0

  x_s  <- geno[match(w_s$varID, rownames(geno)), SampleID, drop = FALSE] |> as.matrix() |> t()
  y_s_hat <- pi * (i_s$Cell1 + x_s %*% w_s$Cell1) + (1 - pi) * (i_s$Cell2 + x_s %*% w_s$Cell2)

  # mixcan new
  w_m=filtered_weights_m_new[match(yName, filtered_weights_m_new$gene),]
  i_m=filtered_intercept_m_new[match(yName, filtered_intercept_m_new$gene),];
  i_m[-1][is.na(i_m[-1])] <- 0
  x=geno[match(w_m$varID, rownames(geno)), SampleID, drop = FALSE] %>% as.matrix %>% t()
  y_m_new_hat= pi*(i_m$Cell1 +x%*% w_m$Cell1) + (1-pi)*(i_m$Cell2 +x%*% w_m$Cell2)

  # mixcan old
   w_m=filtered_weights_m_old[match(yName, filtered_weights_m_old$gene),]
  i_m=filtered_intercept_m_old[match(yName, filtered_intercept_m_old$gene),];
  i_m[-1][is.na(i_m[-1])] <- 0
  x=geno[match(w_m$varID, rownames(geno)), SampleID, drop = FALSE] %>% as.matrix %>% t()
  y_m_old_hat= pi*(i_m$Cell1 +x%*% w_m$Cell1) + (1-pi)*(i_m$Cell2 +x%*% w_m$Cell2)



  # prediXcan
  w_p=filtered_weights_p[match(yName, filtered_weights_p$gene),]
  x=geno[match(w_p$varID, rownames(geno)), SampleID, drop = FALSE] %>% as.matrix %>% t()
  y_p_hat= x%*% w_p$tissue

  # get residuals from the pi.
  residuals <- resid(lm(y ~ pi))
  cor1=c(cor(y_s_hat, y), cor(y_m_new_hat, y),  cor(y_m_old_hat, y), cor(y_p_hat, y))
  cor=rbind(cor, cor1)
  print(g)
}

colnames(cor) <- c("ctOWAS","MiXcan (ctOWAS pi)",  "MiXcan (original pi)", "PrediXcan")
cor[is.na(cor)] <- 0
apply(cor, 2, mean, na.rm=T)
apply(cor[which(cor[,1]==0),], 2, mean, na.rm=T)
apply(cor[which(cor[,1]!=0),], 2, mean, na.rm=T)

save(cor,file = file.path(results_dir, "TCGA_Breast_pi2_Cor_Comparison.RData"))


# --------------------------------------
# -------------------- plot -------------------
# -------------------------------------`-
method_cols <- c(
  "PrediXcan"            = "#B8B8B8",
  "MiXcan (original pi)" = "#4F7CAC",
  "MiXcan (ctOWAS pi)"      = "#7FA6D6",
  "ctOWAS"               = "#F0645E"
)
method_levels <- names(method_cols)

## ---------------- Model size ----------------

count_var <- \(x) {
  x |>
    dplyr::distinct(gene, varID) |>
    dplyr::count(gene, name = "n_varID")
}

n_var_wide <-
  count_var(filtered_weights_p) |>
  dplyr::rename(PrediXcan = n_varID) |>
  dplyr::full_join(
    count_var(filtered_weights_m_old) |>
      dplyr::rename(`MiXcan (original pi)` = n_varID),
    by = "gene"
  ) |>
  dplyr::full_join(
    count_var(filtered_weights_m_new) |>
      dplyr::rename(`MiXcan (ctOWAS pi)` = n_varID),
    by = "gene"
  ) |>
  dplyr::full_join(
    count_var(filtered_weights_s) |>
      dplyr::rename(ctOWAS = n_varID),
    by = "gene"
  ) |>
  dplyr::mutate(
    dplyr::across(-gene, \(x) tidyr::replace_na(x, 0L))
  )

## Number of genes with a fitted model
sapply(n_var_wide[-1], \(x) sum(x > 0))

## Mean and median model size among fitted genes
sapply(n_var_wide[-1], \(x) mean(x[x > 0]))
sapply(n_var_wide[-1], \(x) median(x[x > 0]))

n_var_long <- n_var_wide |>
  tidyr::pivot_longer(
    -gene,
    names_to = "Method",
    values_to = "n_varID"
  ) |>
  dplyr::mutate(
    Method = factor(Method, levels = method_levels)
  )

p1 <- ggplot2::ggplot(
  n_var_long,
  ggplot2::aes(Method, log10(n_varID + 1), fill = Method)
) +
  ggplot2::geom_violin(
    trim = FALSE, alpha = 0.4, color = NA
  ) +
  ggplot2::geom_boxplot(
    width = 0.12, outlier.shape = NA,
    alpha = 0.8, linewidth = 0.6
  ) +
  ggplot2::scale_fill_manual(values = method_cols, drop = FALSE) +
  ggplot2::labs(
    y = expression(log[10] * "(# cis-SNPs per gene + 1)"),
    x = NULL
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )
p1
## ---------------- Correlation ----------------

cor_long <- as.data.frame(cor) |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "Method",
    values_to = "Correlation"
  ) |>
  dplyr::mutate(
    Method = factor(Method, levels = method_levels)
  )

## Method-specific modelled genes
cor_nonzero <- cor_long |>
  dplyr::filter(
    is.finite(Correlation),
    Correlation != 0
  )

mean_nonzero_df <- cor_nonzero |>
  dplyr::group_by(Method) |>
  dplyr::summarise(
    mean_cor = mean(Correlation),
    .groups = "drop"
  )

p2 <- ggplot2::ggplot(
  cor_nonzero,
  ggplot2::aes(Method, Correlation, fill = Method)
) +
  ggplot2::geom_boxplot(
    width = 0.5, outlier.alpha = 0.2,
    linewidth = 0.6
  ) +
  ggplot2::scale_fill_manual(values = method_cols, guide = "none") +
  ggplot2::labs(
    title = "Method-Specific Genes",
    x = NULL,
    y = "Correlation"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )

p2
## ---------------- All genes with at least one model ----------------

modelled_genes <- cor_long |>
  dplyr::group_by(gene) |>
  dplyr::filter(
    any(is.finite(Correlation) & Correlation != 0)
  ) |>
  dplyr::ungroup()

mean_all_df <- modelled_genes |>
  dplyr::group_by(Method) |>
  dplyr::summarise(
    mean_cor = mean(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot2::ggplot(
  modelled_genes,
  ggplot2::aes(Method, Correlation, fill = Method)
) +
  ggplot2::geom_boxplot(
    width = 0.5, outlier.alpha = 0.15,
    linewidth = 0.6, na.rm = TRUE
  ) +
  ggplot2::scale_fill_manual(values = method_cols, guide = "none") +
  ggplot2::labs(
    title = "All Modelled Genes",
    x = NULL,
    y = "Correlation"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )

title_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(
    hjust = 0.5, size = 14
  )
)

p3
p1 <- p1 + title_theme
p2 <- p2 + title_theme
p3 <- p3 + title_theme

combined_plot <- cowplot::plot_grid(
  p1, p2, p3,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 14,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1.1, 1, 1)
)

combined_plot

ggplot2::ggsave(
  file.path(figure_dir, "Figure_TCGA.pdf"),
  combined_plot,
  width = 14,
  height = 4.5,
  units = "in",
  device = "pdf"
)
