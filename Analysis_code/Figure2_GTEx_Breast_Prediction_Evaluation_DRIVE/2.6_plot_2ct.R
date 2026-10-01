# ==============================================================================
# 0. Setup
# ==============================================================================
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(tidyr)
library(Primo)
library(ggforce)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")
paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")
base_font <- 12

# ==============================================================================
# 1. Panel A: Cell Type Composition
# ==============================================================================
load(file.path(results_dir, "GTEx_clean_pi.RData"))

cor_by_marker <- function(marker, p) {
  apply(expr2[intersect(marker, rownames(expr2)), , drop = FALSE], 1, \(x) cor(x, p, use = "pairwise.complete.obs"))
}
make_df <- function(epi, stromal, method) {
  tibble(correlation = c(cor_by_marker(markers$Epithelial, epi), cor_by_marker(markers$Adipocyte, stromal),
                         cor_by_marker(markers$Fibroblast, stromal)),
         cell_type = c(rep("Epithelial", length(intersect(markers$Epithelial, rownames(expr2)))),
                       rep("Stromal", length(intersect(markers$Adipocyte, rownames(expr2))) +
                             length(intersect(markers$Fibroblast, rownames(expr2))))), method = method)
}

plot_df <- bind_rows(make_df(pi2_old$mean_trim_0.05, 1 - pi2_old$mean_trim_0.05, "MiXcan"),
                     make_df(pi2_new[, 1], pi2_new[, 2], "ctOWAS")) %>%
  filter(is.finite(correlation)) %>%
  mutate(cell_type = factor(cell_type, c("Epithelial", "Stromal")), method = factor(method, c("MiXcan", "ctOWAS")))

stats_df <- plot_df %>% group_by(cell_type, method) %>%
  summarise(avg = mean(correlation), med = median(correlation), .groups = "drop") %>%
  mutate(label = sprintf("%s: Mean = %.2f, Median = %.2f", method, avg, med),
         vjust = ifelse(method == "MiXcan", 1.3, 2.8))

stats_df <- plot_df %>% group_by(cell_type, method) %>%
  summarise(Mean = mean(correlation), Median = median(correlation), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = method, values_from = c(Mean, Median)) %>%
  mutate(label = sprintf(
    "Mean: <span style='color:#4C78A8'>%.2f</span> vs <span style='color:#F0645E'>%.2f</span><br>Median: <span style='color:#4C78A8'>%.2f</span> vs <span style='color:#F0645E'>%.2f</span>",
    Mean_MiXcan, Mean_ctOWAS, Median_MiXcan, Median_ctOWAS))

p1 <- ggplot(plot_df, aes(correlation, fill = method, color = method)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  ggtext::geom_richtext(data = stats_df, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
                        hjust = 0, vjust = 1, size = 3.5, fill = NA, label.color = NA, show.legend = FALSE) +
  facet_wrap(~ cell_type, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c(MiXcan = "#4C78A8", ctOWAS = "#F0645E"), breaks = c("MiXcan", "ctOWAS")) +
  scale_color_manual(values = c(MiXcan = "#4C78A8", ctOWAS = "#F0645E"), breaks = c("MiXcan", "ctOWAS")) +
  labs(x = "Marker gene–proportion correlation", y = "Density", fill = NULL, color = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

p1
# ==============================================================================
# 2. Panel B: No of SNPs
# ==============================================================================

load('Results/GTEx_Breast_pi2_ThreeModel_Summary.Rdata')
## ---------------- Model size ----------------
count_var <- \(x) x |> dplyr::distinct(gene, varID) |> dplyr::count(gene, name = "n_varID")

n_var_wide <- count_var(filtered_weights_p) |> dplyr::rename(PrediXcan = n_varID) |>
  dplyr::full_join(count_var(filtered_weights_m_old) |> dplyr::rename(`MiXcan (original pi)` = n_varID), by = "gene") |>
  dplyr::full_join(count_var(filtered_weights_m_new) |> dplyr::rename(`MiXcan (ctOWAS pi)` = n_varID), by = "gene") |>
  dplyr::full_join(count_var(filtered_weights_s) |> dplyr::rename(ctOWAS = n_varID), by = "gene") |>
  dplyr::mutate(dplyr::across(-gene, \(x) tidyr::replace_na(x, 0L)))

## Number of fitted genes; mean and median model size among fitted genes
sapply(n_var_wide[-1], \(x) sum(x > 0))
sapply(n_var_wide[-1], \(x) mean(x[x > 0]))
sapply(n_var_wide[-1], \(x) median(x[x > 0]))
# PrediXcan MiXcan (original pi)   MiXcan (ctOWAS pi)               ctOWAS
# 9.347709             9.439706             9.856737             6.126337

n_var_long <- n_var_wide |> tidyr::pivot_longer(-gene, names_to = "Method", values_to = "n_varID") |>
  dplyr::mutate(Method = factor(Method, levels = method_levels))

p2 <- ggplot2::ggplot(n_var_long, ggplot2::aes(Method, log10(n_varID + 1), fill = Method)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8, linewidth = 0.6) +
  ggplot2::scale_fill_manual(values = method_cols, drop = FALSE) +
  ggplot2::labs(y = expression(log[10] * "(# SNPs + 1)"), x = NULL) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

p2
# ==============================================================================
# 3. Panel C: Correlation with TCGA
# ==============================================================================

load(file.path(results_dir, "TCGA_Breast_pi2_Cor_Comparison.RData"))
method_cols <- c("PrediXcan" = "#B8B8B8", "MiXcan (original pi)" = "#4F7CAC", "MiXcan (ctOWAS pi)" = "#7FA6D6", "ctOWAS" = "#F0645E")
method_levels <- names(method_cols)

cor_long <- tibble::as_tibble(cor, .name_repair = "check_unique") |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  tidyr::pivot_longer(-row_id, names_to = "Method", values_to = "Correlation") |>
  dplyr::mutate(Method = factor(Method, levels = method_levels))


## All genes with at least one model
modelled_genes <- cor_long |> dplyr::group_by(row_id) |>
  dplyr::filter(any(is.finite(Correlation) & Correlation != 0)) |> dplyr::ungroup()

mean_all_df <- modelled_genes |> dplyr::group_by(Method) |>
  dplyr::summarise(mean_cor = mean(Correlation, na.rm = TRUE), .groups = "drop")

p3 <- ggplot2::ggplot(modelled_genes, ggplot2::aes(Method, Correlation, fill = Method)) +
  ggplot2::geom_violin(width = 0.8, alpha = 0.5, trim = TRUE, na.rm = TRUE) +
  ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.15, linewidth = 0.6, na.rm = TRUE) +
  ggplot2::scale_fill_manual(values = method_cols, guide = "none") +
  ggplot2::labs(x = NULL, y = "Pred. Correlation") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

p3


# ==============================================================================
# 4. Panel D: p-value comparison under the same prediction model
# ==============================================================================
drive <- read.csv(file.path(results_dir, "DRIVE_four_models_all_chr.csv"))

plot_df <- dat |> dplyr::filter(is.finite(p_s), is.finite(p_m), p_s > 0, p_s <= 1, p_m > 0, p_m <= 1)
r <- with(plot_df, cor(-log10(p_m), -log10(p_s)))

p4=ggplot2::ggplot(plot_df, ggplot2::aes(x = -log10(p_m), y = -log10(p_s))) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(size = 1.5, alpha = 0.6, color = "#4C78A8") +
  ggplot2::annotate("text", x = -Inf, y = Inf, label = sprintf("r = %.3f", r), hjust = -0.1, vjust = 1.2, size = 4) +
  ggplot2::labs(x = expression(MiXcan~~-log[10](p)), y = expression(ctOWAS[0]~~-log[10](p)))+
  ggplot2::coord_equal() +
  ggplot2::theme_classic(base_size = 14)
p4
# ==============================================================================
# 5. Panel E: Cell Type comparison plots
# ==============================================================================

drive <- read.csv(file.path(results_dir, "DRIVE_four_models_all_chr.csv"))

ensembl38 <- readr::read_csv(file.path(data_dir, "ensembl38.txt"), show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(gene_stable_id = sub("\\..*$", "", gene_stable_id), gene_name = trimws(gene_name)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_stable_id) %>%
  summarise(gene_name = paste(sort(unique(gene_name)), collapse = "/"), .groups = "drop")

plot_data <- drive %>% mutate(
  p_o_1 = pmax(p_o_1, 1e-300), p_o_2 = pmax(p_o_2, 1e-300),
  fdr_o_1 = p.adjust(p_o_1, "BH"), fdr_o_2 = p.adjust(p_o_2, "BH")
  ) %>% mutate(gene_stable_id = sub("\\..*$", "", Gene_ID)) %>%
  left_join(ensembl38, by = "gene_stable_id")


sum(plot_data$fdr_o_1 < 0.1, na.rm = TRUE)
sum(plot_data$fdr_o_2 < 0.1, na.rm = TRUE)

which(plot_data$fdr_o_1 < 0.1)
which(plot_data$fdr_o_2 < 0.1)
intersect(which(plot_data$fdr_o_1 < 0.1), which(plot_data$fdr_o_2 < 0.1))

make_comparison_plot <- function(dat, prefix, title) {
  plot_df <- dat %>%
    transmute(gene_label, p1 = .data[[paste0(prefix, "_1")]], p2 = .data[[paste0(prefix, "_2")]]) %>%
    mutate(across(c(p1, p2), ~ replace(.x, !is.finite(.x) | .x < 0 | .x > 1, NA_real_)),
           fwer1 = p.adjust(p1, "bonferroni"), fwer2 = p.adjust(p2, "bonferroni")) %>%
    filter(!is.na(p1), !is.na(p2)) %>%
    mutate(x = -log10(pmax(p1, 1e-300)), y = -log10(pmax(p2, 1e-300)), significant = fwer1 < 0.05 | fwer2 < 0.05)

  r <- if (nrow(plot_df) > 2 && sd(plot_df$x) > 0 && sd(plot_df$y) > 0) cor(plot_df$x, plot_df$y) else NA_real_
  axis_max <- max(c(plot_df$x, plot_df$y, 1))

  ggplot(plot_df, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7, color = "red") +
    geom_point(alpha = 0.75, size = 2) +
    # geom_text_repel(data = filter(plot_df, significant), aes(label = gene_label), size = 3, fontface = "bold",
    #                max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25, min.segment.length = 0, show.legend = FALSE) +
    annotate("text", x = 0.04 * axis_max, y = 0.9 * axis_max, label = sprintf("r = %.3f", r),
             hjust = 0, vjust = 1, size = 4.2, color = "red") +
    #scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D62728"), limits = c("FALSE", "TRUE"),
    #                   labels = c("FWER ≥ 0.05 in both", "FWER < 0.05 in either"), name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(0, axis_max * 1.05), breaks = pretty(c(0, axis_max))) +
    scale_y_continuous(limits = c(0, axis_max * 1.05), breaks = pretty(c(0, axis_max))) +
    labs(title = title, x = expression(Epithelial~~-log[10](p)), y = expression(Stromal~~-log[10](p))) +
    coord_equal(clip = "off") +
    theme_classic(base_size = 13) +
    theme(axis.line = element_line(linewidth = 0.7), axis.ticks = element_line(linewidth = 0.7),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5), legend.position = "bottom")
}

plot_m <- make_comparison_plot(plot_data, "p_m", "MiXcan")
plot_o <- make_comparison_plot(plot_data, "p_o", "ctOWAS")


# ==============================================================================
# 5. Arrange panels
# ==============================================================================

row1 <- cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.4, 0.3, 0.3), labels = c("A", "B", "C"))
row2 <- cowplot::plot_grid(p4, plot_m, plot_o, nrow = 1,  labels = c("D", "E"))
p_all=cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

ggsave(plot = p_all, filename =file.path(figure_dir, "Figure_Comparison_MiXcan.pdf"), width = 11, height = 6.5)
