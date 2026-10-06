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
# 1. Load data and calculate statistics
# ==============================================================================
dat2 <- read.csv(file.path(results_dir, "DRIVE_four_models_all_chr.csv"))

dat3 <- read.csv(file.path(results_dir, "DRIVE_ctOWAS_3ct_all_chr.csv"))
sum(dat3$fdr_o<0.1, na.rm=T)

# Calculate FWER using all valid tests in each dataset
prep <- function(dat) {
  dat %>% filter(is.finite(p_o), p_o >= 0, p_o <= 1) %>%
    mutate(gene_stable_id = sub("\\..*$", "", Gene_ID), FWER = p.adjust(p_o, "bonferroni")) %>%
    select(gene_stable_id, p_o, FWER)
}
d2 <- prep(dat2)
d3 <- prep(dat3)
stopifnot(!anyDuplicated(d2$gene_stable_id), !anyDuplicated(d3$gene_stable_id))
# Gene annotation
ensembl38 <- readr::read_csv(file.path(data_dir, "ensembl38.txt"), show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(gene_stable_id = sub("\\..*$", "", gene_stable_id), gene_name = trimws(gene_name)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_stable_id) %>%
  summarise(gene_name = paste(sort(unique(gene_name)), collapse = "/"), .groups = "drop")

# ==============================================================================
# 1. Plot QQ
# ==============================================================================
bacon_lambda <- function(p, seed = 123) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  z <- qnorm(p)  # Equivalent to your earlier qnorm(1 - p, lower.tail = FALSE)
  set.seed(seed)
  as.numeric(bacon::inflation(bacon::bacon(z, na.exclude = TRUE)))
}
c(`2ct` = bacon_lambda(dat2$p_o_1), `3ct` = bacon_lambda(dat3$p_o_1))
c(`2ct` = bacon_lambda(dat2$p_o_2), `3ct` = bacon_lambda(dat3$p_o_2))
c(`2ct` = bacon_lambda(dat2$p_o), `3ct` = bacon_lambda(dat3$p_o))

make_qq <- function(dat, title, pvar = "p_o") {
  qq <- dat %>% transmute(Gene_ID, gene_stable_id = sub("\\..*$", "", Gene_ID), p = .data[[pvar]]) %>%
    filter(is.finite(p), p >= 0, p <= 1) %>%
    mutate(FWER = p.adjust(p, "bonferroni")) %>%
    left_join(ensembl38, by = "gene_stable_id") %>%
    arrange(p) %>%
    mutate(p = pmax(p, 1e-300), expected = -log10(ppoints(n())), observed = -log10(p),
           Label = coalesce(na_if(gene_name, ""), Gene_ID))

  if (!nrow(qq)) stop("No valid p-values in ", pvar)
  lambda <- bacon_lambda(qq$p)

  ggplot(qq, aes(expected, observed)) +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = 2) +
    geom_point(color = "#2C7FB8", size = 1.8, alpha = 0.75) +
    ggrepel::geom_text_repel(data = filter(qq, FWER < 0.05), aes(label = Label), size = 3, fontface = "bold",
                             max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25, min.segment.length = 0, seed = 123) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("lambda[GC] == %.3f", lambda),
             parse = TRUE, hjust = -0.1, vjust = 1.2, size = 4.5, color="red") +
    labs(title = title, x = expression(Expected~~-log[10](p)), y = expression(Observed~~-log[10](p))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

qq1_2 <- make_qq(dat2, "Epithelial", pvar = "p_o_1")
qq2_2 <- make_qq(dat2, "Stromal", pvar = "p_o_2")
qq_2 <- make_qq(dat2, "Tissue", pvar = "p_o")

qq1_3 <- make_qq(dat3, "Epithelial", pvar = "p_o_1")
qq2_3 <- make_qq(dat3, "Adipocyte", pvar = "p_o_2")
qq3_3 <- make_qq(dat3, "Fibroblast", pvar = "p_o_3")
qq_3 <- make_qq(dat3, " Tissue ", pvar = "p_o")

title2 <- cowplot::ggdraw() + cowplot::draw_label("Two Cell Type Model", fontface = "bold", size = 16)
title3 <- cowplot::ggdraw() + cowplot::draw_label("Three Cell Type Model", fontface = "bold", size = 16)

row2 <- cowplot::plot_grid(qq1_2, qq2_2, blank, qq_2, ncol = 4, labels = c("A", "B", "", "C"), align = "hv")
row3 <- cowplot::plot_grid(qq1_3, qq2_3, qq3_3, qq_3, ncol = 4, labels = c("D", "E", "F", "G"), align = "hv")

p_qq <- cowplot::plot_grid(title2, row2, title3, row3, ncol = 1, rel_heights = c(0.1, 1, 0.1, 1))
ggplot2::ggsave(file.path(figure_dir, "Figure_3ct_Comparison.pdf"), plot = p_qq,
                width = 12, height = 7, units = "in")
# ==============================================================================
# 1. Plot contrast
# ==============================================================================

# Compare genes with valid joint p-values in both datasets
cmp <- inner_join(d2, d3, by = "gene_stable_id", suffix = c("_2ct", "_3ct")) %>%
  left_join(ensembl38, by = "gene_stable_id") %>%
  mutate(x = -log10(pmax(p_o_2ct, 1e-300)), y = -log10(pmax(p_o_3ct, 1e-300)),
         Group = case_when(FDR_2ct < 0.1 & FDR_3ct < 0.1 ~ "Both",
                           FDR_2ct < 0.1 ~ "2ct only", FDR_3ct < 0.1 ~ "3ct only", TRUE ~ "Neither"),
         Group = factor(Group, levels = c("Neither", "2ct only", "3ct only", "Both")),
         Label = coalesce(gene_name, gene_stable_id)) %>%
  arrange(Group)

axis_max <- max(c(cmp$x, cmp$y, 1)) * 1.08
r <- cor(cmp$x, cmp$y)

p_compare <- ggplot(cmp, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Group), size = 1.8, alpha = 0.75) +
  geom_text_repel(data = filter(cmp, Group != "Neither"), aes(label = Label, color = Group),
                  size = 3, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                  min.segment.length = 0, seed = 123, show.legend = FALSE) +
  annotate("text", x = 0.04 * axis_max, y = 0.96 * axis_max, label = sprintf("r = %.3f", r),
           hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Neither" = "grey75", "2ct only" = "#377EB8", "3ct only" = "#E41A1C", "Both" = "#984EA3"),
                     drop = FALSE, name = "FDR < 10%") +
  coord_equal(xlim = c(0, axis_max), ylim = c(0, axis_max), clip = "off") +
  labs(x = expression("2 cell types: joint "~-log[10](p)), y = expression("3 cell types: joint "~-log[10](p))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

p_compare
table(cmp$Group)
