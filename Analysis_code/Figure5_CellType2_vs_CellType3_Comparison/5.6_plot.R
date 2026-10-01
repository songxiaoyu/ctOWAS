# ==============================================================================
# UPDATE: Figure 1 becomes 4 panels (All + ct1 + ct2 + ct3), like your first script
#   Figure 1 = 4 scatter panels (All / Cell type 1 / Cell type 2 / Cell type 3)
#   Figure 2 = QQ (A) + Venn (B)   [unchanged from prior reply]
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(tidyr)
library(ggforce)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")

base_font <- 12

# ================================= 1. Load Data =============================================
out3 <- read.csv(file.path(results_dir, "bcac2020_result_pi3_annotated.csv"),colClasses = c(MAP_pattern_nonnull = "character"))
out2 <- read.csv(file.path(results_dir, "bcac2020_result_pi2_comparison_annotated.csv"),colClasses = c(MAP_pattern_nonnull = "character"))


# ===================================  2. QQ plot (A) from out3 ===========================================


pvals_q    <- out3$p_join
genename_q <- out3$gene_name
mask_q <- is.finite(pvals_q) & !is.na(pvals_q) & !is.na(genename_q) & pvals_q > 0 & pvals_q < 1
p_clean <- pvals_q[mask_q]
g_clean <- as.character(genename_q[mask_q])

p <- pmin(pmax(p_clean, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(p, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)


ord <- order(p_clean)
fdr_clean <- p.adjust(p_clean, method = "BH")

df_qq <- data.frame(
  obs  = -log10(p_clean[ord]),
  exp  = -log10(ppoints(length(p_clean))),
  gene = g_clean[ord],
  fdr  = fdr_clean[ord]
)
text_genes_df <- df_qq[which(df_qq$fdr<0.1),]
lambda_label <- data.frame(x = 0,y = max(df_qq$obs, na.rm = TRUE),
                           label = sprintf("bold(lambda)[bold(GC)] == bold(%.3f)", lambda_val))

plot_a <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
  geom_text_repel(data = text_genes_df, aes(label = gene),
                  fontface = "bold", size = 3.5) +
  geom_text(data = lambda_label,aes(x = x, y = y, label = label),inherit.aes = FALSE,
    parse = TRUE, hjust = 0, vjust = 1,size = 5) +
  labs(
    x = expression(bold(Expected ~ -log[10](p))),
    y = expression(bold(Observed ~ -log[10](p)))
  ) +
  theme_classic(base_size = base_font) +
  theme(
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

plot_a




# ================================= 3. QQ plot (A) from out2 =============================================

pvals_q    <- out2$p_join
genename_q <- out2$gene_name
mask_q <- is.finite(pvals_q) & !is.na(pvals_q) & !is.na(genename_q) & pvals_q > 0 & pvals_q < 1
p_clean <- pvals_q[mask_q]
g_clean <- as.character(genename_q[mask_q])

p <- pmin(pmax(p_clean, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(p, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)


ord <- order(p_clean)
fdr_clean <- p.adjust(p_clean, method = "BH")

df_qq <- data.frame(
  obs  = -log10(p_clean[ord]),
  exp  = -log10(ppoints(length(p_clean))),
  gene = g_clean[ord],
  fdr  = fdr_clean[ord]
)
text_genes_df <- df_qq[which(df_qq$fdr<0.1),]
lambda_label <- data.frame(x = 0,y = max(df_qq$obs, na.rm = TRUE),label = sprintf("bold(lambda)[bold(GC)] == bold(%.3f)", lambda_val))

plot_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
  geom_text_repel(data = text_genes_df, aes(label = gene),
                  fontface = "bold", size = 3.5) +
  geom_text(data = lambda_label,aes(x = x, y = y, label = label),inherit.aes = FALSE,
            parse = TRUE, hjust = 0, vjust = 1,size = 5) +
  labs(
    x = expression(bold(Expected ~ -log[10](p))),
    y = expression(bold(Observed ~ -log[10](p)))
  ) +
  theme_classic(base_size = base_font) +
  theme(
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

plot_b
# ================================= 3. Venn plot for cell type comparison =============================================
# plot mannually




