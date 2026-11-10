# ==============================================================================
# 0. Setup
# ==============================================================================
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(tidyr)
library(ggforce)

paper_dir <- "/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
data_dir <- file.path(paper_dir, "Data")
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")

# for reproducible please change to:
# data_dir=results_dir=figure_dir='~path to uploaded files'
#data_dir=results_dir=figure_dir='/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Manuscript/GB/Reproducibility'


base_font <- 12

# ==============================================================================
# 1. Load data and calculate statistics
# ==============================================================================
dat2 <- read.csv(file.path(results_dir, "DRIVE_all_models_default_all_chr.csv"))
dat3 <- read.csv(file.path(results_dir, "DRIVE_ctOWAS_3ct_default_all_chr.csv"))




# Gene annotation
ensembl38 <- readr::read_csv(file.path(data_dir, "ensembl38.txt"),show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(gene_stable_id = sub("\\..*$", "", gene_stable_id),
         gene_name = trimws(gene_name),
         cytoband = paste0(chromosome_scaffold_name, karyotype_band)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_stable_id) %>%
  summarise(
    gene_name = paste(sort(unique(gene_name)), collapse = "/"),
    cytoband = first(cytoband),
    .groups = "drop")

d2 <- dat2%>%select(Gene_ID, z_o_1, p_o_1, z_o_2, p_o_2, p_o)
#%>% # TO BE UPDATED!!!!!!
#  mutate(FWER_Epithelial_2ct = p.adjust(p_o_1, method = "bonferroni"),
#         FDR_Epithelial_2ct=p.adjust(p_o_1, method = "BH"),
#         FWER_Stromal_2ct = p.adjust(p_o_2, method = "bonferroni"),
#         FDR_Stromal_2ct=p.adjust(p_o_2, method = "BH"),
#         FWER_Tissue_2ct = p.adjust(p_o, method = "bonferroni"),
#         FDR_Tissue_2ct=p.adjust(p_o, method = "BH"))
d3 <- dat3 %>%
  rename(p_o_1=p_1 ,p_o_2 =p_2 ,p_o_3 = p_3 , p_o=p_joint) %>%
  select(Gene_ID, z_1, p_o_1, z_2, p_o_2, z_3, p_o_3, p_o)
#%>%  mutate(FWER = p.adjust(p_joint, method = "bonferroni"),# FDR=p.adjust(p_joint, method = "BH"))


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
c(`2ct` = bacon_lambda(d2$p_o_1), `3ct` = bacon_lambda(d3$p_o_1))
c(`2ct` = bacon_lambda(d2$p_o_2), `3ct` = bacon_lambda(d3$p_o_2))
c(`2ct` = bacon_lambda(d2$p_o), `3ct` = bacon_lambda(d3$p_o))

make_qq <- function(dat, title, pvar = "p_o") {
  qq <- dat %>% transmute(Gene_ID, gene_stable_id = sub("\\..*$", "", Gene_ID), p = .data[[pvar]]) %>%
    filter(is.finite(p), p >= 0, p <= 1) %>%
    mutate(FWER = p.adjust(p, "bonferroni"),) %>%
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
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    labs(title = title, x = expression(Expected~~-log[10](p)), y = expression(Observed~~-log[10](p))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

qq1_2 <- make_qq(d2, "Epithelial", pvar = "p_o_1")
qq2_2 <- make_qq(d2, "Stromal", pvar = "p_o_2")
qq_2 <- make_qq(d2, "Tissue", pvar = "p_o")

qq1_3 <- make_qq(d3, "Epithelial", pvar = "p_o_1")
qq2_3 <- make_qq(d3, "Adipocyte", pvar = "p_o_2")
qq3_3 <- make_qq(d3, "Fibroblast", pvar = "p_o_3")
qq_3 <- make_qq(d3, " Tissue ", pvar = "p_o")

title2 <- cowplot::ggdraw() + cowplot::draw_label("Two Cell Type Model", fontface = "bold", size = 16)
title3 <- cowplot::ggdraw() + cowplot::draw_label("Three Cell Type Model", fontface = "bold", size = 16)
blank <- cowplot::ggdraw()
row2 <- cowplot::plot_grid(qq1_2, qq2_2, blank, qq_2, ncol = 4, labels = c("A", "B", "", "C"), align = "hv")
row3 <- cowplot::plot_grid(qq1_3, qq2_3, qq3_3, qq_3, ncol = 4, labels = c("D", "E", "F", "G"), align = "hv")

p_qq <- cowplot::plot_grid(title2, row2, title3, row3, ncol = 1, rel_heights = c(0.1, 1, 0.1, 1))
ggplot2::ggsave(file.path(figure_dir, "Figure_3ct_Comparison.pdf"), plot = p_qq,
                width = 12, height = 7, units = "in")



# ==============================================================================
# 2. Table S2
# ==============================================================================


dd2 <- d2 %>% # TO BE UPDATED!!!!!!
 mutate(FWER_Epithelial_2ct = p.adjust(p_o_1, method = "bonferroni"),
        FDR_Epithelial_2ct=p.adjust(p_o_1, method = "BH"),
        FWER_Stromal_2ct = p.adjust(p_o_2, method = "bonferroni"),
        FDR_Stromal_2ct=p.adjust(p_o_2, method = "BH"),
        FWER_Tissue_2ct = p.adjust(p_o, method = "bonferroni"),
        FDR_Tissue_2ct=p.adjust(p_o, method = "BH"))

colnames(dd2)[2:6]=c("z_Epithelial_2ct",
  "p_Epithelial_2ct","z_Stromal_2ct", "p_Stromal_2ct", "p_Tissue_2ct")

dd3 <- d3 %>%
  mutate(FWER_Epithelial_3ct = p.adjust(p_o_1, method = "bonferroni"),
         FDR_Epithelial_3ct=p.adjust(p_o_1, method = "BH"),
         FWER_Adipocyte_3ct = p.adjust(p_o_2, method = "bonferroni"),
         FDR_Adipocyte_3ct=p.adjust(p_o_2, method = "BH"),
         FWER_Fibroblast_3ct = p.adjust(p_o_3, method = "bonferroni"),
         FDR_Fibroblast_3ct=p.adjust(p_o_3, method = "BH"),
         FWER_Tissue_3ct = p.adjust(p_o, method = "bonferroni"),
         FDR_Tissue_3ct=p.adjust(p_o, method = "BH"))
colnames(dd3)[2:8]=c("z_Epithelial_3ct", "p_Epithelial_3ct", "z_Adipocyte_3ct", "p_Adipocyte_3ct",
                     "z_Fibroblast_3ct", "p_Fibroblast_3ct", "p_Tissue_3ct")


dat_merge1=dd2 %>% full_join(dd3, by=c("Gene_ID"))%>%
  mutate(gene_stable_id = sub("\\..*$", "", Gene_ID))
dat_merge <- dat_merge1  %>%
  left_join(ensembl38, by = "gene_stable_id") %>%
  dplyr::select(
    gene_name, Gene_ID, cytoband,
    z_Epithelial_2ct, p_Epithelial_2ct, FWER_Epithelial_2ct, FDR_Epithelial_2ct,
    z_Stromal_2ct, p_Stromal_2ct, FWER_Stromal_2ct, FDR_Stromal_2ct,
    p_Tissue_2ct, FWER_Tissue_2ct, FDR_Tissue_2ct,
    z_Epithelial_3ct, p_Epithelial_3ct, FWER_Epithelial_3ct, FDR_Epithelial_3ct,
    z_Adipocyte_3ct, p_Adipocyte_3ct, FWER_Adipocyte_3ct, FDR_Adipocyte_3ct,
    z_Fibroblast_3ct, p_Fibroblast_3ct, FWER_Fibroblast_3ct, FDR_Fibroblast_3ct,
    p_Tissue_3ct, FWER_Tissue_3ct, FDR_Tissue_3ct
  )

# Table S
# TO BE UPDATED!!!!!!
write.csv(dat_merge,
          file.path("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/result/TableS2.csv"),row.names = FALSE)
dat_merge[1:3,]

#filter if any FDR < 0.1
dat_sig <- dat_merge %>%
  filter(if_any(starts_with("FDR_"), ~ !is.na(.x) & .x < 0.1))
write.csv(dat_merge,
          file.path("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/Github/result/TableS2_Sig.csv"),row.names = FALSE)

sum(dat_merge$FWER_Tissue_2ct<0.05, na.rm=T) #9
sum(dat_merge$FWER_Tissue_3ct<0.05, na.rm=T) #9
a=dat_merge[which(dat_merge$FWER_Tissue_2ct<0.05),"gene_name"];a
b=dat_merge[which(dat_merge$FWER_Tissue_3ct<0.05),"gene_name"];b
intersect(a,b)
# "SRGAP2C"    "AC244021.1" "SLC4A7"     "CDYL2"      "SSBP4"      "DNAAF3"

sum(dat_merge$FDR_Tissue_2ct<0.1, na.rm=T) #29
sum(dat_merge$FDR_Tissue_3ct<0.1, na.rm=T) #35
c=dat_merge[which(dat_merge$FDR_Tissue_2ct<0.1),"gene_name"];c
d=dat_merge[which(dat_merge$FDR_Tissue_3ct<0.1),"gene_name"];d
intersect(c,d)
intersect(a,d)
intersect(b,c)

# Key genes:
dat_merge[which(dat_merge$gene_name=="SLC4A7"),]
dat_merge[which(dat_merge$gene_name=="AC244021.1"),]
dat_merge[which(dat_merge$gene_name=="SRGAP2C"),]
