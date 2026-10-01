# Build the A/B/C bar summary for the heterogeneous SNP-Exp setting.
#
# Panels:
#   A. Valid p-value rate under the null, used as a proxy for how often each
#      method returns a usable prediction/association result.
#   B. Type I error from the 2000-replicate null simulation.
#   C. Power for b1 = 0.5, b2 = 1 across four eta settings.
#
# The power panels use the all-replicate denominator so failed/NA models are
# counted as not detected.
library(data.table)
library(ggplot2)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")

paper_dir <- getwd()
input_dir <- file.path(paper_dir, "Results", "Simulation")
figure_dir <- file.path(paper_dir, "Figure")

power_result= fread(file.path(input_dir, "power_fixed_setting_full_results.csv"))
type1_result= fread(file.path(input_dir, "type1_fixed_setting_full_results.csv"))
full_result=rbind(type1_result,power_result )

# Boxplot: number of selected SNPs --------------------------------------------
method_colors <- c(PrediXcan = "#4C78A8", MiXcan = "#F2A541", ctOWAS = "#D62728")
method_levels <- c("PrediXcan", "MiXcan", "ctOWAS")

snp_plot <- melt(full_result, measure.vars = c("n_snp_predixcan", "n_snp_mixcan","n_snp_ctowas"),
  variable.name = "Method", value.name = "Number_of_SNPs")

snp_plot[, Method := factor( Method,levels = c("n_snp_predixcan", "n_snp_mixcan","n_snp_ctowas"),
                             labels = c( "PrediXcan","MiXcan", "ctOWAS"))]


snp_plot[, Method := factor(Method, levels = method_levels)]


p1 <- ggplot(snp_plot,aes(Method, log2(Number_of_SNPs + 1), fill = Method)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual( values = method_colors,breaks = method_levels,drop = FALSE) +
  labs(
    x = NULL,
    y = expression(log[2]("# SNPs" + 1))
  ) +
  theme_bw() +
  theme(legend.position = "none")

p1

snp_plot[, .(Mean_SNPs = mean(Number_of_SNPs, na.rm = TRUE)), by = Method]


# Tissue level correlation and next cell type level correlation  -----------------------------------------------------------
# Cell type 1 = minor; cell type 2 = major.
cor_plot <- rbindlist(list(
  full_result[, .(Method = "PrediXcan", Level = "Tissue", Correlation = cor_predixcan)],
  full_result[, .(Method = "MiXcan", Level = "Tissue", Correlation = cor_tissue_mixcan)],
  full_result[, .(Method = "ctOWAS", Level = "Tissue", Correlation = cor_tissue_ctowas)],
  full_result[, .(Method = "MiXcan", Level = "Minor Cell Type", Correlation = cor_mixcan_cell1)],
  full_result[, .(Method = "MiXcan", Level = "Major Cell Type", Correlation = cor_mixcan_cell2)],
  full_result[, .(Method = "ctOWAS", Level = "Minor Cell Type", Correlation = cor_ctowas_cell1)],
  full_result[, .(Method = "ctOWAS", Level = "Major Cell Type", Correlation = cor_ctowas_cell2)]
))

method_levels <- c("PrediXcan", "MiXcan", "ctOWAS")
cor_plot[, `:=`(Method = factor(Method, levels = method_levels),
                Level = factor(Level, levels = c("Tissue", "Minor Cell Type", "Major Cell Type")))]

make_cor_plot <- function(level) {
  ggplot(cor_plot[Level == level & is.finite(Correlation)], aes(Method, Correlation, fill = Method)) +
    geom_violin(trim = TRUE, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.alpha = 0.2, fill = "white") +
    scale_fill_manual(values = method_colors) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(title = level, x = NULL, y = "Prediction correlation") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
}

p_tissue <- make_cor_plot("Tissue")
p_major  <- make_cor_plot("Major Cell Type")
p_minor  <- make_cor_plot("Minor Cell Type")

p2 <- cowplot::plot_grid(p_tissue, p_major, p_minor, nrow = 1,align = "hv",
                            rel_widths = c(1.3, 1, 1))
p2

cor_plot[, .(Mean = mean(Correlation, na.rm = TRUE), Median = median(Correlation, na.rm = TRUE)), by = .(Method, Level)]



# Plot values -----------------------------------------------------------

setDT(full_result)
p_columns <- c(PrediXcan = "p_predixcan", MiXcan    = "p_m_join", ctOWAS    = "p_s_join",
               MiXcan_1    = "p_m_join_1", ctOWAS_1    = "p_s_join_1",
               MiXcan_2    = "p_m_join_2", ctOWAS_2    = "p_s_join_2")

method_colors <- c(PrediXcan = "#4C78A8", MiXcan    = "#F2A541",ctOWAS    = "#D62728",
  MiXcan_1  = "#F6C36A",ctOWAS_1  = "#E66B6B",MiXcan_2  = "#C98720",ctOWAS_2  = "#A71919")


scenario_levels <- unique(full_result[order(scenario_id), as.character(scenario)])
scenario_labels <- setNames(c(
  'atop("No mQT-Disease", paste(eta[1], " = ", eta[2], " = 0"))',
  'atop("Homogeneous", paste(eta[1], " = ", eta[2], " = 0.2"))',
  'atop("Minor Cell Type", paste(eta[1], " = 0.2; ", eta[2], " = 0"))',
  'atop("Major Cell Type", paste(eta[1], " = 0; ", eta[2], " = 0.2"))',
  'atop("Opposite Directions", paste(eta[1], " = -0.2; ", eta[2], " = 0.2"))'
), scenario_levels)
scenario_levels <- scenario_levels[c(1, 3, 4, 2, 5)]
scenario_labels <- scenario_labels[scenario_levels]

plot_data <- melt(
  full_result,
  id.vars = c("scenario_id", "scenario"),
  measure.vars = unname(p_columns),
  variable.name = "p_column",
  value.name = "p_value"
)[
  , Method := names(p_columns)[match(p_column, p_columns)]
][
  , Level := fcase(
    Method %chin% c("PrediXcan", "MiXcan", "ctOWAS"), "Tissue Level",
    grepl("_1$", Method), "Minor Cell Type",
    grepl("_2$", Method), "Major Cell Type"
  )
][
  is.finite(p_value),
  .(
    Rate = mean(p_value < 0.05)
  ),
  by = .(scenario_id, scenario, Level, Method)
]

plot_data[, `:=`(
  scenario = factor(scenario, levels = scenario_levels),
  Level = factor(
    Level,
    levels = c("Minor Cell Type", "Major Cell Type", "Tissue Level")
  ),
  Method = factor(Method, levels = names(p_columns)),
  Method_label = sub("_[12]$", "", Method)
)]


plot_data[, `:=`(
  Level = factor(Level, levels = c("Minor Cell Type", "Major Cell Type", "Tissue Level")),
  Method_label = factor(sub("_[12]$", "", as.character(Method)),
                        levels = c("PrediXcan", "MiXcan", "ctOWAS"))
)]

p_power <- ggplot(
  plot_data,
  aes(x = Method_label, y = Rate, fill = Method)
) +
  geom_col(
    width = 0.68,
    color = "grey25",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Rate)),
    vjust = -0.35,
    size = 3
  ) +
  facet_grid(
    rows = vars(Level), cols = vars(scenario),
    labeller = labeller(
      Level = label_value,
      scenario = as_labeller(scenario_labels, label_parsed)
    )
  ) +
  scale_fill_manual(
    values = method_colors,
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.10))
  ) +
  coord_cartesian(
    ylim = c(0, 1.05),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Rejection rate"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90",color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 35,hjust = 1),
    panel.spacing.x = unit(0.7, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

type1_panels <- unique(plot_data[
  as.character(scenario) == scenario_levels[1] |
    (Level == "Minor Cell Type" & as.character(scenario) == scenario_levels[3]) |
    (Level == "Major Cell Type" & as.character(scenario) == scenario_levels[2]),
  .(Level, scenario)
])

p_power <- p_power +
  geom_text(data = type1_panels, aes(x = 2, y = 1.02), label = "Type I Error",
            inherit.aes = FALSE, color = "purple", size = 3.5)

p_power



# Arrange figures ----------------------


library(cowplot)

row1 <- plot_grid(p2, p1, nrow = 1, rel_widths = c(2.5, 1), align = "h", labels = c("A", "B"), label_size = 16)
p_all <- plot_grid(row1, p_power, ncol = 1, rel_heights = c(1, 2.5), labels = c("", "C"), label_size = 16)
ggsave("combined_plot.pdf", p_all, width = 15, height = 10)


ggsave(plot = p_all, filename =  file.path(figure_dir, "Figure_Simulation.pdf"),width = 9, height = 8)
