# ==============================================================================
# Example data generation for ctOWAS
# Save as data-raw/example_data.R and run from the package root:
#   source("data-raw/example_data.R")
# This script creates the small strong-signal datasets used in the README and
# package examples for the 2-cell-type workflow.
# ==============================================================================

set.seed(2026)

# ------------------------------------------------------------------------------
# 1. Basic dimensions
# ------------------------------------------------------------------------------
n <- 20
p <- 10
K <- 2

# ------------------------------------------------------------------------------
# 2. Genotype matrix (samples x SNPs)
# Use a relatively high MAF so the example SNPs have clear variation.
# ------------------------------------------------------------------------------
x_example <- matrix(
  rbinom(n * p, 2, 0.5),
  nrow = n,
  ncol = p
)
colnames(x_example) <- paste0("SNP", 1:p)
rownames(x_example) <- paste0("Sample", 1:n)

# ------------------------------------------------------------------------------
# 3. Cell-type fractions (samples x 2 cell types)
# Use near-pure compositions to make the cell-type contrast easy to detect.
# ------------------------------------------------------------------------------
pi_k <- matrix(NA_real_, nrow = n, ncol = K)

# First half of samples: mostly Cell1.
pi_k[1:(n/2), 1] <- runif(n/2, 0.97, 0.995)
pi_k[1:(n/2), 2] <- 1 - pi_k[1:(n/2), 1]

# Second half of samples: mostly Cell2.
pi_k[(n/2 + 1):n, 2] <- runif(n/2, 0.97, 0.995)
pi_k[(n/2 + 1):n, 1] <- 1 - pi_k[(n/2 + 1):n, 2]

colnames(pi_k) <- c("Cell1", "Cell2")
rownames(pi_k) <- rownames(x_example)

# ------------------------------------------------------------------------------
# 4. Cell-type-specific SNP effects
# Assign opposite effects in the two cell types for the first six SNPs so the
# training example is stably classified as CellTypeSpecific.
# ------------------------------------------------------------------------------
b1 <- c( 8, -8,  6, -6,  5, -5, 0, 0, 0, 0)  # Cell1
b2 <- c(-8,  8, -6,  6, -5,  5, 0, 0, 0, 0)  # Cell2

# ------------------------------------------------------------------------------
# 5. Gene expression outcome
# Build a simple mixture of the two cell-type signals plus a small amount of
# noise.
# ------------------------------------------------------------------------------
signal1 <- as.numeric(x_example %*% b1)
signal2 <- as.numeric(x_example %*% b2)

y_example <- pi_k[, 1] * signal1 +
  pi_k[, 2] * signal2 +
  rnorm(n, sd = 0.05)

names(y_example) <- rownames(x_example)

# ------------------------------------------------------------------------------
# 6. Minimal GWAS summary-statistics example
# This object only illustrates the required input format. It is not intended
# for real inference.
# ------------------------------------------------------------------------------
gwas_example <- list(
  Beta    = rnorm(p, 0, 0.05),
  se_Beta = rep(0.05, p)
)
names(gwas_example$Beta)    <- colnames(x_example)
names(gwas_example$se_Beta) <- colnames(x_example)

# ------------------------------------------------------------------------------
# 7. Simulated xCell enrichment scores (xCell cell types x samples)
# Used as xCellScore in estimate_prior(). Scores are a scaled copy of pi_k plus
# a small amount of noise, truncated at zero like real xCell scores:
#   - "Epithelial cells" tracks Cell1;
#   - "Adipocytes" and "Preadipocytes" track Cell2 (pooled into one compartment
#     in the example);
#   - "Fibroblasts" and "Endothelial cells" are unrelated noise, showing that
#     rows not listed in cell_types are ignored.
# Generated before the set.seed(123) below, so all other datasets are unchanged.
# ------------------------------------------------------------------------------
simulate_xcell_score <- function(fraction, scale, noise_sd = 0.04) {
  pmax(scale * fraction + rnorm(length(fraction), 0, noise_sd), 0)
}

xCellScore_example <- rbind(
  "Epithelial cells"  = simulate_xcell_score(pi_k[, "Cell1"], scale = 0.60),
  "Adipocytes"        = simulate_xcell_score(pi_k[, "Cell2"], scale = 0.50),
  "Preadipocytes"     = simulate_xcell_score(pi_k[, "Cell2"], scale = 0.30),
  "Fibroblasts"       = pmax(rnorm(n, 0.10, 0.03), 0),
  "Endothelial cells" = pmax(rnorm(n, 0.05, 0.02), 0)
)
colnames(xCellScore_example) <- rownames(pi_k)   # Sample1, ..., Sample20
xCellScore_example <- round(xCellScore_example, 4)

# ------------------------------------------------------------------------------
# 8. Example genome-wide results for PRIMO post-processing
# Construct a toy results table with both cell-type-specific and non-specific
# signal patterns.
# ------------------------------------------------------------------------------

set.seed(123)

G <- 300
gene <- paste0("Gene", seq_len(G))

# Split the example genes into cell-type-specific and non-specific groups.
type_ct2 <- c(rep("CellTypeSpecific", G/2),
              rep("NonSpecific",     G/2))

# Start from null p-values.
p_1_ct2    <- runif(G)
p_2_ct2    <- runif(G)
p_join_ct2 <- runif(G)

# Add signal to a subset of cell-type-specific genes.
idx_spec <- which(type_ct2 == "CellTypeSpecific")

# Cell1-only signals.
sig1 <- sample(idx_spec, 20)
p_1_ct2[sig1] <- 10^runif(length(sig1), -8, -5)

# Cell2-only signals.
remaining <- setdiff(idx_spec, sig1)
sig2 <- sample(remaining, 20)
p_2_ct2[sig2] <- 10^runif(length(sig2), -8, -5)

# Shared cell-type-specific signals.
remaining <- setdiff(idx_spec, c(sig1, sig2))
shared_spec <- sample(remaining, 15)
p_1_ct2[shared_spec] <- 10^runif(length(shared_spec), -8, -5)
p_2_ct2[shared_spec] <- 10^runif(length(shared_spec), -8, -5)

# Add non-specific joint signals.
idx_uns <- which(type_ct2 == "NonSpecific")
sig_joint <- sample(idx_uns, 25)
p_join_ct2[sig_joint] <- 10^runif(length(sig_joint), -10, -6)

merged_example <- data.frame(
  gene        = gene,
  type_ct2    = type_ct2,
  p_1_ct2     = p_1_ct2,
  p_2_ct2     = p_2_ct2,
  p_join_ct2  = p_join_ct2,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 9. Save example datasets into package data/
# ------------------------------------------------------------------------------
# Run from the package root. Requires the usethis package.
usethis::use_data(
  x_example, pi_k, y_example, gwas_example, xCellScore_example, merged_example,
  overwrite = TRUE
)
