library(ctOWAS)

# Load example marker signatures and gene expression data
data(markers)
data(expr)
markers

# Extract cell-type names (remove suffixes after "#")
cell_names <- unique(sub("#.*$", "", names(markers)))

# Construct a custom xCell2 reference using marker signatures
custom_xcell2 <- methods::new("xCell2Object", signatures = markers,
                              genes_used = unique(unlist(markers, use.names = FALSE)),
  # The rests are placeholder to run algorithm; do not impact results
  dependencies = list(),
  params = data.frame(celltype = cell_names,a = 1, b = 1, m = 1, n = 0),
  spill_mat = diag( 1, length(cell_names),dimnames = list(cell_names, cell_names))
)

# Calculate cell-type enrichment scores
xCellScore <- xCell2::xCell2Analysis(mix = expr,xcell2object = custom_xcell2,
   minSharedGenes = 0.8,rawScores = TRUE, spillover = FALSE)


# Estimate pi prior ----------


cell_type_3 <- as.list(cell_names); names(cell_type_3)=cell_names
pi_prior <- estimate_prior(xCellScore = xCellScore, cell_types = cell_type_3,
                           mu = c(0.4, 0.5, 0.1),
                           lower = 0.01, upper = 0.98)
pi_prior[1:3,]
# Epithelial Adipocyte Fibroblast
# GTEX-1122O  0.2138548 0.5756515 0.21049366
# GTEX-11EM3  0.6603801 0.2180780 0.12154189
# GTEX-11EMC  0.5025893 0.4516501 0.04576063
# cite: Muse, Meghan E., et al. "Application of novel breast biospecimen cell-type adjustment identifies shared DNA methylation alterations in breast tissue and milk with breast cancer–risk factors." Cancer Epidemiology, Biomarkers & Prevention 32.4 (2023): 550-560.

# Estimate pi   ----------
marker_clean <- setNames(lapply(cell_names, function(ct) unique(unlist(markers[grepl(paste0("^", ct, "#"), names(markers))]))), cell_names)

pi_hat <- pi_estimation_K(
  expr        = expr,    # genes x samples
  n_iteration = 20,               # number of subsampling iterations
  prior       = pi_prior,            # samples x cell types, from step (a)
  marker_list = marker_clean,  # named list, few elements per cell type
  seed        = 1,
  verbose     = FALSE             # passed to deNet_composition_K_once()
)

pi_hat[1:3,]      # samples x cell types; rows sum to 1
# Epithelial Adipocyte Fibroblast
# GTEX-1122O  0.2795854 0.4942176 0.22619695
# GTEX-11EM3  0.6319745 0.2465108 0.12151471
# GTEX-11EMC  0.5620261 0.3855683 0.05240558
