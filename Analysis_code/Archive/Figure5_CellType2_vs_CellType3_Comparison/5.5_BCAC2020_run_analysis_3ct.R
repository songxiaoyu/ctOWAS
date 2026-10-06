# --- Stage 2: run_ctOWAS_analysis.R ---

library(data.table)
library(lme4)
library(glmnet)
library(doParallel)
library(doRNG)
library(ACAT)
library(ctOWAS)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
bcac_input_dir <- file.path(paper_dir, "Results", "bcac2020_input")
ld_input_dir <- file.path(paper_dir, "Results", "bcac2020_filtered_id")
result_dir <- file.path(paper_dir, "Results", "bcac2020_result")

# STEP1: load data--------
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

n1 <- 133384
n0 <- 113789 + 18908


# Set chromosome
for(chr in 1:22){
  print('CHR')
  print(chr)
  # Load pre-merged data
  # Keep the historical workspace folder/file naming used by the archived 3pi run.
  mw_gwas_input_path <- file.path(bcac_input_dir, sprintf("chr%d_mw_gwas_input_bcac2020_oct.rds", chr))
  mw_gwas_input <- readRDS(mw_gwas_input_path)

  # Load LD matrix, ref genome, SNP info
  snp_id <- fread(file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_pi3.snplist", chr)), header = FALSE)
  ref_snp_id <- as.character(snp_id$V1)

  # Read genotype matrix and set proper colnames
  X_ref_dt <- fread(file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_pi3.raw", chr)))
  X_ref <- as.matrix(X_ref_dt[, 7:ncol(X_ref_dt)])
  colnames(X_ref) <- sub("_.*", "", colnames(X_ref))

  # Filter intersecting SNPs
  nrow(mw_gwas_input)
  gwas_ref_snps <- intersect(mw_gwas_input$varID, ref_snp_id)
  length(gwas_ref_snps)
  filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$varID %in% gwas_ref_snps, ]

  # Prepare filtered gene list
  split_df <- split(filtered_mw_gwas_input, filtered_mw_gwas_input$gene)
  filtered_list <- list()
  for (gene in names(split_df)) {
    gene_df <- split_df[[gene]]
    W1 <- gene_df$Cell1
    W2 <- gene_df$Cell2
    W3 <- gene_df$Cell3
    W <- cbind(W1, W2, W3)
    filtered_list[[gene]] <- list(W = W, selected_snp = gene_df)}


  # STEP2 Run ctOWAS----
  G <- length(filtered_list)
  real_result = data.frame(matrix(ncol = 11, nrow = G))
  colnames(real_result) <- c('gene','varID','chr','input_snp_num','Z_1','p_1','Z_2','p_2','Z_3','p_3','p_join')
  real_result$chr <- chr

  for (g in 1:G) {
    gene = names(split_df)[g]
    cat("Processing gene:", gene, "\n")

    W <- filtered_list[[gene]]$W
    selected_snp_id <- filtered_list[[gene]]$selected_snp$varID
    gwas_z_score <- filtered_list[[gene]]$selected_snp$beta.Gwas/filtered_list[[gene]]$selected_snp$SE.Gwas

    X_ref_filtered <- X_ref[, selected_snp_id, drop = FALSE]
    ctOWAS_results <- ctOWAS_assoc_test_K(W=W, gwas_z_score=gwas_z_score, x_g=X_ref_filtered,
                                             n0=n0, n1=n1, family='binomial')
    real_result[g, c('gene','varID','chr', 'input_snp_num')] =
      c(filtered_list[[gene]]$selected_snp[1,c('gene','varID','CHR')], nrow(W))
    real_result[g, c('Z_1','p_1','Z_2','p_2','Z_3','p_3','p_join')] <-
      c(c(rbind(ctOWAS_results$Z_join, ctOWAS_results$p_join_vec)), ctOWAS_results$p_join)
  }

  result_path <- file.path(result_dir, sprintf("bcac2020_chr%d_result_pi3.csv", chr))
  write.csv(real_result, result_path, row.names = FALSE)
}




# # clean results -------------
combined3 <- do.call(rbind, lapply(1:22, function(chr) {
  result_path <- file.path(result_dir, sprintf("bcac2020_chr%d_result_pi3.csv", chr))
  read.csv(result_path)
}))
# combined_path <- file.path(result_dir, "bcac2020_result_pi3.csv")
# write.csv(combined3, combined_path, row.names = FALSE)
# #---Input----
# combined3 = read.csv(file.path(results_dir, "bcac2020_result_pi3.csv"))

#---Add anotations----
ensembl_ref = read.csv(file.path(data_dir, "ensembl38.txt"))
setDT(combined3)
combined3[, gene_id := sub("\\..*$", "", gene)]  # drop any ENSG version
setDT(ensembl_ref)

# Build ENSG -> cytoband map (deduplicate just in case)
ensembl_cyto <- unique(ensembl_ref[, .(ENSG = sub("\\..*$", "", Gene.stable.ID),
                                       CYTOBAND = Karyotype.band, gene_name = Gene.name)],
                       by = "ENSG")

# Join cytoband by ENSG
combined3j <- ensembl_cyto[combined3, on = .(ENSG = gene_id)]


# Calculate fwer cutoff
combined3j$fwer_p_join <- p.adjust(combined3j$p_join, method = "bonferroni")
combined3j$fdr_p_join <- p.adjust(combined3j$p_join, method = "BH")
length(which(combined3j$fwer_p_join < 0.05))  #
0.05 / nrow(combined3j) # fwer cutoff

length(which(combined3j$fdr_p_join < 0.1))
fdr_cutoff <- max(combined3j$p_join[combined3j$fdr_p_join  < 0.1], na.rm=TRUE)
fdr_cutoff # fdr cutoff




# Write table S1
out3_rename <- combined3j %>%
  dplyr::rename(
    Z_epithelial = Z_1,
    p_epithelial = p_1,
    Z_adipocyte = Z_2,
    p_adipocyte = p_2,
    Z_fibroblast= Z_3,
    p_fibroblast = p_3,
    p_joint     = p_join,
    fwer_joint = fwer_p_join,
    fdr_joint = fdr_p_join
  ) %>%
  dplyr::select(
    gene_name, ENSG, chr, CYTOBAND, input_snp_num, Z_epithelial, p_epithelial,
    Z_adipocyte, p_adipocyte, Z_fibroblast, p_fibroblast,
    p_joint, fwer_joint, fdr_joint
  )

write.csv(out3_rename, file.path(results_dir, "bcac2020_result_pi3_annotated.csv"), row.names = FALSE)
write.csv(out3_rename, file.path(results_dir, "tableS3.csv"), row.names = FALSE)
