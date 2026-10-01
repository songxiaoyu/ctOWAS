rm(list=ls())
library(dplyr)
library(ctOWAS)
library(MiXcan)
library(purrr)
library(data.table)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()
results_dir <- file.path(paper_dir, "Results")
drive_pheno_dir <- file.path(paper_dir, "Data", "DRIVE")
drive_geno_dir <- file.path('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_MiXcan/Data/OncoArray')




# Weights
load(file.path(results_dir, "GTEx_Breast_pi2_ThreeModel_Summary.RData"))
merge_keys <- c("gene", "rsid", "varID", "ref_allele", "eff_allele")

mw_input <- list(
  filtered_weights_s %>% dplyr::rename(ctOWAS_Cell1 = Cell1, ctOWAS_Cell2 = Cell2),
  filtered_weights_m_new %>% dplyr::rename(MiXcan_Cell1 = Cell1, MiXcan_Cell2 = Cell2),
  filtered_weights_p %>% dplyr::rename(PrediXcan_tissue = tissue)
) %>%
  purrr::reduce(full_join, by = merge_keys) %>%
  mutate(
    across(c(ctOWAS_Cell1, ctOWAS_Cell2, MiXcan_Cell1, MiXcan_Cell2, PrediXcan_tissue), ~ coalesce(.x, 0)),
    CHR = as.integer(sub("^chr([0-9]+)_.*$", "\\1", varID)),
    POS = as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", varID))
  ) %>%
  select(
    gene, rsid, varID, CHR, POS, ref_allele, eff_allele,
    ctOWAS_Cell1, ctOWAS_Cell2,
    MiXcan_Cell1, MiXcan_Cell2,
    PrediXcan_tissue
  ) %>%
  arrange(CHR, POS, gene)


# DRIVE
drive_pheno <- read.csv(file.path(drive_pheno_dir, "oncoarray-drive.pheno.csv"))
pca_result <- fread(file.path(drive_geno_dir,"oncoarray_principal_components_fastmode_2020-10-17.txt")) %>%
  dplyr::rename(FID = subject_ID)
names(pca_result)[2:11] <- paste0("PC", 1:10)
covariates <- drive_pheno %>%
  select(subject_index, FID, genotyping_center, age_int) %>%
  inner_join(pca_result) %>%
  select(-FID) %>%
  mutate(subject_index = as.character(subject_index)) %>%
  mutate(age_int = case_when(age_int == 888 ~ NA_real_,
                             TRUE ~ age_int)) %>%
  filter(PC1 <= 0.006)

setDT(mw_input)



for (chr in seq_len(22L)) {

  message("\nProcessing chromosome ", chr, " of 22")

  # STEP 1 READ DATA-------------

  # load genomic data that's in prediction models
  genotype_file <- file.path(drive_geno_dir, sprintf("oncoarray_dosages_chr%02d.txt.gz", chr))
  key_file <- file.path(results_dir,sprintf("variant_keys_chr%02d.txt", chr))
  subset_file <- tempfile(fileext = ".txt")

  variant_keys <- unique(mw_input[CHR == chr, .(CHR_key = paste0("chr", CHR),POS,
                                                REF = ref_allele,ALT = eff_allele)])

  fwrite(variant_keys, key_file, sep = "\t", col.names = FALSE)

  cmd <- sprintf(
    "gzip -dc %s | awk 'NR==FNR {k[$1 FS $2 FS $3 FS $4]=1; next} FNR==1 || (($1 FS $2 FS $4 FS $5) in k)' %s - > %s",
    shQuote(genotype_file),
    shQuote(key_file),
    shQuote(subset_file)
  )

  status <- system(cmd)

  header_names <- scan(text = readLines(subset_file, n = 1L), what = "", quiet = TRUE)

  drive_genome <- fread(file = subset_file, header = FALSE, skip = 1L, data.table = TRUE, check.names = FALSE)

  setnames(drive_genome, c("CHR", header_names[-1L]))

  drive_genome[, `:=`(CHR = as.integer(sub("^chr", "", CHR)),
    POS = as.integer(POS),
    REF = toupper(REF),
    ALT = toupper(ALT))]

  unlink(c(key_file, subset_file))

  # merge two datasets
  mw_drive_input <- merge(mw_input, drive_genome,
                        by.x = c("CHR", "POS", "ref_allele", "eff_allele"),
                        by.y = c("CHR", "POS", "REF", "ALT"))
  # dim(mw_drive_input)

  X_genome=t(mw_drive_input[,18:ncol(mw_drive_input)])
  dim(X_genome)
  # Sample ID
  SampleID=rownames(X_genome)

  #STEP 2 RUN GWAS------------------

  # Match rows of drive_pheno to X_genome by subject_index
  D <- drive_pheno[match(SampleID, drive_pheno$subject_index), "affection_status" ]
  cov<- covariates[match(SampleID, drive_pheno$subject_index),-1]
  # run gwas
  drive_gwas_result = run_gwas(X=X_genome, D=D, family0 =stats::binomial(), covar=cov,
                             method_binomial = "score")
  drive_gwas_result = data.frame(drive_gwas_result)

  # combine gwas results with inpit
  drive_input_total = cbind(drive_gwas_result, mw_drive_input)


  #STEP 3 RUN four models (MiXcan, PrediXcan, ctOWAS, ctOWAS with MiXcan weights) --------------
  # Note: Here are helper function below. Run them first.

  #split DATA BY GENE
  drive_gene <-unique(mw_drive_input$gene)
  drive_split_df <- split(drive_input_total, drive_input_total$gene)

  drive_result = data.frame(matrix(ncol = 10, nrow = length(drive_gene)))
  colnames(drive_result) <- c('p_o_1','p_o_2','p_o',
                              'p_s_1','p_s_2','p_s','p_m_1','p_m_2','p_m', 'p_predixcan')
  print(c("Total number of genes", length(drive_gene)))
  # by gene
  for (g in 1:length(drive_gene)) {
    # for (g in 1:3) {
    gene = drive_gene[g]
    gene_df <- drive_split_df[[gene]]  # Access each dataframe
    print(paste("Processing gene:", gene))
    new_x=gene_df[,22:ncol(gene_df)] %>% t()
    dim(new_x)

    # run ctOWAS
    W <- as.matrix(gene_df[, c("ctOWAS_Cell1", "ctOWAS_Cell2")])

    if (any(W!=0)) {
     gwas_results = gene_df$Z
     ctOWAS_res <- ctOWAS_assoc_test_K(W = W, gwas_z_score = gwas_results, x_g = new_x,
                                       n0 = sum(D == 0, na.rm = TRUE), n1 = sum(D, na.rm = TRUE),
                                       family = "binomial")
     drive_result[g, c('p_o_1','p_o_2','p_o')] =c(ctOWAS_res$p_join_vec[1:2], ctOWAS_res$p_join)
    }

    # run ctOWAS with MiXcan prediction model
    W <- as.matrix(gene_df[, c("MiXcan_Cell1", "MiXcan_Cell2")])
    if (any(W!=0)) {
      gwas_results = gene_df$Z
      ctOWAS_res <- ctOWAS_assoc_test_K(W = W, gwas_z_score = gwas_results,x_g = new_x,
                                        n0 = sum(D == 0, na.rm = TRUE), n1 = sum(D, na.rm = TRUE), family = "binomial")
      drive_result[g, c("p_s_1", "p_s_2", "p_s")] =c(ctOWAS_res$p_join_vec[1:2], ctOWAS_res$p_join)
    }

   # Run MiXcan
    W <- as.matrix(gene_df[, c("MiXcan_Cell1", "MiXcan_Cell2")])
    if (any(W!=0)) {
     pred=MiXcan::MiXcan_prediction(weight=W, new_x=new_x)
     mixcan_res=MiXcan_association(new_y=pred, new_outcome=D, new_cov=cov, family = "binomial")
     drive_result[g, c("p_m_1", "p_m_2", "p_m")] <- c(mixcan_res$cell1_p, mixcan_res$cell2_p, mixcan_res$p_combined)
    }

    # Run PrediXcan
    W <- as.matrix(gene_df[, c("PrediXcan_tissue")])
    if (any(W!=0)) {
    pred=new_x %*%W
    PrediXcan_res=glm(D~pred+cov, family = "binomial")
    drive_result[g, 'p_predixcan'] = summary(PrediXcan_res)$coefficients[2,4]
    }

  }
  rownames(drive_result) <- drive_gene
  # STEP 4: Save chromosome-specific results
  drive_result <- drive_result %>% tibble::rownames_to_column(var = "Gene_ID")
  fwrite(drive_result,file = file.path(results_dir, sprintf("DRIVE_four_models_chr%02d.csv", chr)))

  rm(drive_genome, mw_drive_input, X_genome, drive_gwas_result, drive_input_total, drive_split_df)
  gc()
}


# Merge results from 22 chr


files <- file.path(results_dir, sprintf("DRIVE_four_models_chr%02d.csv", 1:22))
drive_result_all <- data.table::rbindlist(lapply(files[file.exists(files)], data.table::fread), use.names = TRUE, fill = TRUE)
data.table::fwrite(drive_result_all, file.path(results_dir, "DRIVE_four_models_all_chr.csv"))

View(drive_result_all)

plot(-log10(drive_result_all$p_m), -log10(drive_result_all$p_s))
abline(0,1)
cor(-log10(drive_result_all$p_m), -log10(drive_result_all$p_s), use="pairwise")


plot(-log10(drive_result_all$p_m), -log10(drive_result_all$p_o))
abline(0,1)
cor(-log10(drive_result_all$p_m), -log10(drive_result_all$p_o), use="pairwise")



plot(-log10(drive_result_all$p_o_1), -log10(drive_result_all$p_o_2))
abline(0,1)
cor(-log10(drive_result_all$p_o_1), -log10(drive_result_all$p_o_2), use="pairwise")

plot(-log10(drive_result_all$p_m_1), -log10(drive_result_all$p_m_2))
abline(0,1)
cor(-log10(drive_result_all$p_m_1), -log10(drive_result_all$p_m_2), use="pairwise")
