# Power simulation: generate one dataset, train models, run associations,
# record the result, and then proceed to the next dataset.

library(data.table)
library(foreach)
library(doParallel)
library(doRNG)
library(glmnet)
library(ctOWAS)
library(MiXcan)
library(ggplot2)


# Generate one training/test dataset ------------------------------------------
# n_train = 300; n_test = 3000; b0 = 1; b1 = 1; b2 = 1; eta1 = 0.2; eta2 = 0; var1 = 0.25; var2 = 0.25; group = "heter xy"; seed = 123
generate_one_dataset <- function(X_pool, n_train, n_test, b0, b1, b2, eta1, eta2,
                                 var1 = 0.25, var2 = 0.25,
                                 group = "heter xy", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (!group %in% c("heter xy", "homo xy")) {stop("`group` must be 'heter xy' or 'homo xy'.")}

  n <- n_train + n_test
  p <- ncol(X_pool) - 1L

  sampled_rows <- sample.int(nrow(X_pool), n, replace = FALSE)
  x <- as.matrix(X_pool[sampled_rows, -1L, drop = FALSE])
  storage.mode(x) <- "double"

  snp_names <- colnames(x)
  if (is.null(snp_names)) snp_names <- paste0("SNP", seq_len(p))

  beta1 <- beta2 <- numeric(p + 1L)
  beta1[1L] <- b0

  causal_snp <- sample.int(p, 2L)
  effect_sign <- sample(c(1, 1), 2L, replace = TRUE)

  causal_snp_cell1 <- causal_snp[1L]
  causal_snp_cell2 <- if (group == "heter xy") causal_snp[2L] else causal_snp[1L]

  direction_cell1 <- effect_sign[1L]
  direction_cell2 <- if (group == "heter xy") effect_sign[2L] else effect_sign[1L]

  effect_cell1 <- direction_cell1 * b1
  effect_cell2 <- direction_cell2 * b2

  beta1[causal_snp_cell1 + 1L] <- effect_cell1
  beta2[causal_snp_cell2 + 1L] <- effect_cell2

  design <- cbind(Intercept = 1, x)

  y1 <- as.numeric(design %*% beta1 + rnorm(n, sd = sqrt(var1)))
  y2 <- as.numeric(design %*% beta2 + rnorm(n, sd = sqrt(var2)))

  pi1 <- rbeta(n, 2, 3)
  y <- pi1 * y1 + (1 - pi1) * y2

  disease_lp <- eta1 * y1 + eta2 * y2
  disease_prob <- plogis(disease_lp - mean(disease_lp))
  disease <- rbinom(n, size = 1L, prob = disease_prob)

  train <- seq_len(n_train)
  test <- n_train + seq_len(n_test)

  list(
    y_train = y[train],
    x_train = x[train, , drop = FALSE],
    pi_train = pi1[train],
    x_test = x[test, , drop = FALSE],
    y_test = y[test],
    y1_test=y1[test],
    y2_test=y2[test],
    pi_test = pi1[test],
    outcome_test = disease[test],

    causal_snp_cell1 = causal_snp_cell1,
    causal_snp_name_cell1 = snp_names[causal_snp_cell1],
    effect_cell1 = effect_cell1,

    causal_snp_cell2 = causal_snp_cell2,
    causal_snp_name_cell2 = snp_names[causal_snp_cell2],
    effect_cell2 = effect_cell2
  )
}


# Generate, predict, and test one dataset -------------------------------------
# scenario <- data.table::data.table(scenario_id = "S3", scenario = "eta1 = 0.2, eta2 = 0", eta1 = 0.2, eta2 = 0)

# replicate_id = 1L; scenario_index = 1L; n_train = 300L; n_test = 3000L; b0 = 1; b1 = 1; b2 = 1; group = "heter xy"; nfolds = 10L; prediction_alpha = 0.5; regularization = "fixed"; reg_scale = 0.05; condition_max = 100; reg_min_scale = 0.001; reg_max_scale = 0.1; ld_reference = "test"

run_one_replication <- function(replicate_id, scenario_index, scenario, X_pool, n_train, n_test,
                                b0, b1, b2, group, nfolds, prediction_alpha,
                                regularization, reg_scale, condition_max,
                                reg_min_scale, reg_max_scale, ld_reference = "test") {
  seed <- 100000L * scenario_index + replicate_id

  # generate dat
  sim <- generate_one_dataset(
    X_pool = X_pool, n_train = n_train, n_test = n_test,
    b0 = b0, b1 = b1, b2 = b2, eta1 = scenario$eta1, eta2 = scenario$eta2,
    var1 = 0.25, var2 = 0.25, group = group, seed = seed)

  x_train <- sim$x_train
  x_test <- sim$x_test
  outcome <- sim$outcome_test
  pi_k <- cbind(Cell1 = sim$pi_train, Cell2 = 1 - sim$pi_train)

  p <- ncol(x_train)
  n <-nrow(x_train)
  snp_ids <- colnames(x_train)

  if (is.null(snp_ids)) snp_ids <- paste0("SNP", seq_len(p))
  x_name_matrix <- data.frame(varID = snp_ids, stringsAsFactors = FALSE)

  set.seed(seed)
  foldid <- sample(rep(seq_len(10), length.out = n))
  # 1. Train prediction models on this dataset.
  fit_ctowas <- ctOWAS::ctOWAS_train_K(y = sim$y_train, x = x_train, cov = NULL, pi_k = pi_k, xNameMatrix = x_name_matrix,
    yName = "simulation", foldid=foldid)

  fit_mixcan <- MiXcan::MiXcan(y = sim$y_train, x = x_train,
                               pi = pi_k[, "Cell1"], cov = NULL,
                               xNameMatrix = x_name_matrix, yName = "simulation", foldid = foldid)

  ctowas_weights_all <- as.matrix(fit_ctowas$W)
  mixcan_weights_all <- cbind(Cell1 = as.numeric(unlist(fit_mixcan$beta.SNP.cell1$weight)), Cell2 = as.numeric(unlist(fit_mixcan$beta.SNP.cell2$weight)))
  predixcan_weights_all <- as.numeric(fit_mixcan$beta.all.models[seq_len(p) + 1L, "Tissue"])

  keep_ctowas <- rowSums(abs(ctowas_weights_all), na.rm = TRUE) > 1e-12
  keep_mixcan <- rowSums(abs(mixcan_weights_all), na.rm = TRUE) > 1e-12
  keep_predixcan <- is.finite(predixcan_weights_all) & abs(predixcan_weights_all) > 1e-12

  mean_causal_snp_ctowas <- mean(
    abs(ctowas_weights_all[
      cbind(
        c(sim$causal_snp_cell1, sim$causal_snp_cell2),
        c(1L, 2L)
      )
    ]) > 1e-12
  )

  mean_causal_snp_mixcan <- mean(
    abs(mixcan_weights_all[
      cbind(
        c(sim$causal_snp_cell1, sim$causal_snp_cell2),
        c(1L, 2L)
      )
    ]) > 1e-12
  )

  mean_causal_snp_predixcan <- mean(
    keep_predixcan[
      unique(c(
        sim$causal_snp_cell1,
        sim$causal_snp_cell2
      ))
    ]
  )


  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) < 3L) return(NA_real_)

    sx <- sd(x[ok])
    sy <- sd(y[ok])

    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
      return(NA_real_)
    }

    cor(x[ok], y[ok])
  }


  pred_s <- x_test %*% ctowas_weights_all
  pred_m <- x_test %*% mixcan_weights_all
  pred_tissue_s <- sim$pi_test * (fit_ctowas$intercept[1]+ pred_s[, 1L]) +
            (1 - sim$pi_test) * (fit_ctowas$intercept[2]+pred_s[, 2L])
  pred_tissue_m <- sim$pi_test * (fit_mixcan$beta.all.models[1,2]+ pred_m[, 1L]) +
    (1 - sim$pi_test) * (fit_mixcan$beta.all.models[1,3]+ pred_m[, 2L])
  pred_p <- as.numeric(x_test %*% predixcan_weights_all)

  cor_ct_ctowas <- c( safe_cor(pred_s[, 1L], sim$y1_test),safe_cor(pred_s[, 2L],  sim$y2_test))
  cor_ct_mixcan <- c(safe_cor(pred_m[, 1L], sim$y1_test),safe_cor(pred_m[, 2L],  sim$y2_test))
  cor_tissue_ctowas <- safe_cor(pred_tissue_s, sim$y_test)
  cor_tissue_mixcan <- safe_cor(pred_tissue_m, sim$y_test)
  cor_predixcan <- safe_cor(pred_p, sim$y_test)


  # 2. Run association analyses immediately for this same dataset.
  p_m1 <- p_m2 <- p_m_joint <- NA_real_
  if (any(keep_mixcan)) {
    predicted_cell_expression <- x_test[, keep_mixcan, drop = FALSE] %*% mixcan_weights_all[keep_mixcan, , drop = FALSE]
    fit_m_assoc <- tryCatch(
      MiXcan::MiXcan_association(new_y = predicted_cell_expression, new_cov = NULL,
                                 new_outcome = outcome, family = "binomial"),
      error = function(e) NULL
    )
    if (!is.null(fit_m_assoc) && (is.null(fit_m_assoc$mode) || identical(fit_m_assoc$mode, "joint"))) {
      p_m1 <- fit_m_assoc$cell1_p
      p_m2 <- fit_m_assoc$cell2_p
      p_m_joint <- fit_m_assoc$p_combined
    }
  }

  p_predixcan <- NA_real_
  if (any(keep_predixcan)) {
    predicted_tissue_expression <- as.numeric(x_test[, keep_predixcan, drop = FALSE] %*% predixcan_weights_all[keep_predixcan])
    if (is.finite(stats::sd(predicted_tissue_expression)) && stats::sd(predicted_tissue_expression) > 0) {
      fit_p_assoc <- tryCatch(stats::glm(outcome ~ predicted_tissue_expression, family = stats::binomial()), error = function(e) NULL)
      if (!is.null(fit_p_assoc)) p_predixcan <- summary(fit_p_assoc)$coefficients["predicted_tissue_expression", "Pr(>|z|)"]
    }
  }

  p_s1 <- p_s2 <- p_s_joint <- NA_real_
  if (any(keep_ctowas)) {
    selected_snp <- which(keep_ctowas)
    gwas <- ctOWAS::run_gwas(
      X = x_test[, selected_snp, drop = FALSE], D = outcome,
      family0 = stats::binomial(), method_binomial = "glm"
    )
    x_ld <- if (ld_reference == "test") x_test[, selected_snp, drop = FALSE] else x_train[, selected_snp, drop = FALSE]
    W <- ctowas_weights_all[selected_snp, , drop = FALSE]
    colnames(W) <- c("Cell1", "Cell2")

    fit_s_assoc <- tryCatch(
      ctOWAS::ctOWAS_assoc_test_K(
        W = W, gwas_z_score = gwas$Z, x_g = x_ld,
        n0 = sum(outcome == 0), n1 = sum(outcome == 1), family = "binomial",
        regularization = regularization, reg_scale = reg_scale,
        condition_max = condition_max, reg_min_scale = reg_min_scale, reg_max_scale = reg_max_scale
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_s_assoc)) {
      p_s1 <- fit_s_assoc$p_join_vec[1L]
      p_s2 <- fit_s_assoc$p_join_vec[2L]
      p_s_joint <- fit_s_assoc$p_join
      selected_reg_scale <- fit_s_assoc$reg_scale_selected
      reg_condition <- fit_s_assoc$reg_condition
    }

  }

  data.table(
    scenario_id = scenario$scenario_id, scenario = scenario$scenario,
    eta1 = scenario$eta1, eta2 = scenario$eta2, replicate = replicate_id,
    causal_snp_cell1 = sim$causal_snp_cell1, causal_snp_cell2 = sim$causal_snp_cell2,
    n_snp_ctowas = sum(keep_ctowas), n_snp_mixcan = sum(keep_mixcan), n_snp_predixcan = sum(keep_predixcan),
    mean_causal_snp_ctowas=mean_causal_snp_ctowas,
    mean_causal_snp_mixcan = mean_causal_snp_mixcan,
    mean_causal_snp_predixcan = mean_causal_snp_predixcan,
    cor_tissue_ctowas = cor_tissue_ctowas,
    cor_tissue_mixcan = cor_tissue_mixcan,
    cor_predixcan = cor_predixcan,
    cor_ctowas_cell1 = cor_ct_ctowas[1L],
    cor_ctowas_cell2 = cor_ct_ctowas[2L],
    cor_mixcan_cell1 = cor_ct_mixcan[1L],
    cor_mixcan_cell2 = cor_ct_mixcan[2L],
    p_s_join_1 = p_s1, p_s_join_2 = p_s2, p_s_join = p_s_joint,
    p_m_join_1 = p_m1, p_m_join_2 = p_m2, p_m_join = p_m_joint,
    p_predixcan = p_predixcan
  )
}

