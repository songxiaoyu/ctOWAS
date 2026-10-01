#' Train a Bayesian K-cell-type ctOWAS prediction model
#'
#' This version replaces split-sample screening and elastic-net fitting with a
#' hierarchical mean-plus-cell-type-contrast model fitted by SuSiE.
#'
#' For sample i, the bulk molecular phenotype is modeled as
#'
#'   y_i = sum_k pi_ik * {a_k + x_i' beta_k} + cov_i' gamma + error_i.
#'
#' The K cell-type coefficient vectors are reparameterized into one shared mean
#' effect and K - 1 Helmert cell-type contrasts. SuSiE estimates posterior SNP
#' effects in this transformed space. Cell-type intercepts are unpenalized,
#' whereas both SNP effects and supplied covariates are included in the
#' penalized Bayesian regression. Posterior genetic effects are then transformed
#' back to the original cell types.
#'
#' Cell-type fractions are needed only in the training data. In an application
#' cohort, cell-type-specific genetically regulated expression is predicted by
#'
#'   predicted_expression = new_x %*% W.
#'
#' @param y Numeric vector or one-column matrix of length N.
#' @param x Numeric N by P genotype dosage matrix.
#' @param cov Optional N by Q covariate matrix. Covariates are penalized. Do not
#'   include an intercept.
#' @param pi_k Numeric N by K cell-type composition matrix.
#' @param xNameMatrix Optional P-row SNP annotation.
#' @param yName Optional gene or protein identifier.
#' @param alpha Retained for compatibility; not used by the Bayesian model.
#' @param num_split Retained for compatibility; not used by the Bayesian model.
#' @param seed Optional random seed.
#' @param L Maximum number of SuSiE single-effect components. If NULL, uses at
#'   most 10 components.
#' @param pip_cutoff Default = 0.5. Posterior inclusion probability cutoff used to summarize
#'   selected effects. It is also used for hard thresholding when requested.
#' @param shared_prior_weight Relative prior inclusion weight for the shared
#'   genetic component. If NULL, uses K, giving the shared component more prior
#'   weight than each individual contrast.
#' @param cov_prior_weight Relative prior inclusion weight assigned to every
#'   covariate. The default gives each covariate the same weight as one genetic
#'   contrast effect.
#' @param hard_threshold If TRUE, transformed effects with PIP below pip_cutoff
#'   are set to zero. The default FALSE retains posterior mean effects, which is
#'   generally preferable for prediction.
#' @param coverage Credible-set coverage passed to susieR::susie().
#' @param max_iter Maximum number of SuSiE iterations.
#' @param tol SuSiE convergence tolerance.
#'
#' @return A list containing the original ctOWAS fields plus Bayesian summaries:
#'   PIP, PIP_contrast, beta.contrast, contrast_matrix, susie_fit, pipcut,
#'   Num_SNP_model, and Num_effect_model.
#'
#' @importFrom susieR susie
#' @export

ctOWAS_train_K <- function(
    y,
    x,
    cov = NULL,
    pi_k,
    xNameMatrix = NULL,
    yName = NULL,
    alpha = NULL,
    num_split = NULL,
    seed = NULL,
    L = NULL,
    pip_cutoff = 0.1,
    shared_prior_weight = 1,
    cov_prior_weight = 1,
    hard_threshold = TRUE,
    coverage = 0.99,
    max_iter = 1000L
) {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("Package `susieR` is required. Install it with install.packages('susieR').")
  }

  if (!is.null(seed)) set.seed(seed)

  y <- as.numeric(y)
  x <- as.matrix(x)
  pi_k <- as.matrix(pi_k)
  if (!is.null(cov)) cov <- as.matrix(cov)

  n <- nrow(x)
  p <- ncol(x)
  K <- ncol(pi_k)

  if (length(y) != n) stop("`y` must have length nrow(x).")
  if (nrow(pi_k) != n) stop("`pi_k` must have nrow(x) rows.")
  if (!is.null(cov) && nrow(cov) != n) stop("`cov` must have nrow(x) rows.")
  if (K < 2L) stop("`pi_k` must contain at least two cell types.")
  if (p < 1L) stop("`x` must contain at least one SNP.")
  if (!is.numeric(x) || !is.numeric(pi_k)) stop("`x` and `pi_k` must be numeric.")
  if (!is.numeric(pip_cutoff) || length(pip_cutoff) != 1L || pip_cutoff <= 0 || pip_cutoff >= 1) {
    stop("`pip_cutoff` must be a single number strictly between 0 and 1.")
  }

  # work on x name
  if (is.null(xNameMatrix)) {
    snp_names <- colnames(x)
    if (is.null(snp_names)) snp_names <- paste0("SNP", seq_len(p))
    xNameMatrix <- data.frame(SNP = snp_names, stringsAsFactors = FALSE)
  }
  if (nrow(xNameMatrix) != p) stop("`xNameMatrix` must contain one row per SNP.")
  snp_names <- as.character(xNameMatrix[[1L]])

  # work on missing data
  complete <- is.finite(y) & apply(pi_k, 1L, function(z) all(is.finite(z)))
  if (!is.null(cov)) complete <- complete & apply(cov, 1L, function(z) all(is.finite(z)))
  y <- y[complete]
  x <- x[complete, , drop = FALSE]
  pi_k <- pi_k[complete, , drop = FALSE]
  if (!is.null(cov)) cov <- cov[complete, , drop = FALSE]
  n <- nrow(x)
  if (n < 20L) stop("Fewer than 20 complete training samples remain.")

  # Mean-impute missing genotype dosages. Completely missing SNPs are rejected.
  for (j in seq_len(p)) {
    missing_j <- !is.finite(x[, j])
    if (any(missing_j)) {
      observed_mean <- mean(x[!missing_j, j])
      if (!is.finite(observed_mean)) stop("SNP ", snp_names[j], " contains no finite dosages.")
      x[missing_j, j] <- observed_mean
    }
  }

  if (is.null(colnames(pi_k))) colnames(pi_k) <- paste0("Cell", seq_len(K))
  cell_names <- colnames(pi_k)

  # T maps shared/contrast effects back to the original K cell types:
  # beta_cell = T %*% beta_transformed.
  contrast_matrix <- stats::contr.helmert(K)
  T_matrix <- cbind(Shared = 1, contrast_matrix)
  rownames(T_matrix) <- cell_names
  colnames(T_matrix) <- c("Shared", paste0("Contrast", seq_len(K - 1L)))

  # z_i = pi_i' T. Because compositions are normalized, z[, 1] equals one.
  z <- pi_k %*% T_matrix

  # Only the K cell-type intercept components are unpenalized.
  nuisance_design <- z
  nuisance_qr <- qr(nuisance_design)

  if (nuisance_qr$rank < ncol(nuisance_design)) {
    stop("The cell-type-intercept design is rank deficient.")
  }

  # Genetic design, ordered as P shared SNP effects followed by P effects for
  # each of the K - 1 cell-type contrasts.
  genetic_design <- do.call(
    cbind,
    lapply(seq_len(K), function(k) x * z[, k])
  )

  colnames(genetic_design) <- unlist(
    lapply(colnames(T_matrix), function(component) paste0(component, ":", snp_names)),
    use.names = FALSE
  )

  # Covariates join the SNP terms in the penalized Bayesian design.
  model_design <- if (is.null(cov)) genetic_design else cbind(genetic_design, cov)

  if (!is.null(cov)) {
    cov_names <- colnames(cov)
    if (is.null(cov_names)) cov_names <- paste0("Cov", seq_len(ncol(cov)))
    colnames(model_design)[p * K + seq_len(ncol(cov))] <- paste0("Covariate:", cov_names)
  }

  # Frisch-Waugh-Lovell residualization leaves only cell-type intercepts
  # unpenalized. SNP effects and covariates are fitted jointly by SuSiE.
  y_residual <- qr.resid(nuisance_qr, y)
  model_residual <- qr.resid(nuisance_qr, model_design)

  residual_sd <- apply(model_residual, 2L, stats::sd)
  estimable <- is.finite(residual_sd) & residual_sd > sqrt(.Machine$double.eps)

  if (!any(estimable)) stop("No penalized predictor remains variable after intercept adjustment.")

  if (is.null(shared_prior_weight)) shared_prior_weight <- 1
  if (!is.finite(shared_prior_weight) || shared_prior_weight <= 0) {
    stop("`shared_prior_weight` must be positive.")
  }
  if (!is.finite(cov_prior_weight) || cov_prior_weight <= 0) {
    stop("`cov_prior_weight` must be positive.")
  }

  prior_weights_all <- rep(c(shared_prior_weight, rep(1, K - 1L)), each = p)
  if (!is.null(cov)) prior_weights_all <- c(prior_weights_all, rep(cov_prior_weight, ncol(cov)))
  prior_weights <- prior_weights_all[estimable]
  prior_weights <- prior_weights / sum(prior_weights)

  if (is.null(L)) L <- min(10L, sum(estimable))
  L <- min(as.integer(L), sum(estimable), n - K - 1L)
  if (!is.finite(L) || L < 1L) stop("`L` must be at least one and smaller than the residual degrees of freedom.")

  susie_fit <- susieR::susie(
    X = model_residual[, estimable, drop = FALSE],
    y = y_residual,
    L = L,
    scaled_prior_variance = 0.2,
    residual_variance = NULL,
    estimate_residual_variance = TRUE,
    prior_weights = prior_weights,
    standardize = TRUE,
    intercept = FALSE,
    coverage = coverage,
    max_iter = max_iter,
    tol = 1e-3,
    verbose = FALSE
  )

  # coef() returns effects on the original input scale even though SuSiE
  # standardizes predictors internally. Depending on the susieR version, the
  # returned vector may include the (zero) intercept as its first element.
  posterior_mean_estimable <- as.numeric(stats::coef(susie_fit))
  if (length(posterior_mean_estimable) == sum(estimable) + 1L) {
    posterior_mean_estimable <- posterior_mean_estimable[-1L]
  }
  if (length(posterior_mean_estimable) != sum(estimable)) {
    stop("Unexpected coefficient length returned by `susieR::susie()`.")
  }

  pip_estimable <- as.numeric(susie_fit$pip)

  posterior_mean <- numeric(ncol(model_design))
  pip <- numeric(ncol(model_design))
  posterior_mean[estimable] <- posterior_mean_estimable
  pip[estimable] <- pip_estimable

  genetic_index <- seq_len(p * K)
  covariate_index <- if (is.null(cov)) integer(0) else p * K + seq_len(ncol(cov))
  selected_effect <- pip[genetic_index] >= pip_cutoff
  selected_covariate <- pip[covariate_index] >= pip_cutoff


  if (hard_threshold) {
    posterior_mean[genetic_index[!selected_effect]] <- 0
    posterior_mean[covariate_index[!selected_covariate]] <- 0
  }

  posterior_mean_genetic <- posterior_mean[genetic_index]
  pip_genetic <- pip[genetic_index]
  cov_coef <- posterior_mean[covariate_index]
  cov_pip <- pip[covariate_index]

  beta_contrast <- matrix(
    posterior_mean_genetic,
    nrow = p,
    ncol = K,
    dimnames = list(snp_names, colnames(T_matrix))
  )

  pip_contrast <- matrix(
    pip_genetic,
    nrow = p,
    ncol = K,
    dimnames = list(snp_names, colnames(T_matrix))
  )

  # Transform posterior genetic effects back to the original cell types.
  W <- beta_contrast %*% t(T_matrix)
  colnames(W) <- cell_names
  rownames(W) <- snp_names

  # Estimate unpenalized cell-type intercepts conditional on all penalized
  # posterior mean effects, including covariates.
  nuisance_coef <- qr.coef(nuisance_qr, y - as.numeric(model_design %*% posterior_mean))
  intercept_contrast <- nuisance_coef[seq_len(K)]
  intercept <- as.numeric(T_matrix %*% intercept_contrast)
  names(intercept) <- cell_names

  if (is.null(cov)) {
    beta_all_models <- rbind(intercept = intercept, W)
  } else {
    cov_block <- matrix(
      rep(cov_coef, K),
      nrow = length(cov_coef),
      ncol = K,
      dimnames = list(cov_names, cell_names)
    )
    beta_all_models <- rbind(intercept = intercept, W, cov_block)
  }

  selected_contrast <- matrix(selected_effect, nrow = p, ncol = K)
  selected_snp <- rowSums(selected_contrast) > 0L

  list(
    intercept = intercept,
    W = W,
    beta.all.models = beta_all_models,
    yName = yName,
    xNameMatrix = xNameMatrix,
    pipcut = pip_cutoff,
    PIP = apply(pip_contrast, 1L, max),
    PIP_contrast = pip_contrast,
    beta.contrast = beta_contrast,
    contrast_matrix = T_matrix,
    Num_SNP_input = p,
    Num_SNP_model = sum(selected_snp),
    Num_effect_model = sum(selected_effect),
    Num_covariate_model = sum(selected_covariate),
    covariate_coef = cov_coef,
    PIP_covariate = cov_pip,
    hard_threshold = hard_threshold,
    susie_fit = susie_fit
  )
}
