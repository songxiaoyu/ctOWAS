#' Test K cell-type-specific molecular-trait associations
#'
#' @description
#' Performs ctOWAS association testing for K cell types by combining
#' cell-type-specific SNP prediction weights, GWAS summary statistics and
#' genotype data from a linkage-disequilibrium (LD) reference panel. The
#' function returns separate (marginal) cell-type statistics, and joint
#' statistics adjusted for correlations among the genetically predicted
#' cell-type-specific molecular traits.
#'
#' @param W A numeric P by K matrix of cell-type-specific SNP prediction
#'   weights (for example, \code{W} from \code{\link{ctOWAS_train_K}}). Rows
#'   are SNPs and columns are cell types. SNP order must match the columns of
#'   \code{x_g} and the elements of \code{gwas_z_score}.
#' @param gwas_z_score A numeric vector of length P containing GWAS
#'   z-scores, calculated as the effect-size estimate divided by its standard
#'   error (\eqn{\hat\beta / \mathrm{SE}(\hat\beta)}). Effect alleles must be
#'   aligned with the dosage coding in \code{x_g} and \code{W}.
#' @param x_g A numeric matrix of reference-panel genotype dosages, with
#'   individuals in rows and the same P SNPs in columns.
#' @param n0 Number of controls in the GWAS. Required and must be positive
#'   when \code{family = "binomial"}; not used when
#'   \code{family = "gaussian"}.
#' @param n1 Number of cases in the GWAS. Required and must be positive when
#'   \code{family = "binomial"}; not used when \code{family = "gaussian"}.
#' @param family Outcome family: \code{"binomial"} (default) for a
#'   case-control trait or \code{"gaussian"} for a continuous trait.
#' @param regularization How the diagonal regularization scale is chosen:
#'   \code{"fixed"} (default) uses \code{reg_scale}; \code{"estimate"}
#'   chooses the smallest scale that brings the condition number of the
#'   regularized correlation matrix down to \code{condition_max}, bounded by
#'   \code{reg_min_scale} and \code{reg_max_scale}.
#' @param reg_scale Nonnegative regularization scale used when
#'   \code{regularization = "fixed"}. Default \code{0.1}.
#' @param condition_max Target condition number when
#'   \code{regularization = "estimate"}. Default \code{100}.
#' @param reg_min_scale Lower bound for the estimated scale. This value is
#'   also used when the unregularized matrix already meets
#'   \code{condition_max}, so some regularization is always applied.
#'   Default \code{0.001}.
#' @param reg_max_scale Upper bound for the estimated scale; also used if
#'   the eigenvalues are not finite. Default \code{0.1}.
#'
#' @details
#' \strong{Separate statistics.} The predicted traits in the reference panel
#' are \eqn{\widehat{Y} = X_g W}. For cell type \eqn{k}, the separate
#' statistic is
#'
#' \deqn{Z^{\mathrm{sep}}_k = \frac{\sum_{l} w_{lk}\, \hat\sigma_l\, z_l}
#'                                 {\hat\sigma_{k}},}
#'
#' where \eqn{z_l} is the GWAS z-score of SNP \eqn{l}, \eqn{\hat\sigma_l} is
#' its standard deviation in \code{x_g}, and \eqn{\hat\sigma_k} is the
#' standard deviation of the \eqn{k}th column of \eqn{\widehat{Y}}.
#' Cell types whose predicted trait has zero or non-finite variance (for
#' example, a column of \code{W} that is all zero) get \code{NA}.
#'
#' \strong{Joint statistics.} The columns of \eqn{\widehat{Y}} are
#' standardized to give \eqn{\tilde{Y}}. For a binary outcome, a column of
#' ones is added for the intercept, whose statistic is
#' \eqn{Z_0 = \log(n_1/n_0) / \sqrt{1/n_1 + 1/n_0}}. For a continuous
#' outcome, only the K predicted traits are used. With
#' \eqn{A = \tilde{Y}^\top \tilde{Y}}, \eqn{\Omega = \mathrm{diag}(A)} and
#' \eqn{R} the correlation form of \eqn{A}, the regularized inverse is
#'
#' \deqn{S = \Omega^{-1/2} (R + sI)^{-1} \Omega^{-1/2},}
#'
#' with scale \eqn{s} chosen according to \code{regularization}. The joint
#' statistics are
#'
#' \deqn{Z^{\mathrm{join}} = \mathrm{diag}(S)^{-1/2}\, S\, \Omega^{1/2} Z,}
#'
#' where \eqn{Z} stacks \eqn{Z_0} (binary only) and \eqn{Z^{\mathrm{sep}}}.
#' The intercept element is dropped from the result. Two-sided p-values
#' come from the standard normal distribution.
#'
#' \strong{Fallback.} The joint step is skipped, and the separate
#' statistics are returned in its place, if any cell type has an invalid
#' predicted-trait variance (\code{mode = "invalid_variance_separate"}) or
#' any pair of predicted traits has absolute correlation above 0.999999
#' (\code{mode = "collinear_separate"}). Note that a single cell type with
#' all-zero weights triggers the fallback for the whole test.
#'
#' \strong{Combined p-value.} The K p-values in \code{p_join_vec} are
#' combined with the aggregated Cauchy association test (ACAT) via
#' \code{\link{safe_ACAT}}, which drops missing or out-of-range p-values.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Z_join}{Numeric vector of length K of joint,
#'     correlation-adjusted cell-type z-statistics, or the separate
#'     z-statistics when the joint step is skipped.}
#'   \item{p_join_vec}{Numeric vector of length K of two-sided p-values for
#'     \code{Z_join}.}
#'   \item{p_join}{ACAT-combined p-value across the K cell types, or
#'     \code{NA} if no valid p-values remain or ACAT fails.}
#'   \item{Z_sep}{Numeric vector of length K of separate cell-type
#'     z-statistics (\code{NA} for cell types with invalid variance).}
#'   \item{p_sep}{Numeric vector of length K of two-sided p-values for
#'     \code{Z_sep}.}
#'   \item{reg_scale_selected}{Regularization scale used in the joint step,
#'     or \code{NA} if it was skipped.}
#'   \item{reg_condition}{Approximate condition number (from
#'     \code{kappa(exact = FALSE)}) of the regularized correlation matrix, or
#'     \code{NA} if the joint step was skipped.}
#'   \item{mode}{Character string: \code{"joint"},
#'     \code{"invalid_variance_separate"} or \code{"collinear_separate"}.}
#' }
#'
#' @seealso \code{\link{ctOWAS_train_K}} for training the weights \code{W};
#'   \code{\link{safe_ACAT}}.
#'
#' @examples
#' \dontrun{
#' fit <- ctOWAS_train_K(y = expr_gene, x = geno, pi_k = cell_fractions)
#'
#' res <- ctOWAS_assoc_test_K(
#'   W            = fit$W,          # P x K
#'   gwas_z_score = gwas$beta / gwas$se,
#'   x_g          = ref_geno,       # reference individuals x P
#'   n0           = 20000,
#'   n1           = 5000,
#'   family       = "binomial",
#'   regularization = "estimate"
#' )
#' res$mode
#' res$p_join
#' }
#'
#' @export
ctOWAS_assoc_test_K <- function(W,
                                 gwas_z_score,
                                 x_g,
                                 n0,
                                 n1,
                                 family = c("binomial", "gaussian"),
                                 regularization = c("fixed", "estimate"),
                                 reg_scale = 0.1,
                                 condition_max = 100,
                                 reg_min_scale = 0.001,
                                 reg_max_scale = 0.1) {
  family <- match.arg(family)
  regularization <- match.arg(regularization)

  gwas_z_score <- as.numeric(gwas_z_score)
  if (!is.matrix(W)) {stop("W must be a numeric matrix of dimension p * K.")}
  p <- nrow(W)
  K <- ncol(W)

  if (ncol(x_g) != p) {stop("Number of SNPs (columns) in x_g must match nrow(W).")}
  # calculate separate analysis
  sig_l <- .ctOWAS_col_sds(x_g)
  Yhat <- x_g %*% W
  sig2_g <- .ctOWAS_col_vars(Yhat)
  num <- colSums(W * (gwas_z_score * sig_l))

  Z_sep <- rep(NA_real_, K)
  ok <- sig2_g > 0 & is.finite(sig2_g)
  Z_sep[ok] <- num[ok] / sqrt(sig2_g[ok])
  p_sep <- ifelse(is.na(Z_sep), NA_real_, 2 * stats::pnorm(abs(Z_sep), lower.tail = FALSE))

  Z_join <- Z_sep; p_join_vec <- p_sep
  mode <- "separate"; reg_scale_selected <- NA_real_; reg_condition <- NA_real_

  choose_reg_scale <- function(YtY) {
    if (regularization == "estimate") {
      .estimate_reg_scale_from_cov(YtY, target_condition = condition_max,
        min_scale = reg_min_scale, max_scale = reg_max_scale)
    } else {reg_scale}}

  if (family == "binomial") {
    if (is.null(n0) || is.null(n1) || is.na(n0) || is.na(n1) || n0 <= 0 || n1 <= 0) {
      stop("n0 and n1 must be positive for binomial family.")}

    Z0 <- log(n1 / n0) / sqrt(1 / n1 + 1 / n0)
    if (all(ok)) {
      Y_scaled <- scale(Yhat)
      Y <- cbind(1, Y_scaled)
      colnames(Y) <- c("Y0", paste0("Y", seq_len(K)))
      YtY <- crossprod(Y)
      Omega <- diag(YtY)
      corY <- stats::cov2cor(YtY)
    } else {
      corY <- matrix(NA_real_, nrow = K + 1L, ncol = K + 1L)
    }

    if (any(is.na(corY))) {
      mode <- "invalid_variance_separate"
    } else if (any(abs(corY[upper.tri(corY)]) > 0.999999)) {
      mode <- "collinear_separate"
    } else {
      reg_scale_selected <- choose_reg_scale(YtY)
      inv <- regularized_inverse_cov(YtY, reg_scale = reg_scale_selected)
      v <- c(sqrt(Omega[1]) * Z0, sqrt(Omega[-1]) * Z_sep)
      fit <- .joint_p_from_inverse(inv$inv, v, drop_intercept = TRUE)
      Z_join <- fit$Z_join
      p_join_vec <- fit$p_join_vec
      reg_condition <- inv$condition
      mode <- "joint"
    }
  } else if (family == "gaussian") {
    if (all(ok)) {
      Y_scaled <- scale(Yhat)
      colnames(Y_scaled) <- paste0("Y", seq_len(K))
      YtY <- crossprod(Y_scaled)
      Omega <- diag(YtY)
      corY <- stats::cov2cor(YtY)
    } else {
      corY <- matrix(NA_real_, nrow = K, ncol = K)
    }

    if (any(is.na(corY))) {
      mode <- "invalid_variance_separate"
    } else if (any(abs(corY[upper.tri(corY)]) > 0.999999)) {
      mode <- "collinear_separate"
    } else {
      reg_scale_selected <- choose_reg_scale(YtY)
      inv <- regularized_inverse_cov(YtY, reg_scale = reg_scale_selected)
      v <- sqrt(Omega) * Z_sep
      fit <- .joint_p_from_inverse(inv$inv, v, drop_intercept = FALSE)
      Z_join <- fit$Z_join
      p_join_vec <- fit$p_join_vec
      reg_condition <- inv$condition
      mode <- "joint"
    }
  } else {
    stop("family must be 'binomial' or 'gaussian'")
  }

  p_join <- safe_ACAT(p_join_vec)

  list(
    Z_join = Z_join,
    p_join_vec = p_join_vec,
    p_join = p_join,
    Z_sep = Z_sep,
    p_sep = p_sep,
    reg_scale_selected = reg_scale_selected,
    reg_condition = reg_condition,
    mode = mode
  )
}

# ================= Helper functions ============
regularized_inverse_cov <- function(X, reg_scale = 0.1) {
  X_cor <- stats::cov2cor(X)
  X_reg <- X_cor + reg_scale * diag(nrow(X_cor))
  X_cor_inv <- tryCatch(
    chol2inv(chol(X_reg)),
    error = function(e) solve(X_reg)
  )
  D <- diag(X)
  X_inv <- diag(1 / sqrt(D)) %*% X_cor_inv %*% diag(1 / sqrt(D))

  list(
    inv = X_inv,
    lambda = reg_scale,
    condition = kappa(X_reg, exact = FALSE)
  )
}

.ctOWAS_col_sds <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colSds(x))
  }
  sqrt(colSums((x - matrix(colMeans(x), nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2) /
         (nrow(x) - 1))
}

.ctOWAS_col_vars <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colVars(x))
  }
  colSums((x - matrix(colMeans(x), nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2) /
    (nrow(x) - 1)
}

.estimate_reg_scale_from_cov <- function(X, target_condition = 100,
                                         min_scale = 0.001,
                                         max_scale = 0.1) {
  X_cor <- stats::cov2cor(X)
  eig <- eigen(X_cor, symmetric = TRUE, only.values = TRUE)$values
  lambda_max <- max(eig, na.rm = TRUE)
  lambda_min <- min(eig, na.rm = TRUE)
  if (!is.finite(lambda_max) || !is.finite(lambda_min)) {
    return(max_scale)
  }
  if (lambda_min > 0 && lambda_max / lambda_min <= target_condition) {
    return(min_scale)
  }

  # For X_cor + sI, condition = (lambda_max + s) / (lambda_min + s).
  # Solve for the smallest s that makes condition <= target_condition.
  estimated <- (lambda_max - target_condition * lambda_min) / (target_condition - 1)
  estimated <- max(min_scale, estimated, na.rm = TRUE)
  min(estimated, max_scale)
}

.joint_p_from_inverse <- function(S, v, drop_intercept = FALSE) {
  Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
  Z_join <- as.numeric(Z_full)
  if (drop_intercept) {
    Z_join <- Z_join[-1]
  }
  p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
  list(
    Z_join = Z_join,
    p_join_vec = p_join_vec,
    p_join = safe_ACAT(p_join_vec)
  )
}

safe_ACAT <- function(p_values) {
  p_values <- as.numeric(p_values)
  p_values <- p_values[is.finite(p_values) & p_values >= 0 & p_values <= 1]
  if (length(p_values) == 0L) {
    return(NA_real_)
  }
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)
  })
}
