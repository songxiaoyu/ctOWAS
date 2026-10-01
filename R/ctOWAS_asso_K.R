#' Compute a regularized inverse covariance matrix
#'
#' Stabilizes inversion of a covariance matrix by converting it to a correlation
#' matrix, adding a nonnegative constant to its diagonal, and transforming the
#' inverse back to the original covariance scale. This procedure is useful when
#' the input matrix is ill-conditioned or nearly singular.
#'
#' @param X A symmetric numeric covariance or cross-product matrix with strictly
#'   positive diagonal entries.
#' @param reg_scale A nonnegative numeric value specifying the diagonal
#'   regularization added to the correlation matrix. Larger values produce more
#'   stable but more strongly regularized inverses. The default is \code{0.1}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{inv}{The regularized inverse covariance matrix, transformed back to
#'     the original scale of \code{X}.}
#'   \item{lambda}{The diagonal regularization value used.}
#'   \item{condition}{The approximate condition number of the regularized
#'     correlation matrix.}
#' }
#'
#' @details
#' Let \eqn{D} denote the diagonal matrix formed from the diagonal of \eqn{X},
#' and let
#'
#' \deqn{R = D^{-1/2} X D^{-1/2}}
#'
#' be the corresponding correlation matrix. The function calculates
#'
#' \deqn{D^{-1/2}(R + \lambda I)^{-1}D^{-1/2},}
#'
#' where \eqn{\lambda} is specified by \code{reg_scale}. Cholesky inversion is
#' attempted first; if it fails, the function falls back to a general matrix
#' inverse.
#'
#' @export
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

#' Test K cell-type-specific molecular-trait associations
#'
#' Performs ctOWAS association testing for K cell types by combining
#' cell-type-specific SNP prediction weights, GWAS summary statistics, and
#' genotype data from a linkage-disequilibrium reference panel. The function
#' returns both separate cell-type statistics and joint statistics adjusted for
#' correlations among the genetically predicted cell-type-specific molecular
#' traits.
#'
#' @param W A numeric P by K matrix of cell-type-specific SNP prediction
#'   weights. Rows represent SNPs, columns represent cell types, and the SNP
#'   ordering must match \code{x_g} and \code{gwas_z_score}.
#' @param gwas_z_score {Beta/Se_Beta}{A numeric vector of length P containing GWAS z score, calculated
#' from effect-size estimates over standard errors.}
#' @param x_g A numeric reference-panel genotype dosage matrix with individuals
#'   in rows and the same P SNPs in columns.
#' @param n0 Number of controls in the GWAS. This argument must be positive when
#'   \code{family = "binomial"} and is ignored for a Gaussian outcome.
#' @param n1 Number of cases in the GWAS. This argument must be positive when
#'   \code{family = "binomial"} and is ignored for a Gaussian outcome.
#' @param family Outcome family. Use \code{"binomial"} for a case-control trait
#'   or \code{"gaussian"} for a continuous trait.
#' @param regularization Method used to select the diagonal regularization
#'   scale. Use \code{"fixed"} to use \code{reg_scale}, or \code{"estimate"} to
#'   estimate a scale targeting the condition number specified by
#'   \code{condition_max}.
#' @param reg_scale Nonnegative fixed regularization scale used when
#'   \code{regularization = "fixed"}. The default is \code{0.1}.
#' @param condition_max Target upper bound for the condition number when
#'   \code{regularization = "estimate"}. The default is \code{100}.
#' @param reg_min_scale Minimum permitted estimated regularization scale.
#' @param reg_max_scale Maximum permitted estimated regularization scale.
#'
#' @details
#' For each cell type, the function first calculates a separate association
#' statistic using the GWAS SNP statistics and the variance of the genetically
#' predicted molecular trait in the reference panel.
#'
#' The predicted traits are calculated as
#'
#' \deqn{\widehat{Y} = X_g W,}
#'
#' where \eqn{X_g} is the reference genotype matrix and \eqn{W} is the
#' cell-type-specific weight matrix. Correlations among the columns of
#' \eqn{\widehat{Y}} are used to adjust the K association statistics.
#'
#' For a binary outcome, the joint model additionally incorporates an intercept
#' statistic derived from the case-control ratio. For a continuous outcome, the
#' joint model contains only the K predicted cell-type-specific traits.
#'
#' If the predicted traits have invalid variances or are nearly perfectly
#' collinear, the function returns the separate statistics instead of attempting
#' joint adjustment. Otherwise, it regularizes and inverts their covariance
#' matrix to obtain the joint cell-type statistics. The joint cell-type
#' p-values are combined using the aggregated Cauchy association test.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Z_join}{A numeric vector of length K containing the joint,
#'     correlation-adjusted cell-type Z-statistics. When joint adjustment cannot
#'     be performed, this contains the separate Z-statistics.}
#'   \item{p_join_vec}{A numeric vector of length K containing p-values
#'     corresponding to \code{Z_join}.}
#'   \item{p_join}{The ACAT-combined p-value across the K cell types.}
#'   \item{Z_sep}{A numeric vector of length K containing the separate
#'     cell-type Z-statistics.}
#'   \item{p_sep}{A numeric vector of length K containing the corresponding
#'     separate p-values.}
#'   \item{reg_scale_selected}{The regularization scale used for joint
#'     adjustment. This is \code{NA} when joint adjustment is not performed.}
#'   \item{reg_condition}{The condition number after regularization. This is
#'     \code{NA} when joint adjustment is not performed.}
#'   \item{mode}{A character string describing the analysis performed:
#'     \code{"joint"}, \code{"invalid_variance_separate"}, or
#'     \code{"collinear_separate"}.}
#' }
#'
#' @export
#'
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

#' Safely combine p-values with ACAT
#'
#' This helper wraps \code{ACAT::ACAT()} and returns \code{NA} instead of
#' stopping if ACAT fails because of numerical or input issues.
#'
#' @param p_values Numeric vector of p-values to combine.
#'
#' @return A single combined p-value, or \code{NA} if ACAT fails.
#'
#' @importFrom ACAT ACAT
#' @export
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
