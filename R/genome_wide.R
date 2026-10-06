#' Run a single-SNP GWAS scan
#'
#' @description
#' Performs univariate association testing of each SNP (column) in \code{X}
#' with a phenotype, optionally adjusting for covariates. Supports
#' quantitative traits (Gaussian) and case-control traits (binomial).
#'
#' @param X Numeric genotype matrix (or object coercible with
#'   \code{as.matrix()}) with samples in rows and SNPs in columns. Should not
#'   contain missing values (see Details).
#' @param D Numeric phenotype vector of length \code{nrow(X)}. For
#'   \code{family0 = binomial()}, coded 0 (control) and 1 (case). Samples
#'   with non-finite values are removed.
#' @param covar Optional matrix or data frame of covariates with
#'   \code{nrow(X)} rows. Character columns are trimmed and converted to
#'   factors; factors are dummy-coded with \code{model.matrix()}. Samples with
#'   any missing covariate are removed, and constant covariates are dropped.
#'   An intercept is always included. Default \code{NULL} (intercept only).
#' @param family0 A family object: \code{stats::gaussian()} (default) or
#'   \code{stats::binomial()}. Must be a family object, not a character
#'   string.
#' @param method_binomial Method used when \code{family0 = binomial()}:
#'   \code{"score"} (default) or \code{"glm"}. The score test fits the null
#'   model once and is much faster for many SNPs, but returns only z-scores
#'   and p-values (\code{Beta} and \code{se_Beta} are \code{NA}). The
#'   \code{"glm"} option fits a separate logistic regression for each SNP and
#'   returns Wald estimates. Ignored for Gaussian traits.
#'
#' @details
#' \strong{Sample filtering.} Samples with a non-finite phenotype, or with
#' any missing covariate, are removed before testing. Missing genotypes are
#' not filtered: for Gaussian traits and the score test, a SNP with any
#' missing genotype returns \code{NA}; for \code{"glm"}, the affected
#' samples are dropped for that SNP only. Impute or filter genotypes
#' beforehand for consistent sample sets.
#'
#' \strong{Gaussian traits.} The phenotype and each SNP are residualized on
#' the covariate matrix \eqn{Z} (including the intercept), and the effect is
#' \eqn{\hat\beta_j = \tilde{x}_j^\top \tilde{y} / \tilde{x}_j^\top
#' \tilde{x}_j}. The residual variance is estimated once from the null
#' (covariates-only) model with \eqn{n - q - 1} degrees of freedom, where
#' \eqn{q = \mathrm{ncol}(Z)}, and \eqn{\mathrm{SE}(\hat\beta_j) =
#' \sqrt{\hat\sigma^2 / \tilde{x}_j^\top \tilde{x}_j}}. P-values use the
#' normal approximation rather than the t distribution, which is accurate
#' for large samples.
#'
#' \strong{Binomial traits, \code{"glm"}.} Each SNP is tested by Wald test
#' in \code{glm(D ~ x + covariates, family = family0)}. SNPs whose fit fails
#' return \code{NA}.
#'
#' \strong{Binomial traits, \code{"score"}.} The null model is fitted once
#' with \code{glm.fit()}, giving fitted probabilities \eqn{\hat\mu_i} and
#' weights \eqn{w_i = \hat\mu_i(1 - \hat\mu_i)}. Each SNP is projected off the
#' covariates in the weighted metric, giving \eqn{\tilde{x}_j}, and
#'
#' \deqn{Z_j = \frac{\sum_i \tilde{x}_{ij}(D_i - \hat\mu_i)}
#'                  {\sqrt{\sum_i w_i \tilde{x}_{ij}^2}}.}
#'
#' These weights assume the canonical logit link, so the score test should
#' be used only with \code{binomial(link = "logit")}.
#'
#' Monomorphic SNPs (zero variance after covariate adjustment) return
#' \code{NA} or \code{NaN} statistics. All p-values are two-sided.
#'
#' @return A data frame with one row per SNP, in the column order of
#'   \code{X}, and columns:
#' \describe{
#'   \item{Beta}{Estimated SNP effect (log odds ratio for binomial
#'     \code{"glm"}); \code{NA} for the score test.}
#'   \item{se_Beta}{Standard error of \code{Beta}; \code{NA} for the score
#'     test.}
#'   \item{Z}{Z-statistic.}
#'   \item{P}{Two-sided p-value.}
#' }
#'
#' @seealso \code{\link{ctOWAS_assoc_test_K}}, which takes the \code{Z}
#'   column as \code{gwas_z_score}.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' p <- 10
#' X <- matrix(rbinom(n * p, 2, 0.3), n, p)
#' covar <- data.frame(age = rnorm(n, 50, 10),
#'                     sex = sample(c("F", "M"), n, replace = TRUE))
#'
#' # Quantitative trait
#' D <- 0.3 * X[, 1] + rnorm(n)
#' res_q <- run_gwas(X, D, covar = covar, family0 = gaussian())
#' head(res_q)
#'
#' # Case-control trait, score test
#' D_bin <- rbinom(n, 1, plogis(-1 + 0.5 * X[, 1]))
#' res_b <- run_gwas(X, D_bin, covar = covar, family0 = binomial(),
#'                   method_binomial = "score")
#' head(res_b)
#'
#' @importFrom stats complete.cases model.matrix lm.fit glm glm.fit pnorm
#'   gaussian binomial
#' @importFrom utils type.convert
#' @export
run_gwas <- function(X, D, covar = NULL, family0 = stats::gaussian(), method_binomial = c("score", "glm")) {

  X <- as.matrix(X)
  D <- as.numeric(D)
  n <- nrow(X)
  p <- ncol(X)

  if (length(D) != n) {stop("Length of D must equal nrow(X).")}

  # clean based on covariates
  if (!is.null(covar)) {

    Z_df <- as.data.frame(covar, stringsAsFactors = FALSE, check.names = FALSE)

    if (nrow(Z_df) != n) stop("Covariates must have the same number of rows as X.")

    Z_df[] <- lapply(Z_df, function(x) {
      if (is.character(x)) x <- trimws(x)
      type.convert(x, as.is = FALSE)
    })

    # Remove subjects with missing phenotype or covariates
    valid_rows <- is.finite(D) & complete.cases(Z_df)

    X <- X[valid_rows, , drop = FALSE]
    D <- D[valid_rows]
    Z_df <- Z_df[valid_rows, , drop = FALSE]
    n <- nrow(X)

    # Remove constant covariates
    keep_columns <- vapply(Z_df, function(x) length(unique(x)) > 1L, logical(1))
    Z_df <- Z_df[, keep_columns, drop = FALSE]

    # Convert categorical variables to numeric dummy variables
    Z <- if (ncol(Z_df)) stats::model.matrix(~ ., data = Z_df) else matrix(1, nrow = n, ncol = 1L)

    if (ncol(Z) == 1L) colnames(Z) <- "Intercept"

    storage.mode(Z) <- "double"

    if (anyNA(Z) || any(!is.finite(Z))) stop("Invalid values remain in the covariate model matrix.")

  } else {

    valid_rows <- is.finite(D)

    X <- X[valid_rows, , drop = FALSE]
    D <- D[valid_rows]
    n <- nrow(X)

    Z <- matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "Intercept"))
  }

  is_gaussian <- identical(family0$family, "gaussian")
  is_binomial <- identical(family0$family, "binomial")

  if (!is_gaussian && !is_binomial) {stop("Only gaussian() and binomial() are supported.")}

  # =============================
  # Gaussian Trait (Fast)
  # =============================
  if (is_gaussian) {

    fit_y <- stats::lm.fit(Z, D)
    ry <- fit_y$residuals

    ZtZ_inv <- solve(crossprod(Z))
    Px <- Z %*% (ZtZ_inv %*% crossprod(Z, X))
    Xr <- X - Px

    XtX <- colSums(Xr * Xr)
    XtY <- as.numeric(crossprod(Xr, ry))

    beta <- XtY / XtX
    df <- n - ncol(Z) - 1
    sigma2 <- sum(ry^2) / max(df, 1)
    se <- sqrt(sigma2 / XtX)

    z <- beta / se
    pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

    return(data.frame(
      Beta = beta,
      se_Beta = se,
      Z = z,
      P = pval
    ))
  }

  # =============================
  # Binomial Trait
  # =============================
  method_binomial <- match.arg(method_binomial)

  if (method_binomial == "glm") {

    results <- matrix(NA_real_, p, 4)

    for (j in seq_len(p)) {

      xj <- X[, j]
      dat <- data.frame(D = D, x = xj)

      if (ncol(Z) > 1) {
        dat <- cbind(dat, Z[, -1, drop = FALSE])
      }

      fit <- tryCatch(
        stats::glm(D ~ ., data = dat, family = family0),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        cs <- summary(fit)$coefficients
        if ("x" %in% rownames(cs)) {
          b <- cs["x", 1]
          se <- cs["x", 2]
          z <- b / se
          pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
          results[j, ] <- c(b, se, z, pval)
        }
      }
    }

    colnames(results) <- c("Beta", "se_Beta", "Z", "P")
    return(as.data.frame(results))
  }

  # =============================
  # Binomial Score Test (Fast)
  # =============================
  if (method_binomial == "score") {
    null_fit <- stats::glm.fit(x = Z, y = D, family = family0)
    mu <- null_fit$fitted.values
    W <- mu * (1 - mu)
    r <- D - mu

    ZtWZ_inv <- solve(crossprod(Z, Z * W))

    Zmat <- Z
    results <- matrix(NA_real_, p, 2)

    for (j in seq_len(p)) {

      x <- X[, j]
      ZtWx <- crossprod(Zmat, W * x)
      adj <- Zmat %*% (ZtWZ_inv %*% ZtWx)
      xt <- x - adj

      U <- sum(xt * r)
      V <- sum(W * xt * xt)

      if (is.finite(V) && V > 0) {
        z <- U / sqrt(V)
        pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        results[j, ] <- c(z, pval)
      }
    }
    colnames(results) <- c("Z", "P")

    return(data.frame(Beta = NA_real_,se_Beta = NA_real_, Z = results[, 1],P = results[, 2]))
  }


}
