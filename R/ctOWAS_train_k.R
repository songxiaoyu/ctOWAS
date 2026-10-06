#' Train a K-cell-type ctOWAS prediction model
#'
#' @description
#' Trains a cell-type-specific molecular prediction model for one gene or
#' protein using genotype data, estimated cell-type fractions and optional
#' covariates. Repeated split-sample stability screening assigns each
#' predictor an empirical FDR score, which is then used as a predictor-specific
#' penalty factor in a final elastic-net fit.
#'
#' @param y A numeric vector or one-column matrix of length N containing the
#'   molecular phenotype for one gene or protein.
#' @param x A numeric N by P genotype dosage matrix, with samples in rows and
#'   cis-SNPs in columns.
#' @param cov An optional numeric N by Q covariate matrix. Column names are
#'   required when supplied. If \code{NULL} (default), no covariates are
#'   included.
#' @param pi_k A numeric N by K matrix of estimated cell-type fractions, with
#'   samples in rows (in the same order as \code{x}) and cell types in
#'   columns. Rows would typically sum to approximately one.
#' @param xNameMatrix Optional data frame of SNP annotation with P rows. The
#'   first column must contain SNP identifiers in the same order as the columns
#'   of \code{x}. If \code{NULL} (default), identifiers are taken from
#'   \code{colnames(x)}, or generated as \code{SNP1}, \code{SNP2}, and so on.
#' @param yName Optional gene or protein identifier, returned unchanged with
#'   the fitted model.
#' @param foldid Optional integer vector of length N assigning samples to
#'   cross-validation folds. If \code{NULL} (default), random 10-fold
#'   assignments are generated.
#' @param No_split Number of splits for evaulation FDR. A small integer should
#' give fairly robust results. Default No_split = 11.
#' @param alpha Elastic-net mixing parameter passed to \code{glmnet}; \code{0}
#'   is ridge and \code{1} is lasso. Default \code{0.5}.
#' @param seed Optional integer random seed. It is applied with
#'   \code{set.seed()} only when \code{foldid} is \code{NULL}, before the
#'   folds are generated, and so also fixes the subsequent sample splits.
#'   When \code{foldid} is supplied, \code{seed} is ignored and the splits use
#'   the current random-number state. With \code{foldid = NULL} and
#'   \code{seed = NULL}, the generator is re-initialised randomly.
#'
#' @details
#' \strong{Design matrix.} For K cell types, the design matrix contains a
#' separate intercept and a separate set of P SNP coefficients for each cell
#' type. For sample \eqn{i} and cell type \eqn{k}, the design block is
#'
#' \deqn{\pi_{ik}(1, x_{i1}, \ldots, x_{iP}),}
#'
#' where \eqn{\pi_{ik}} is the fraction of cell type \eqn{k} and
#' \eqn{x_{ij}} is the dosage of SNP \eqn{j}. Covariates, if any, are appended
#' as shared columns. All models are fitted with \code{intercept = FALSE},
#' because the cell-type intercepts are part of the design.
#'
#' \strong{Stability screening.} A penalty \eqn{\lambda} is first chosen by
#' cross-validated elastic net (\code{lambda.1se}). Then, in each of 11
#' iterations, samples are randomly split into two halves. In half 1, an
#' elastic-net model is fitted at \eqn{\lambda}, giving
#' \eqn{\hat\beta_1}, and a near-unpenalised fit (ridge,
#' \eqn{\lambda = 10^{-6}}) gives the direction estimate
#' \eqn{\tilde\beta_1}. In half 2, a near-unpenalised fit gives
#' \eqn{\hat\beta_2}. For every column of the design matrix, the agreement
#' statistic is
#'
#' \deqn{M = \mathrm{sign}(\tilde\beta_1)\,\mathrm{sign}(\hat\beta_2)
#'           \left(|\hat\beta_1| + |\hat\beta_2|\right).}
#'
#' Positive values indicate directionally consistent effects. For a
#' predictor with \eqn{M = m > 0}, the empirical FDR is
#'
#' \deqn{\widehat{\mathrm{FDR}}(m) =
#'   \frac{\#\{M_j < 0,\ |M_j| \ge m\}}{\max(\#\{M_j \ge m\}, 1)},}
#'
#' and predictors with \eqn{M \le 0} receive an FDR of 1. FDR values are
#' averaged across the 11 splits and capped at 1.
#'
#' \strong{Final model.} The averaged FDR values are used as
#' \code{penalty.factor} in a cross-validated elastic net (\code{lambda.1se}),
#' so stable predictors are penalised less and predictors with FDR of 1
#' receive the full penalty. There is no hard FDR threshold: selection is
#' done by the penalised fit itself. The FDR-based penalty applies to every
#' column, including cell-type intercepts and covariates.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{intercept}{Numeric vector of length K with the fitted
#'     cell-type-specific intercepts.}
#'   \item{W}{Numeric P by K matrix of cell-type-specific SNP weights, with SNP
#'     identifiers as row names and \code{Cell1}, \ldots, \code{CellK} as
#'     column names. Unselected coefficients are zero.}
#'   \item{beta.all.models}{Coefficient matrix with K columns and
#'     \eqn{1 + P + Q} rows: the intercept, the P SNP coefficients and, when
#'     supplied, the Q covariate coefficients. Covariate coefficients are
#'     shared, so they are repeated in every column.}
#'   \item{yName}{The identifier supplied through \code{yName}.}
#'   \item{xNameMatrix}{The SNP annotation supplied to, or generated by, the
#'     function.}
#'   \item{fdr.mean}{Numeric vector of length \eqn{K(P + 1) + Q} with the
#'     averaged empirical FDR for each design column (ordered by cell type:
#'     intercept then SNPs, followed by covariates), used as the penalty
#'     factors.}
#'   \item{Num_SNP_input}{The number of input SNPs, P.}
#'   \item{Num_SNP_model}{The number of unique SNPs with a nonzero weight in
#'     at least one cell type.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- ctOWAS_train_K(
#'   y           = expr_gene,        # length N
#'   x           = geno,             # N x P dosages
#'   cov         = covariates,       # N x Q, or NULL
#'   pi_k        = cell_fractions,   # N x K
#'   yName       = "GENE1",
#'   alpha       = 0.5,
#'   seed        = 2026
#' )
#' fit$Num_SNP_model
#' head(fit$W)
#' }
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @export

ctOWAS_train_K <- function(y, x, cov = NULL, pi_k,xNameMatrix = NULL, yName = NULL,
                            foldid = NULL,  No_split=11, alpha=0.5, seed=NULL) {

  ## ---- Input checks and setup ----
  y    <- as.matrix(y)
  x    <- as.matrix(x)
  pi_k <- as.matrix(pi_k)

  n <- nrow(x)
  p <- ncol(x)
  K <- ncol(pi_k)
  if (is.null(cov)) {pcov <- 0L} else {pcov <- ncol(cov)}

  if (nrow(pi_k) != n) { stop("pi_k must have the same number of rows as x (samples).")}

  # Use supplied SNP annotation, or create a simple fallback label set.
  if (is.null(xNameMatrix)) {
    if (!is.null(colnames(x))) {xNameMatrix <- data.frame(SNP = colnames(x), stringsAsFactors = FALSE)} else {xNameMatrix <- data.frame(SNP = paste0("SNP", seq_len(p)),stringsAsFactors = FALSE)}
  }

  # fix randomness

  if (is.null(foldid)) {  set.seed(seed); foldid <- sample(rep(1:10, length.out = n))}


  ## ---- Create cell-type-level model design matrix ----
  XK <- do.call(cbind, lapply(seq_len(K), \(i) cbind(1, x) * pi_k[, i]))
  if (is.null(cov)) {XX <- XK} else {cov=as.matrix(cov); XX <- cbind(XK, cov)}

  # ------- add random split of the data into two components to identify important predictors ----------

  ft11 <- glmnet::cv.glmnet(x = XX, y = y,family= "gaussian",foldid= foldid, alpha = alpha,intercept = FALSE)
  lambda_cell <- ft11$lambda.1se

  fdr=NULL
  for (split in 1:No_split) {

    idx=sample(rep(1:2, length.out = n))
    ft1 <- glmnet::glmnet(
      x = XX[idx==1,], y = y[idx==1,],
      family         = "gaussian",
      lambda         = lambda_cell,
      alpha          = alpha,
      intercept = FALSE
    )
    ft1d <- glmnet::glmnet(
      x = XX[idx==1,], y = y[idx==1,],
      family         = "gaussian",
      lambda         = 0.000001,
      alpha          = 0,
      intercept = FALSE
    ) # direction
    ft2 <- glmnet::glmnet(
      x = XX[idx==2,], y = y[idx==2,],
      family         = "gaussian",
      lambda         = 0.000001,
      alpha          = 0,
      intercept = FALSE
    ) # OLS

    M1 = sign(ft1d$beta)*sign(ft2$beta)*(abs(ft1$beta) + abs(ft2$beta))
    # m1 <- as.numeric(M1)[1: ((p+1L)*K)]

    fdr1 <- vapply(M1, \(x) {
      if (x <= 0) 1
      else sum(M1 < 0 & abs(M1) >= x) / max(sum(M1 >= x),1)
    }, numeric(1))

    fdr=cbind(fdr, fdr1)
  }
  fdr_mean=apply(fdr, 1, mean)
  fdr_mean[which(fdr_mean>1)]=1
  # use fdr_mean as penalty to refit lasso

  ft_w <- glmnet::cv.glmnet(x = XX, y = y,family= "gaussian",foldid= foldid,
                            penalty.factor = fdr_mean, alpha = alpha,intercept = FALSE)
  lambda_cell <- ft_w$lambda.1se
  ft <- glmnet::glmnet( x = XX, y = y, family = "gaussian",
                        lambda = lambda_cell, alpha = alpha, penalty.factor = fdr_mean,
                        intercept = FALSE)

  ## ---- Harvest results ----
  coef_no_intercept <- as.numeric(ft$beta)
  coef_cell=matrix(coef_no_intercept[1:(K*(1L+p))], ncol=K)
  colnames(coef_cell) <- paste0("Cell", seq_len(K))
  rownames(coef_cell) <- c("intercept", xNameMatrix[[1]])

  ## Append common covariate coefficients to every cell-type column
  if (!is.null(cov)) {
    cov_coef <- coef_no_intercept[(K * (p + 1L) + 1L):ncol(XX)]
    coef_cell <- rbind(coef_cell, matrix(cov_coef, ncol = K, nrow = pcov))
    rownames(coef_cell)[(p + 2L):(p + 1L + pcov)] <- colnames(cov)
  }


  a_cell <- coef_cell[1,] # length K
  B_mat  <- coef_cell[2:(p + 1L), , drop = FALSE]

  idx_x_selected <- which(rowSums(abs(B_mat)) > 0)

  list(
    intercept       = a_cell,
    W               = B_mat,
    beta.all.models = coef_cell,
    yName           = yName,
    xNameMatrix     = xNameMatrix,
    fdr.mean        = fdr_mean,
    Num_SNP_input   = p,
    Num_SNP_model   = length(idx_x_selected)
  )

}

