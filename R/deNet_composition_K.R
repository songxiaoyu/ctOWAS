#' Estimate Prior Cell-Type Composition from xCell Scores
#'
#' Converts xCell enrichment scores into sample-specific prior cell-type
#' proportions with specified mean proportions and lower/upper boundaries.
#'
#' @param xCellScore Numeric matrix of xCell enrichment scores, with xCell cell
#'   types in rows and samples in columns.
#' @param cell_types Named list defining the target cell types. Each element
#'   contains one or more row names from `xCellScore`. When multiple xCell types
#'   are provided for one target cell type, their scores are averaged.
#' @param mu Numeric vector of target mean proportions. For multiple cell types,
#'   its length must equal `length(cell_types)` and its values must sum to one.
#'   When only one cell type is supplied, `mu` is its target mean proportion,
#'   and an `"Others"` compartment with mean `1 - mu` is created automatically.
#' @param lower Numeric lower boundary for the estimated proportions. It may be
#'   a single value or one value per cell type. When only one cell type is
#'   supplied, this is the lower boundary for that cell type, and the upper
#'   boundary for `"Others"` is set to `1 - lower`.
#' @param upper Numeric upper boundary for the estimated proportions. It may be
#'   a single value or one value per cell type. When only one cell type is
#'   supplied, this is the upper boundary for that cell type, and the lower
#'   boundary for `"Others"` is set to `1 - upper`.
#' @param tol Numeric convergence tolerance. Default is `1e-8`.
#' @param max_iter Maximum number of calibration iterations. Default is `10000`.
#' @param eps Small positive constant added to the xCell scores to prevent zero
#'   values and division by zero. Default is `1e-8`.
#'
#' @return A numeric sample-by-cell-type matrix. Each row sums to one, every
#'   estimated proportion satisfies its specified boundaries, and the column
#'   means approximately equal `mu`. For a single target cell type, the output
#'   contains the target cell type and `"Others"`.
#'
#' @details
#' xCell scores are relative enrichment scores, not direct cell-type
#' proportions. This function calibrates their sample-to-sample variation using
#' externally specified mean proportions.
#'
#' For each target cell type, the mean xCell score is calculated across its
#' specified xCell populations. Negative scores are truncated to zero.
#' Iterative proportional fitting and bounded normalization are then used to:
#'
#' \enumerate{
#'   \item match the specified mean composition across samples;
#'   \item make the proportions sum to one within each sample; and
#'   \item constrain each proportion between `lower` and `upper`.
#' }
#'
#' The bounds and target means must be mutually compatible. Specifically, each
#' target mean must lie within its corresponding bounds, the lower bounds must
#' sum to no more than one, and the upper bounds must sum to at least one.
#'
#' @examples
#' ## Three cell types
#' cell_types <- list(
#'   Epithelial = "Epithelial cells",
#'   `Adipocytes/Preadipocytes` = c("Adipocytes", "Preadipocytes"),
#'   Fibroblasts = "Fibroblasts"
#' )
#'
#' pi_prior <- estimate_prior(
#'   xCellScore = xCellScore,
#'   cell_types = cell_types,
#'   mu = c(0.37, 0.54, 0.09),
#'   lower = c(0.10, 0.20, 0.02),
#'   upper = c(0.70, 0.85, 0.30)
#' )
#'
#' colMeans(pi_prior)
#' apply(pi_prior, 2, range)
#' rowSums(pi_prior)
#'
#' ## One target cell type versus all other cells
#' pi_epithelial <- estimate_prior(
#'   xCellScore = xCellScore,
#'   cell_types = list(Epithelial = "Epithelial cells"),
#'   mu = 0.37,
#'   lower = 0.10,
#'   upper = 0.70
#' )
#'
#' colMeans(pi_epithelial)
#' apply(pi_epithelial, 2, range)
#' rowSums(pi_epithelial)
#'
#' @export
estimate_prior <- function(xCellScore, cell_types, mu, lower = 0.01, upper = 0.98,
                           tol = 1e-8, max_iter = 10000, eps = 1e-8) {

  K <- length(cell_types)

  if (K == 1L) {
    stopifnot(length(mu) == 1L, mu > 0, mu < 1,
              length(lower) == 1L, length(upper) == 1L,
              lower <= mu, mu <= upper)

    mu <- c(mu, 1 - mu)
    lower <- c(lower, 1 - upper)
    upper <- c(upper, 1 - lower[1])
  } else {
    lower <- rep_len(lower, K)
    upper <- rep_len(upper, K)

    stopifnot(length(mu) == K, abs(sum(mu) - 1) < tol, all(mu > 0),
              all(lower <= mu), all(mu <= upper),
              sum(lower) <= 1, sum(upper) >= 1)
  }

  missing <- setdiff(unlist(cell_types), rownames(xCellScore))
  if (length(missing))
    stop("Missing xCell types: ", paste(missing, collapse = ", "))

  P <- do.call(rbind, lapply(cell_types, function(x)
    colMeans(xCellScore[x, , drop = FALSE], na.rm = TRUE)))

  nm <- names(cell_types)
  if (is.null(nm)) nm <- paste0("Cell", seq_along(cell_types))
  rownames(P) <- nm

  # For a single specified type, add its complementary compartment
  if (K == 1L) {
    P <- rbind(P, Others = 1)
    K <- 2L
  }

  P <- pmax(P, 0) + eps

  normalize_bounded <- function(x) {
    lo <- 0
    hi <- 1

    while (sum(pmin(pmax(hi * x, lower), upper)) < 1)
      hi <- hi * 2

    for (j in seq_len(100)) {
      mid <- (lo + hi) / 2
      if (sum(pmin(pmax(mid * x, lower), upper)) < 1)
        lo <- mid
      else
        hi <- mid
    }

    pmin(pmax(hi * x, lower), upper)
  }

  error <- Inf

  for (i in seq_len(max_iter)) {
    P <- P * (ncol(P) * mu / rowSums(P))
    P <- apply(P, 2, normalize_bounded)

    error <- max(abs(rowMeans(P) - mu), abs(colSums(P) - 1))
    if (error < tol) break
  }

  if (error >= tol)
    warning("Did not fully converge; final error = ", signif(error, 3))

  t(P)
}



#' Refine cell-type composition estimates using a K-component EM algorithm
#'
#' Refines prior sample-level cell-type composition estimates using bulk
#' transcriptomic data. Each gene is represented as a weighted mixture of
#' K latent cell-type-specific expression components.
#'
#' @param expr Numeric N by G matrix of bulk expression values, with samples
#'   in rows and genes in columns.
#' @param composition Numeric N by K matrix containing initial cell-type
#'   composition estimates. Each row must sum to one.
#' @param max_outer Maximum number of outer EM iterations.
#' @param max_inner Maximum number of gene-level EM iterations.
#' @param tol_outer Convergence tolerance for the refined compositions.
#' @param tol_inner Convergence tolerance for gene-level parameters.
#' @param min_variance Lower bound for component variances.
#' @param verbose Logical; print iteration information.
#'
#' @return A list containing:
#' \describe{
#'   \item{Ecomposition}{Refined N by K cell-type composition matrix.}
#'   \item{paraM}{Array of dimension G by K by 2 containing the fitted
#'     component means and variances.}
#'   \item{concentration}{Estimated Dirichlet concentration parameter.}
#'   \item{convergence}{Character string describing convergence.}
#'   \item{iterations}{Number of outer iterations performed.}
#'   \item{time}{Processing time reported by \code{proc.time}.}
#' }
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom stats optim var
#' @export






## ============================================================
## Helper functions
## ============================================================

# Single-core K-cell-type composition estimation

Estep.K <- function(mu, sigma2, X, P, P2 = P^2, min_variance = 1e-8) {
  N <- nrow(P)
  sigma_mat <- rep(sigma2, each = N)
  U <- drop(P %*% mu)
  V <- pmax(drop(P2 %*% sigma2), min_variance)
  adjustment <- P * sigma_mat
  EY <- adjustment * ((X - U) / V) + rep(mu, each = N)
  conditional_variance <- pmax(sigma_mat - adjustment^2 / V, min_variance)
  list(EY = EY, EY2 = conditional_variance + EY^2)
}

Mstep.K <- function(EY, EY2, min_variance = 1e-8) {
  mu <- colMeans(EY)
  sigma2 <- pmax(colMeans(EY2) - mu^2, min_variance)
  list(mu = mu, sigma2 = sigma2)
}

additive_logistic <- function(theta) {
  z <- c(theta, 0)
  z <- exp(z - max(z))
  z / sum(z)
}

# max_outer = 50L; max_inner = 200L;
# tol_outer = 0.01; tol_inner = 0.01;
# min_variance = 1e-8; prior_strength = 10;
# minor_power = 0.5; max_cell_weight = 3;
# marker_margin = 0.05; optim_maxit = 100L;
# verbose = TRUE
deNet_composition_K_once <- function( expr, composition, marker_list, max_outer = 50L, max_inner = 200L, tol_outer = 0.01, tol_inner = 0.01,
    min_variance = 1e-8, prior_strength = 10, minor_power = 0.5, max_cell_weight = 3, marker_margin = 0.05, optim_maxit = 100L, verbose = TRUE) {

  start_time <- proc.time()
  expr <- as.matrix(expr)                 # samples x genes
  composition <- as.matrix(composition)   # samples x cell types
  storage.mode(expr) <- storage.mode(composition) <- "double"

  # check input
  if (nrow(expr) != nrow(composition)) stop("composition and expr must have the same number of rows.")
  if (is.null(colnames(expr))) stop("expr must have gene names as column names.")

  K <- ncol(composition);  if (K < 2L) stop("At least two cell types are required.")
  if (max_outer < 1L || max_inner < 1L || optim_maxit < 1L) stop("Iteration limits must be positive.")

  if (!is.null(rownames(expr)) && !is.null(rownames(composition))) {
    if (!setequal(rownames(expr), rownames(composition)))  stop("expr and composition contain different samples.")
    composition <- composition[rownames(expr), , drop = FALSE]
  } else if (!is.null(rownames(expr))) {
    rownames(composition) <- rownames(expr)
  } else if (!is.null(rownames(composition))) {
    rownames(expr) <- rownames(composition)
  }

  if (any(!is.finite(composition)) || any(composition < 0))  stop("composition must contain finite, nonnegative values.")
  if (any(rowSums(composition) <= 0)) stop("Every composition row must have a positive sum.")
  composition <- composition / rowSums(composition)

  cell_names <- colnames(composition)
  if (is.null(cell_names)) {cell_names <- paste0("Cell", seq_len(K)); colnames(composition) <- cell_names}

  if (!is.list(marker_list) || is.null(names(marker_list)))  stop("marker_list must be a named list.")
  if (!setequal(names(marker_list), cell_names))  stop("names(marker_list) must match colnames(composition).")
  marker_list <- lapply(marker_list[cell_names], unique)


  all_markers <- unlist(marker_list, use.names = FALSE) # Retain markers assigned to exactly one modeled cell type.
  marker_frequency <- table(all_markers)
  unique_markers <- names(marker_frequency[marker_frequency == 1L])
  marker_list <- lapply(marker_list, intersect, y = unique_markers)
  marker_list <- lapply(marker_list, intersect, y = colnames(expr))

  n_markers <- lengths(marker_list)
  if (any(n_markers == 0L)) stop("No unique marker genes found for: ", paste(names(n_markers)[n_markers == 0L], collapse = ", "))
  if (verbose) message("Unique markers retained: ", paste(names(n_markers), n_markers, sep = " = ", collapse = "; "))

  gene_names <- unique(unlist(marker_list, use.names = FALSE))
  expr <- expr[, gene_names, drop = FALSE]
  marker_owner <- integer(length(gene_names))
  for (k in seq_len(K)) marker_owner[gene_names %in% marker_list[[k]]] <- k

  finite_gene <- colSums(!is.finite(expr)) == 0L
  gene_sd <- apply(expr, 2L, stats::sd)
  keep <- finite_gene & is.finite(gene_sd) & gene_sd > min_variance
  if (verbose && any(!keep)) message("Removing ", sum(!keep), " non-finite or zero-variance markers.")
  expr <- expr[, keep, drop = FALSE]
  gene_names <- gene_names[keep]
  marker_owner <- marker_owner[keep]
  if (!ncol(expr)) stop("No valid marker genes remain.")

  # Balance total marker evidence across types, then upweight minor types.
  marker_count <- tabulate(marker_owner, nbins = K)
  mean_composition <- colMeans(composition)
  cell_weight <- pmin(mean_composition^(-minor_power), max_cell_weight)
  cell_weight <- cell_weight / mean(cell_weight)
  gene_weight <- cell_weight[marker_owner] / marker_count[marker_owner]
  gene_weight <- gene_weight / mean(gene_weight)
  names(gene_weight) <- gene_names
  if (verbose)    message("Cell-type weights: ", paste(cell_names, signif(cell_weight, 3), sep = " = ", collapse = "; "))

  expr <- t(scale(expr))         # genes x samples
  if (any(!is.finite(expr))) stop("Standardized expression is non-finite.")

  # start the analysis
  N <- ncol(expr)
  G <- nrow(expr)

  prior_composition <- composition
  cur.P <- composition
  muM <- varianceM <- NULL
  converged <- FALSE

  for (j in seq_len(max_outer)) {
    P2 <- cur.P^2

    para_list <- lapply(seq_len(G), function(g) {
      x <- expr[g, ]
      target <- marker_owner[g]

      if (j == 1L) {
        mu <- drop(crossprod(cur.P, x)) / pmax(colSums(cur.P), .Machine$double.eps)
        sigma2 <- rep(max(stats::var(x), min_variance), K)
      } else {
        mu <- muM[g, ]
        sigma2 <- varianceM[g, ]
      }

      # Positive-marker constraint without swapping component labels.
      other_max <- max(mu[-target])
      if (mu[target] < other_max + marker_margin)   {mu[target] <- other_max + marker_margin}

      for (ii in seq_len(max_inner)) {
        old_mu <- mu
        old_sigma2 <- sigma2
        # E-step: conditional moments of the unobserved cell-type-level QTs.
        estep <- Estep.K(mu, sigma2, x, cur.P, P2, min_variance)
        # M-step: update component means and variances from E-step moments.
        mstep <- Mstep.K(estep$EY, estep$EY2, min_variance)
        mu <- mstep$mu
        sigma2 <- mstep$sigma2

        other_max <- max(mu[-target])
        if (mu[target] < other_max + marker_margin) {mu[target] <- other_max + marker_margin}
        if (max(abs(mu - old_mu), abs(sigma2 - old_sigma2)) < tol_inner) break
      }
      list(mu = mu, sigma2 = sigma2)
    })

    muM <- do.call(rbind, lapply(para_list, `[[`, "mu"))
    varianceM <- do.call(rbind, lapply(para_list, `[[`, "sigma2"))
    dimnames(muM) <- dimnames(varianceM) <- list(gene_names, cell_names)
    before.P <- cur.P

    new.P <- t(vapply(seq_len(N), function(n) {
      q <- prior_composition[n, ]

      evaluate <- local({
        previous_theta <- previous_result <- NULL
        function(theta) {
          if (!is.null(previous_theta) && identical(theta, previous_theta))
            return(previous_result)

          p <- additive_logistic(theta)
          U <- drop(muM %*% p)
          V <- pmax(drop(varianceM %*% p^2), min_variance)
          residual <- expr[, n] - U
          # Weighted Gaussian negative log-likelihood, omitting constants.
          data_nll <- 0.5 * sum(gene_weight * (log(V) + residual^2 / V))
          # Cross-entropy penalty:
          # -lambda_pi * sum_k pi_ik^(0) * log(pi_ik).
          # Equivalent to lambda_pi * KL(pi_i^(0) || pi_i), up to a fixed constant.
          # prior_strength defaults to 10; it is not estimated.
          prior_nll <- -prior_strength * sum(q * log(pmax(p, 1e-12)))
          value <- data_nll + prior_nll

          weight_mu <- gene_weight * (-residual / V)
          weight_variance <- gene_weight * (1 / V - residual^2 / V^2)
          gradient_p <- drop(crossprod(muM, weight_mu)) +
            p * drop(crossprod(varianceM, weight_variance)) -
            prior_strength * q / pmax(p, 1e-12)

          weighted_gradient <- sum(p * gradient_p)
          gradient_theta <- p[-K] * (gradient_p[-K] - weighted_gradient)

          if (!is.finite(value) || any(!is.finite(gradient_theta))) {
            value <- .Machine$double.xmax
            gradient_theta <- rep(0, K - 1L)
          }
          previous_theta <<- theta
          previous_result <<- list(value = value, gradient = gradient_theta)
          previous_result
        }
      })

      theta0 <- log(pmax(cur.P[n, -K], 1e-8) / pmax(cur.P[n, K], 1e-8))
      fit <- try(stats::optim(
        theta0,
        fn = function(theta) evaluate(theta)$value,
        gr = function(theta) evaluate(theta)$gradient,
        method = "BFGS",
        control = list(maxit = optim_maxit, reltol = 1e-6)
      ), silent = TRUE)

      if (inherits(fit, "try-error") || any(!is.finite(fit$par))) cur.P[n, ]
      else additive_logistic(fit$par)
    }, numeric(K)))

    dimnames(new.P) <- list(rownames(composition), cell_names)
    cur.P <- new.P
    mean_change <- mean(abs(cur.P - before.P))
    max_change <- max(abs(cur.P - before.P))

    if (verbose)
      message("Iteration ", j, ": mean change = ", signif(mean_change, 4),
              "; maximum change = ", signif(max_change, 4))
    if (max_change < tol_outer) {
      converged <- TRUE
      break
    }
  }

  paraM <- array(NA_real_, c(G, K, 2L),
                 dimnames = list(gene_names, cell_names, c("mean", "variance")))
  paraM[, , "mean"] <- muM
  paraM[, , "variance"] <- varianceM

  list(
    Ecomposition = cur.P,
    paraM = paraM,
    marker_owner = setNames(cell_names[marker_owner], gene_names),
    gene_weight = gene_weight,
    cell_weight = setNames(cell_weight, cell_names),
    convergence = if (converged) "P converged" else "P not converged",
    iterations = j,
    time = proc.time() - start_time
  )
}

# Bootstrap aggregation. Input expr is genes x samples; output is samples x types.
pi_estimation_K <- function(
    expr, n_iteration, prior, marker_list, seed = 1,
    sample_fraction = 0.8, trim = 0.05, ...) {

  set.seed(seed)
  expr <- as.matrix(expr)
  prior <- as.matrix(prior)
  if (is.null(colnames(expr))) stop("expr must have sample column names.")

  if (!is.null(rownames(prior))) {
    if (!all(colnames(expr) %in% rownames(prior)))
      stop("Some expression samples are missing from prior.")
    prior <- prior[colnames(expr), , drop = FALSE]
  } else {
    if (ncol(expr) != nrow(prior)) stop("ncol(expr) must equal nrow(prior).")
    rownames(prior) <- colnames(expr)
  }

  estimates <- vector("list", n_iteration)
  n_sample <- round(ncol(expr) * sample_fraction)

  for (i in seq_len(n_iteration)) {
    index <- sample(seq_len(ncol(expr)), n_sample, replace = FALSE)
    expr_i <- expr[, index, drop = FALSE]
    keep <- apply(expr_i, 1L, function(x) all(is.finite(x)) && sd(x) > 1e-8)
    expr_i <- expr_i[keep, , drop = FALSE]

    fit <- try(deNet_composition_K_once(
      expr = t(expr_i),
      composition = prior[index, , drop = FALSE],
      marker_list = marker_list,
      ...
    ), silent = TRUE)

    if (inherits(fit, "try-error")) {
      warning("Iteration ", i, " failed: ", as.character(fit))
      next
    }
    estimates[[i]] <- fit$Ecomposition
    message("Finished iteration ", i)
  }

  estimates <- estimates[!vapply(estimates, is.null, logical(1))]
  if (!length(estimates))
    stop("deNet_composition_K_once() failed in every iteration.")

  result <- matrix(NA_real_, nrow = ncol(expr), ncol = ncol(prior),
                   dimnames = list(colnames(expr), colnames(prior)))
  for (s in rownames(result)) {
    for (k in colnames(result)) {
      values <- vapply(estimates, function(P)
        if (s %in% rownames(P)) P[s, k] else NA_real_, numeric(1))
      values <- values[is.finite(values)]
      if (length(values)) result[s, k] <- mean(values, trim = trim)
    }
  }

  valid <- complete.cases(result) & rowSums(result) > 0
  result[valid, ] <- result[valid, , drop = FALSE] /
    rowSums(result[valid, , drop = FALSE])
  result
}



pi_estimation_K <- function(expr, n_iteration, prior, marker_list,
                            seed = 1, sample_fraction = 0.8,
                            trim = 0.05, ...) {
  set.seed(seed)

  expr <- as.matrix(expr)
  prior <- as.matrix(prior)

  stopifnot(ncol(expr) == nrow(prior))

  if (is.null(colnames(expr))) stop("expr must have sample names as column names.")

  if (!is.null(rownames(prior))) {
    prior <- prior[colnames(expr), , drop = FALSE]
  } else { rownames(prior) <- colnames(expr)}

  if (is.null(colnames(prior)))  colnames(prior) <- paste0("Cell", seq_len(ncol(prior)))

  estimates <- vector("list", n_iteration)
  n_sample <- round(ncol(expr) * sample_fraction)

  for (i in seq_len(n_iteration)) {

    sample_index <- sample(seq_len(ncol(expr)), n_sample, replace = FALSE)

    fit <- try(
      deNet_composition_K_once(
        expr = t(expr[, sample_index, drop = FALSE]),
        composition = prior[sample_index, , drop = FALSE],
        marker_list=marker_list,
        ...
      ),
      silent = F
    )

    if (inherits(fit, "try-error")) {
      warning("Iteration ", i, " failed.")
      next
    }

    estimates[[i]] <- fit$Ecomposition
    message("Finished iteration ", i)
  }

  estimates <- estimates[!vapply(estimates, is.null, logical(1))]

  if (!length(estimates))
    stop("deNet_composition_K_once() failed in every iteration.")

  sample_names <- colnames(expr)
  cell_names <- colnames(prior)

  result <- matrix(
    NA_real_,
    nrow = length(sample_names),
    ncol = length(cell_names),
    dimnames = list(sample_names, cell_names)
  )

  for (s in sample_names) {
    for (k in cell_names) {
      values <- vapply(estimates, function(P) {
        if (s %in% rownames(P)) P[s, k] else NA_real_
      }, numeric(1))

      values <- values[is.finite(values)]

      if (length(values))
        result[s, k] <- mean(values, trim = trim)
    }
  }

  # Renormalize each sample to sum to one
  valid <- complete.cases(result) & rowSums(result) > 0
  result[valid, ] <- result[valid, , drop = FALSE] / rowSums(result[valid, , drop = FALSE])

  result
}
