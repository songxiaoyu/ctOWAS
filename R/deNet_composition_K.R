#' Estimate Prior Cell-Type Composition from xCell Scores
#'
#' Converts xCell enrichment scores into sample-specific prior cell-type
#' proportions with user-specified mean proportions and lower and upper bounds.
#'
#' @param xCellScore A numeric matrix of xCell enrichment scores, with xCell
#'   cell types in rows and samples in columns. The matrix must have row names.
#' @param cell_types A list defining the target cell-type compartments. Each
#'   element contains one or more row names of `xCellScore`. When multiple xCell
#'   types are assigned to the same compartment, their scores are averaged
#'   within each sample. If the list is named, its names are used as the output
#'   column names; otherwise, names of the form `"Cell1"`, `"Cell2"`, and so on
#'   are generated.
#' @param mu A numeric vector specifying the target mean proportions. When
#'   multiple compartments are supplied, its length must equal
#'   `length(cell_types)`, all values must be positive, and the values must sum
#'   to one within `tol`. When a single compartment is supplied, `mu` must be a
#'   scalar strictly between zero and one; a complementary `"Others"`
#'   compartment with target mean `1 - mu` is added automatically.
#' @param lower Numeric lower bounds for the estimated proportions. For multiple
#'   compartments, a scalar is applied to all compartments, and a vector is
#'   recycled to `length(cell_types)`. For a single supplied compartment,
#'   `lower` must be a scalar and the lower bound for `"Others"` is set to
#'   `1 - upper`. Default is `0.01`.
#' @param upper Numeric upper bounds for the estimated proportions. For multiple
#'   compartments, a scalar is applied to all compartments, and a vector is
#'   recycled to `length(cell_types)`. For a single supplied compartment,
#'   `upper` must be a scalar and the upper bound for `"Others"` is set to
#'   `1 - lower`. Default is `0.98`.
#' @param tol Numeric convergence tolerance for the maximum deviation of the
#'   estimated column means from `mu` and the row sums from one. Default is
#'   `1e-8`.
#' @param max_iter A positive integer specifying the maximum number of
#'   calibration iterations. Default is `10000`.
#' @param eps A small positive constant added to all nonnegative enrichment
#'   scores to prevent zero values and division by zero. Default is `1e-8`.
#'
#' @return A numeric matrix with samples in rows and cell-type compartments in
#'   columns. Each row sums to one, and each estimated proportion lies within
#'   its corresponding lower and upper bounds. Upon convergence, the column
#'   means equal the target proportions `mu` within `tol`.
#'
#'   When one target compartment is supplied, the output contains two columns:
#'   the specified compartment and `"Others"`. Their target means are `mu` and
#'   `1 - mu`, respectively.
#'
#' @details
#' xCell scores are relative enrichment scores rather than direct estimates of
#' cell-type proportions. This function uses their sample-to-sample variation
#' to construct prior composition estimates calibrated to externally specified
#' population means.
#'
#' For each target compartment and sample, the function averages the scores of
#' the corresponding xCell populations using `na.rm = TRUE`. Negative scores
#' are truncated to zero, and `eps` is added to all scores. When only one target
#' compartment is supplied, the initial score for the complementary `"Others"`
#' compartment is set to one in every sample.
#'
#' The function then alternates between:
#'
#' \enumerate{
#'   \item scaling each compartment across samples toward its target mean; and
#'   \item normalizing each sample to sum to one while enforcing the specified
#'   lower and upper bounds.
#' }
#'
#' Iteration stops when both the maximum deviation of the compartment means
#' from `mu` and the maximum deviation of the sample totals from one are less
#' than `tol`, or when `max_iter` is reached. A warning is issued if convergence
#' is not achieved.
#'
#' The target means and bounds must be mutually compatible: every target mean
#' must lie within its corresponding bounds, the lower bounds must sum to no
#' more than one, and the upper bounds must sum to at least one. All xCell types
#' listed in `cell_types` must be present in `rownames(xCellScore)`.
#'
#' @examples
#' xCellScore <- matrix(
#'   c(
#'     0.8, 0.6, 0.7, 0.9,
#'     0.4, 0.5, 0.3, 0.6,
#'     0.5, 0.4, 0.6, 0.3,
#'     0.2, 0.3, 0.2, 0.4
#'   ),
#'   nrow = 4,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("Epithelial cells", "Adipocytes",
#'       "Preadipocytes", "Fibroblasts"),
#'     paste0("Sample", 1:4)
#'   )
#' )
#'
#' ## Three target compartments
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
#' ## One target compartment versus all other cells
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
#' @description
#' Refines prior sample-level cell-type composition estimates using bulk
#' transcriptomic data and cell-type marker genes. Each marker gene is
#' modelled as a mixture of K latent cell-type-specific expression
#' components, whose means and variances are fitted by EM. Sample
#' compositions are then updated by penalised maximum likelihood, with a
#' cross-entropy penalty that shrinks them toward the prior. The two steps
#' alternate until the compositions converge.
#'
#' @param expr Numeric N by G matrix of bulk expression values, with samples
#'   in rows and genes in columns. Column names (gene IDs) are required.
#' @param composition Numeric N by K matrix of initial (prior) cell-type
#'   compositions, with samples in rows and cell types in columns. Values
#'   must be finite and nonnegative, and every row must have a positive sum.
#'   Rows are rescaled to sum to one. If both \code{expr} and
#'   \code{composition} have row names, they must contain the same samples,
#'   and \code{composition} is reordered to match \code{expr}. If column
#'   names are missing, cell types are named \code{"Cell1"},
#'   \code{"Cell2"}, and so on. At least two cell types are required.
#' @param marker_list Named list of character vectors giving the marker
#'   genes for each cell type. Names must match \code{colnames(composition)}
#'   (in any order). Only markers that are assigned to exactly one cell type
#'   and are present in \code{colnames(expr)} are used; every cell type must
#'   retain at least one such marker.
#' @param max_outer Maximum number of outer iterations (alternating
#'   gene-parameter and composition updates). Default \code{50L}.
#' @param max_inner Maximum number of gene-level EM iterations within each
#'   outer iteration. Default \code{200L}.
#' @param tol_outer Convergence tolerance for the compositions: the
#'   algorithm stops when the maximum absolute change in any composition
#'   value between outer iterations falls below this. Default \code{0.01}.
#' @param tol_inner Convergence tolerance for gene-level EM: the maximum
#'   absolute change in component means and variances. Default \code{0.01}.
#' @param min_variance Lower bound for component variances, also used as
#'   the threshold for removing zero-variance genes. Default \code{1e-8}.
#' @param prior_strength Nonnegative weight \eqn{\lambda} of the
#'   cross-entropy penalty \eqn{-\lambda \sum_k \pi^{(0)}_{ik}
#'   \log \pi_{ik}} that shrinks each sample's composition toward its prior.
#'   Larger values keep estimates closer to \code{composition}. It is fixed,
#'   not estimated. Default \code{10}.
#' @param minor_power Exponent used to upweight markers of low-abundance
#'   cell types: cell type weights are proportional to
#'   \code{mean_composition^(-minor_power)}. Use \code{0} for no
#'   upweighting. Default \code{0.5}.
#' @param max_cell_weight Upper cap on the raw cell type weight before
#'   normalisation. Default \code{3}.
#' @param marker_margin Minimum amount by which a marker's own cell-type
#'   component mean must exceed the largest mean of the other components
#'   (on the standardised scale). Enforced after every M-step. Default
#'   \code{0.05}.
#' @param optim_maxit Maximum number of BFGS iterations in each per-sample
#'   composition update. Default \code{100L}.
#' @param verbose Logical; if \code{TRUE}, print the retained marker counts,
#'   cell type weights and per-iteration composition changes. Default
#'   \code{TRUE}.
#'
#' @details
#' \strong{Gene selection and weighting.} Only marker genes are modelled.
#' Markers shared by more than one cell type, markers absent from
#' \code{expr}, and genes with non-finite values or standard deviation
#' \eqn{\le} \code{min_variance} are removed. Each remaining gene is
#' standardised across samples. Gene weights balance total marker evidence
#' across cell types (each type's weight is divided by its number of
#' markers) and upweight minor cell types via \code{minor_power} and
#' \code{max_cell_weight}; weights are normalised to mean one.
#'
#' \strong{Gene-level step.} For each gene, the bulk value in sample
#' \eqn{n} is modelled as \eqn{\sum_k \pi_{nk} Y_{nk}} with
#' \eqn{Y_{nk} \sim N(\mu_k, \sigma^2_k)}. Component means and variances
#' are fitted by EM, initialised in the first outer iteration from
#' composition-weighted means and the gene's overall variance, and warm
#' started thereafter. The marker's own component is constrained to have
#' the largest mean (by at least \code{marker_margin}).
#'
#' \strong{Composition step.} Each sample's composition is updated by
#' minimising the weighted Gaussian negative log-likelihood (mean
#' \eqn{\sum_k \pi_k \mu_k}, variance \eqn{\sum_k \pi_k^2 \sigma^2_k})
#' plus the cross-entropy prior penalty, using BFGS on an additive
#' log-ratio parameterisation so estimates stay on the simplex. If the
#' optimisation fails for a sample, its previous composition is kept.
#'
#' This function relies on the internal helpers \code{Estep.K()},
#' \code{Mstep.K()} and \code{additive_logistic()}.
#'
#' @return A list containing:
#' \describe{
#'   \item{Ecomposition}{Refined N by K cell-type composition matrix; rows
#'     sum to one, with sample and cell type names as dimnames.}
#'   \item{paraM}{Array of dimension G by K by 2 containing the fitted
#'     component means (\code{"mean"}) and variances (\code{"variance"}) on
#'     the standardised scale, for the G retained marker genes.}
#'   \item{marker_owner}{Named character vector giving the cell type each
#'     retained marker gene belongs to.}
#'   \item{gene_weight}{Named numeric vector of gene weights used in the
#'     likelihood.}
#'   \item{cell_weight}{Named numeric vector of normalised cell type
#'     weights.}
#'   \item{convergence}{\code{"P converged"} or \code{"P not converged"}.}
#'   \item{iterations}{Number of outer iterations performed.}
#'   \item{time}{Processing time reported by \code{proc.time}.}
#' }
#'
#' @seealso \code{\link{pi_estimation_K}} for repeated subsampling and
#'   averaging of these estimates.
#'
#' @examples
#' \dontrun{
#' fit <- deNet_composition_K_once(
#'   expr        = t(expr_matrix),   # samples x genes
#'   composition = prior_matrix,     # samples x cell types
#'   marker_list = list(
#'     Tcell  = c("CD3D", "CD3E"),
#'     Bcell  = c("CD79A", "MS4A1"),
#'     Mono   = c("CD14", "LYZ")
#'   ),
#'   prior_strength = 10
#' )
#' head(fit$Ecomposition)
#' fit$convergence
#' }
#'
#' @importFrom stats optim sd var setNames
#' @export
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

#' Estimate cell-type composition by repeated subsampling
#'
#' @description
#' Estimates the cell-type composition (proportions) of each sample by
#' repeatedly fitting \code{deNet_composition_K_once()} to random subsets of
#' samples, then combining the per-iteration estimates with a trimmed mean.
#' The combined proportions are renormalised so that each sample's row sums
#' to 1.
#'
#' @param expr A numeric matrix (or object coercible with \code{as.matrix()})
#'   of expression values with genes in rows and samples in columns. Column
#'   names (sample IDs) are required.
#' @param n_iteration Integer. Number of subsampling iterations to run.
#' @param prior A numeric matrix (or coercible object) of prior cell-type
#'   compositions with samples in rows and cell types in columns. If it has
#'   row names, they must include every sample in \code{colnames(expr)};
#'   rows are reordered to match \code{expr}, and extra samples are dropped.
#'   If it has no row names, it must have exactly \code{ncol(expr)} rows in
#'   the same order as the columns of \code{expr}.
#' @param marker_list A list of marker genes for each cell type, passed to
#'   \code{deNet_composition_K_once()}.
#' @param seed Integer. Random seed set with \code{set.seed()} before
#'   subsampling, so results are reproducible. Default \code{1}.
#' @param sample_fraction Numeric in (0, 1]. Fraction of samples drawn
#'   (without replacement) in each iteration; the subset size is
#'   \code{round(ncol(expr) * sample_fraction)}. Default \code{0.8}.
#' @param trim Numeric in [0, 0.5). Fraction of observations trimmed from
#'   each end when averaging a sample's estimates across iterations, passed
#'   to \code{mean(..., trim = trim)}. Default \code{0.05}.
#' @param ... Additional arguments passed to
#'   \code{deNet_composition_K_once()}.
#'
#' @details
#' In each iteration the function:
#' \enumerate{
#'   \item draws a random subset of samples;
#'   \item removes genes that contain non-finite values or have near-zero
#'     variance (standard deviation \eqn{\le 10^{-8}}) within that subset;
#'   \item calls \code{deNet_composition_K_once()} with the transposed
#'     expression subset (samples in rows, genes in columns) and the matching
#'     rows of \code{prior};
#'   \item stores the returned \code{Ecomposition} matrix.
#' }
#' Iterations in which \code{deNet_composition_K_once()} throws an error are
#' skipped with a warning; the function stops only if every iteration fails.
#'
#' For each sample and cell type, the final value is the trimmed mean of the
#' finite estimates from the iterations in which that sample was drawn.
#' Samples never drawn, or with no finite estimates, are returned as
#' \code{NA}. Rows with no missing values and a positive sum are rescaled to
#' sum to 1; other rows are returned unscaled.
#'
#' Note that \code{set.seed()} changes the global random number generator
#' state as a side effect.
#'
#' @return A numeric matrix with one row per sample (row names from
#'   \code{colnames(expr)}) and one column per cell type (column names from
#'   \code{colnames(prior)}), containing the estimated cell-type proportions.
#'
#' @seealso \code{\link{deNet_composition_K_once}}
#'
#' @examples
#' \dontrun{
#' est <- pi_estimation_K(
#'   expr            = expr_matrix,   # genes x samples
#'   n_iteration     = 50,
#'   prior           = prior_matrix,  # samples x cell types
#'   marker_list     = markers,
#'   seed            = 123,
#'   sample_fraction = 0.8,
#'   trim            = 0.05
#' )
#' head(est)
#' rowSums(est)
#' }
#'
#' @importFrom stats sd complete.cases
#' @export
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

