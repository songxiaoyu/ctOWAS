#' Estimate K-cell-type fractions using BayesDeBulk
#'
#' Estimates the proportions of K cell types in bulk expression samples using
#' \code{BayesDeBulk}. The function constructs the marker-gene index required by
#' \code{BayesDeBulk}, fits the model, and returns both the estimated cell-type
#' fractions and cell-type-specific expression profiles.
#'
#' @param expr A numeric gene-by-sample matrix or data frame of bulk expression
#'   values. Rows represent genes and columns represent samples. Gene identifiers
#'   should be supplied as row names and sample identifiers as column names.
#' @param markers A named list of marker-gene vectors. Each list element
#'   corresponds to one cell type, and its name is used as the cell-type label.
#' @param seed An integer random seed used to ensure reproducibility.
#' @param n.iter A positive integer specifying the total number of MCMC
#'   iterations passed to \code{BayesDeBulk}.
#' @param burn.in A nonnegative integer specifying the number of burn-in
#'   iterations passed to \code{BayesDeBulk}.
#' @param ... Additional arguments passed to \code{BayesDeBulk}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{cell_fraction}{A data frame of estimated cell-type fractions, with
#'     one row per sample and one column per cell type.}
#'   \item{cell_expression}{The fitted cell-type-specific expression profiles
#'     returned by \code{BayesDeBulk}.}
#' }
#'
#' @importFrom BayesDeBulk BayesDeBulk
#' @export

pi_estimation_K <- function(expr,
                            markers,
                            seed   = 1,
                            n.iter = 10000,
                            burn.in = 1000,
                            ...) {
  cell.type<-names(markers)
  index.matrix<-NULL
  for (s in 1:length(cell.type)){
    for (k in 1:length(cell.type)){
      if (s!=k){
        mg<-match(markers[[s]],markers[[k]])
        index.matrix<-rbind(index.matrix,cbind(rep(cell.type[s],sum(is.na(mg))),
                                               rep(cell.type[k],sum(is.na(mg))),markers[[s]][is.na(mg)]))
      }
    }
  }
  set.seed(seed)

  fit <- BayesDeBulk(
    n.iter  = n.iter,
    burn.in = burn.in,
    Y       = list(expr),
    markers = index.matrix,
    ...
  )

  cell_fraction <- as.data.frame(fit$cell.fraction)

  list(
    cell_fraction  = cell_fraction,
    cell_expression = fit$cell.expression
  )
}


#' Estimate robust cell-type proportions using an ensemble of TSNet models
#'
#' Refines prior cell-type proportion estimates by repeatedly fitting
#' \code{TSNet} to random subsets containing 80 percent of the samples. The
#' resulting sample-level estimates are aggregated across successful ensemble
#' iterations using a 5-percent trimmed mean, providing estimates that are less
#' sensitive to individual subsamples and extreme fitted values. This strategy
#' in MiXcan.
#'
#' @param expr A numeric gene-by-sample matrix or data frame of bulk expression
#'   values. Rows represent genes and columns represent tissue samples. Gene
#'   identifiers should be supplied as row names and sample identifiers as
#'   column names.
#' @param n_iteration A positive integer specifying the number of TSNet ensemble
#'   iterations.
#' @param prior A numeric vector of prior cell-type proportion estimates. Its
#'   length and ordering must correspond to the columns of \code{expr}.
#' @param seed An integer random seed used to make sample subsampling
#'   reproducible.
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{sample}{Sample identifier.}
#'   \item{mean_trim_0.05}{Robust cell-type proportion estimate, calculated as
#'     the 5-percent trimmed mean across successful TSNet ensemble fits.}
#' }
#'
#' @importFrom dplyr bind_rows group_by summarise ungroup
#' @importFrom tibble tibble
#' @export

pi_estimation_2 <- function(expr, n_iteration, prior, seed=1){
  set.seed(seed)
  TSNetB_prop <- NULL
  for (i in 1:n_iteration){
    sample_index <- sample(sample(1:dim(expr)[2], round(dim(expr)[2] * 0.8), replace = FALSE))
    TSNetB=try(deNet_purity(exprM=t(expr[,sample_index]),purity=prior[sample_index]))
    if(inherits(TSNetB, "try-error")) {
      #error handling code, maybe just skip this iteration using
      next }
    TSNetB_prop_once <- tibble(sample = colnames(expr[,sample_index]), prop = TSNetB[[1]], rep = i)
    TSNetB_prop <- TSNetB_prop %>% bind_rows(TSNetB_prop_once)
    print(c("Finished iteration", i))
  }
  TSNetB_prop_trim_mean <- TSNetB_prop %>%
    group_by(sample) %>%
    summarise(mean_trim_0.05 = mean(prop, trim = 0.05)) %>%
    ungroup

  return(TSNetB_prop_trim_mean)
}



