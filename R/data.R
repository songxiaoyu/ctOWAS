#' Cell-type-specific marker gene reference
#'
#' A total of 9 marker gene sets for 3 cell types used in ctOWAS
#' for estimating cell-type proportions from bulk gene expression data.
#' Each cell type has multiple marker gene sets.
#'
#' @format A named list of 9 character vectors containing marker gene symbols
#'   for 3 cell types. Each list element represents a marker gene set,
#'   with names formatted as \code{CellType#customN}, where
#'   \code{CellType} indicates the cell type and \code{N} identifies
#'   the marker gene set.
#' @source Cell-type-specific marker gene sets compiled for ctOWAS analysis.
"markers"



#' Example bulk gene expression data
#'
#' A bulk gene expression dataset used in ctOWAS examples for
#' estimating cell-type proportions using cell-type-specific marker genes.
#'
#' @format A data frame with 24,919 genes in rows and 125 samples
#'   in columns. Row names represent gene symbols, and column names
#'   represent sample identifiers.
#' @source Example bulk gene expression data provided with the ctOWAS package.
"expr"



#' Example genotype matrix
#'
#' A small simulated genotype dosage matrix used in the README and examples.
#'
#' @format An integer matrix of dosages (0, 1, 2) with 20 samples in rows
#'   (\code{Sample1}, \ldots, \code{Sample20}) and 10 SNPs in columns
#'   (\code{SNP1}, \ldots, \code{SNP10}).
#' @source Simulated example data.
"x_example"




#' Example molecular trait vector
#'
#' A simulated molecular trait (e.g. expression of one gene) for the 20
#' samples in \code{\link{x_example}}, for use as \code{y} in
#' \code{\link{ctOWAS_train_K}}.
#'
#' @format A named numeric vector of length 20 (names \code{Sample1}, \ldots,
#'   \code{Sample20}).
#' @source Simulated example data.
"y_example"




#' Example GWAS summary statistics
#'
#' GWAS effect sizes and standard errors for the 10 SNPs in
#' \code{\link{x_example}}. Use \code{Beta / se_Beta} as
#' \code{gwas_z_score} in \code{\link{ctOWAS_assoc_test_K}}.
#'
#' @format A list with two named numeric vectors of length 10 (names
#'   \code{SNP1}, \ldots, \code{SNP10}):
#' \describe{
#'   \item{Beta}{Estimated SNP effect sizes.}
#'   \item{se_Beta}{Standard errors of \code{Beta}.}
#' }
#' @source Simulated example data.
"gwas_example"





#' Example cell-type fraction matrix (2 cell types)
#'
#' A simulated cell-type fraction matrix used in the README and examples,
#' for use as \code{pi_k} in \code{\link{ctOWAS_train_K}}.
#'
#' @format A numeric matrix with 20 samples in rows (\code{Sample1}, \ldots,
#'   \code{Sample20}) and 2 cell types in columns (\code{Cell1},
#'   \code{Cell2}). Rows sum to 1.
#'
#' @details
#' This dataset is artificially generated for demonstration purposes and is
#' not intended for biological interpretation.
#'
#' @source Simulated example data.
"pi_k"
