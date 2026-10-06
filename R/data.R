#' Cell-type marker gene reference
#'
#' Marker-gene sets for 489 cell-type signatures, compiled from Aran et al.
#' (2017), the reference used by xCell.
#'
#' @format A data frame with 489 rows and 202 columns:
#' \describe{
#'   \item{Celltype_Source_ID}{Signature identifier (cell type, source and
#'     replicate, e.g. \code{"aDC_HPCA_1"}).}
#'   \item{# of genes}{Number of marker genes in the signature.}
#'   \item{V3, \ldots, V202}{Marker gene symbols; unused cells are empty
#'     strings.}
#' }
#' @source Aran D, Hu Z, Butte AJ (2017). xCell: digitally portraying the
#'   tissue cellular heterogeneity landscape. \emph{Genome Biology} 18, 220.
"GeneXCell"

#' Example bulk expression matrix
#'
#' A small simulated bulk expression matrix used in the README and examples.
#'
#' @format A numeric matrix with 50 genes in rows (\code{Gene1}, \ldots,
#'   \code{Gene50}) and 20 samples in columns (\code{Sample1}, \ldots,
#'   \code{Sample20}).
#' @source Simulated example data.
"exprB_example"

#' Example marker gene list
#'
#' Marker genes for two cell types, for use as \code{marker_list} in
#' \code{\link{pi_estimation_K}}.
#'
#' @format A named list of two character vectors, \code{CellType1}
#'   (\code{Gene1}--\code{Gene3}) and \code{CellType2}
#'   (\code{Gene4}--\code{Gene6}). Genes are row names of
#'   \code{\link{exprB_example}}.
#' @source Simulated example data.
"markers_example"

#' Example simulated xCell enrichment scores
#'
#' Simulated xCell-like enrichment scores for the 20 samples in
#' \code{\link{exprB_example}}, for use as \code{xCellScore} in
#' \code{\link{estimate_prior}}.
#'
#' @format A numeric matrix with 5 xCell cell types in rows
#'   (\code{"Epithelial cells"}, \code{"Adipocytes"}, \code{"Preadipocytes"},
#'   \code{"Fibroblasts"}, \code{"Endothelial cells"}) and 20 samples in
#'   columns (\code{Sample1}, \ldots, \code{Sample20}). Values are
#'   nonnegative.
#'
#' @details
#' Scores were generated as a scaled copy of the true fractions in
#' \code{\link{pi_k}} plus Gaussian noise, truncated at zero:
#' \code{"Epithelial cells"} tracks \code{pi_k[, "Cell1"]}, and
#' \code{"Adipocytes"} and \code{"Preadipocytes"} track
#' \code{pi_k[, "Cell2"]}. \code{"Fibroblasts"} and
#' \code{"Endothelial cells"} are unrelated noise. The generating code is in
#' \code{data-raw/example_data.R} in the source repository.
#'
#' @source Simulated example data.
#'
#' @examples
#' data(xCellScore_example)
#' prior <- estimate_prior(
#'   xCellScore = xCellScore_example,
#'   cell_types = list(CellType1 = "Epithelial cells",
#'                     CellType2 = c("Adipocytes", "Preadipocytes")),
#'   mu = c(0.5, 0.5), lower = 0.01, upper = 0.99
#' )
#' head(prior)
"xCellScore_example"

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

#' Example merged results table
#'
#' A toy results table with marginal p-values for each cell type, a joint
#' p-value and a cell-type specificity label.
#'
#' @format A data frame with 300 rows (genes) and 5 columns:
#' \describe{
#'   \item{gene}{Gene identifier.}
#'   \item{type_ct2}{Cell-type specificity label (e.g.
#'     \code{"CellTypeSpecific"}).}
#'   \item{p_1_ct2, p_2_ct2}{Marginal p-values for cell types 1 and 2.}
#'   \item{p_join_ct2}{Joint p-value across the two cell types.}
#' }
#' @source Simulated example data.
"merged_example"

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
