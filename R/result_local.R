#########################################################
#' result_local data for SpaCCI
#'
#' This example dataset is the result of running the \code{run_SpaCCI} function.
#' It contains the inferred cell-cell interactions across the global scale.
#'
#' These objects can be used for testing and running example analyses with the SpaCCI package.
#' @docType data
#' @usage data(result_local)
#'
#' @format
#' A list containing:
#' \describe{
#'   \item{dataframelist}{A list of data frame of the p-value reults of each spatial neighborhood.}
#'   \item{RegionIDs_matrix}{A list of matrix contains the spot IDs of each spatial neighborhood.}
#' }
#' @keywords datasets
#' @examples
#'   data(result_local)
#'   print(result_local)
"result_local"
