#########################################################
#' result_local data and spatial coordinates for SpaCCI
#' Local Regional Data
#'
#' This dataset contains the results of spatial neighborhood analyses.
#' These objects can be used for testing and running example analyses with the SpaCCI package.
#' @docType data
#' @usage data(result_local)
#'
#' @format
#' A list containing:
#' \describe{
#'   \item{dataframelist}{A list of data frames, each containing p-value results for spatial neighborhoods.}
#'   \item{RegionIDs_matrix}{A list of matrices, each containing the spot IDs for specific neighborhoods.}
#' }
#' @keywords datasets
#' @examples
#'   data(result_local)
#'   print(result_local)
"result_local"
