#' Local Spatial Coordinates Data Frame
#'
#' A dataset containing spatial coordinates for local analysis, typically used in spatial transcriptomics.
#'
#' @docType data
#' @usage data(result_local)
#' 
#' @format A data frame with X rows and 3 columns:
#' \describe{
#'   \item{imagecol}{The x-coordinate of the spot in the image.}
#'   \item{imagerow}{The y-coordinate of the spot in the image.}
#'   \item{Spot_ID}{The unique identifier for each spatial spot.}
#' }
#'
#' @keywords datasets
#' @examples
#' data(result_local_spatial_coords_df)
#' head(result_local_spatial_coords_df)
"result_local_spatial_coords_df"

