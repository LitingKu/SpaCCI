#########################################################
#' Load Example Data for SpaCCI
#'
#' The example data includes a Seurat object, a cell type proportion data frame, and a list of interest region Spot IDs.
#' These objects can be used to run example analyses with the SpaCCI package.
#'
#' @examples
#' \dontrun{
#'   # Load the example data
#'   load(system.file("extdata", "Tutorial_example_data.rda", package = "SpaCCI"))
#'   
#'   # Example usage
#'   print(Example_Seurat)
#'   head(cell_prop_df)
#'   print(interest_region_Spot_IDs)
#' }
#'
#' @export