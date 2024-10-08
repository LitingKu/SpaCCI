% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visiualization.R
\name{plot_localized}
\alias{plot_localized}
\title{Plot Localized Hotspot Pattern}
\usage{
plot_localized(
  spatial_coord,
  resultdf_list,
  RegionIDs_matrix,
  celltype_ligand,
  celltype_receptor,
  plot_size,
  L_R_pair_name = NULL,
  alpha = 0.05
)
}
\arguments{
\item{spatial_coord}{A data frame of the spatial coordinates. The columns should include \code{"Spot_ID"}, \code{"imagerow"}, and \code{"imagecol"}. And the row names must be the names of \code{"Spot_ID"}, which is the same as the rownames in cell type proportion data frame or the colnames of the gene* spot expression data frame}

\item{resultdf_list}{A result of data frame list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`dataframelist`}}

\item{RegionIDs_matrix}{A result of matrix list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`RegionIDs_matrix`}}

\item{celltype_ligand}{Ligand cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.}

\item{celltype_receptor}{Receptor cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.}

\item{plot_size}{As this function incorporate with \code{Seurat}'s \code{`SpatialFeaturePlot`}, this parameter could control the plotting size of the each spot.}

\item{L_R_pair_name}{Initially this is set to \code{NULL}, if one is interested in a specific Ligand-Receptor pair, then one could specify the L_R_pair_name here. Note: the input name should match the L-R pair name exists in the dataframe in the output of SpaCCI_local "dataframelist".}

\item{alpha}{This is the significant cutoff for the adjusted-p-value of thr permutation test. Initially this is set to \code{0.05}, one could adjust the cutoff.}
}
\value{
The localized plot from the inferred cell-cell interaction on the local scale.
}
\description{
Visualize the inferred cell-cell interaction localized pattern if NOT using Seurat Object
}
\examples{
\donttest{
# Run localized hotspot plot
Result <- run_SpaCCI(..., analysis_scale = "local",...)
local_plot <- plot_localized(spatial_coord = spatial_coords_df,
                            resultdf_list = Result$dataframelist,
                            RegionIDs_matrix = Result$RegionIDs_matrix,
                            celltype_ligand = "Beta_cells",
                            celltype_receptor = "T_ells",
                            plot_size = 3)
}

}
