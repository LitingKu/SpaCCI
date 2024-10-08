% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visiualization.R
\name{plot_SpaCCI_heatmap}
\alias{plot_SpaCCI_heatmap}
\title{Plot SpaCCI Results on the heatmap
Visualize inferred significant cell-cell interactions using a heatmap}
\usage{
plot_SpaCCI_heatmap(
  SpaCCI_Result_List,
  specific_celltypes = NULL,
  pathways = NULL,
  interaction = NULL,
  log1p_transform = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  border_color = "white",
  fontsize_row = 11,
  fontsize_col = 11,
  family = "Arial",
  main = "",
  treeheight_col = 0,
  treeheight_row = 0,
  low_col = "dodgerblue4",
  mid_col = "peachpuff",
  high_col = "deeppink4",
  alpha = 0.05,
  return_tables = FALSE,
  symmetrical = FALSE,
  ...
)
}
\arguments{
\item{SpaCCI_Result_List}{A list containing the results from a SpaCCI \code{"regional"} or \code{"global"} analysis. This list should include \code{pvalue_df}, which are the outputs from \code{run_SpaCCI(..., analysis_scale = "regional",...)} or  \code{run_SpaCCI(..., analysis_scale = "global",...)}.}

\item{specific_celltypes}{A vector of cell types to include in the heatmap, i.e c("Celltype_A","Celltype_B"). NOTE: the cell type names should match the names input in the SpaCCI analysis.}

\item{pathways}{A vector of pathways to filter the interactions. Initially set to \code{NULL}, if not, then it will aggregate the results of the selected pathways.}

\item{interaction}{A vector of interactions to filter. Initially set to \code{NULL}, if not, then it will aggregate the results of the selected interactions.}

\item{log1p_transform}{Logical; whether to apply a log(1 + x) transformation to the count matrix.}

\item{show_rownames}{Logical; whether to show row names in the heatmap.}

\item{show_colnames}{Logical; whether to show column names in the heatmap.}

\item{scale}{Character; whether to scale the data ("row", "column", "none").}

\item{cluster_cols}{Logical; whether to cluster columns.}

\item{cluster_rows}{Logical; whether to cluster rows.}

\item{border_color}{Character; color of the heatmap borders.}

\item{fontsize_row}{Numeric; font size for row names.}

\item{fontsize_col}{Numeric; font size for column names.}

\item{family}{Character; font family for text in the heatmap.}

\item{main}{Character; title of the heatmap.}

\item{treeheight_col}{Numeric; height of the column dendrogram.}

\item{treeheight_row}{Numeric; height of the row dendrogram.}

\item{low_col}{Character; color for low values in the heatmap.}

\item{mid_col}{Character; color for mid values in the heatmap.}

\item{high_col}{Character; color for high values in the heatmap.}

\item{alpha}{Numeric; significance threshold for p-values, initailly set to \code{0.05}.}

\item{return_tables}{Logical; whether to return the count matrix and summary tables.}

\item{symmetrical}{Logical; whether to make the heatmap symmetrical.}

\item{...}{Additional arguments passed to `pheatmap`.}
}
\value{
If `return_tables` is \code{FALSE} (default), the function returns a heatmap object created by \code{pheatmap}, showing the count of significant cell-cell interactions. If `return_tables` is \code{TRUE}, the function returns a list containing:
\describe{
  \item{\strong{heatmap}}{The heatmap object showing the significant cell-cell interactions.}
  \item{\strong{heatmap_countmatrix}}{The matrix used to generate the heatmap, with cell types as rows and columns, and counts of significant interactions as values.}
  \item{\strong{table}}{A data frame summarizing the counts of significant interactions between each ligand and receptor cell type combination.}
}
}
\description{
Plot SpaCCI Results on the heatmap
Visualize inferred significant cell-cell interactions using a heatmap
}
\examples{
library(SpaCCI)
library(dplyr)
library(reshape2)
library(grDevices)
library(pheatmap)
data(result_global)
celltypes <- c("beta" , "delta" , "ductal","macrophage",
                "activated_stellate", "quiescent_stellate")
plot_SpaCCI_heatmap(SpaCCI_Result_List = result_global,
                    symmetrical = FALSE, cluster_cols = FALSE, return_tables = FALSE,
                    cluster_rows = FALSE, #cellheight = 10, cellwidth = 10,
                    specific_celltypes = c(celltypes),
                    main= "Cell-Cell Interaction Count")

}
