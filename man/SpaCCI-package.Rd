\name{SpaCCI-package}
\alias{SpaCCI-package}
\docType{package}
\title{SpaCCI: Spatially Aware Cell-Cell Interaction Analysis}
\description{
  Provides tools for analyzing spatial cell-cell interactions based on ligand-receptor pairs.
  Functions include local, regional, and global analysis using spatial transcriptomics data.
}
\author{
  Li-Ting Ku
  Maintainer: Li-Ting Ku <lku@mdanderson.org>
}
\references{
  Some references that are relevant to the methods used in the package.
  - Author1 et al. (Year). Title of the paper. Journal Name, Volume(Issue), Pages.
}
\seealso{
  \code{\link{LR_database}}, \code{\link{run_SpaCCI}}
}
\examples{
\dontrun{
  # Example of how to use the main functions
  result <- LR_database(
    species = "Human",
    database_name = "CellChat",
    gene_spot_expression_dataframe = gene_spot_expression_dataframe
  )
}
}


