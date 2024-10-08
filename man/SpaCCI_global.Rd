% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling.R
\name{SpaCCI_global}
\alias{SpaCCI_global}
\title{Infer Cell-Cell Interactions on a Global Scale}
\usage{
SpaCCI_global(
  gene_spot_df,
  spot_cell_prop_df,
  matching_L_R_pairs,
  matching_L_R_pairs_info
)
}
\arguments{
\item{gene_spot_df}{A data frame where the rows are genes and the columns are spots (Spot_IDs), representing gene expression levels across spatial spots.}

\item{spot_cell_prop_df}{A data frame of cell type proportions for each spot. The rows represent spots (Spot_IDs), and the columns represent different cell types.}

\item{matching_L_R_pairs}{A data frame containing matching ligand-receptor pairs. Each row corresponds to a ligand-receptor pair, with columns for \code{ligand_vector} and \code{receptor_vector}.}

\item{matching_L_R_pairs_info}{A data frame providing additional information for each ligand-receptor pair, such as pathway information.}
}
\value{
A list containing:
\describe{
  \item{pvalue_df}{A data frame of inferred interactions across the global scale, including information on ligand and receptor cell types, interaction strength, P-values, and adjusted P-values.}
}
}
\description{
This function infers cell-cell interactions on a global scale using spatial transcriptomics data. It applies permutation testing to identify significant ligand-receptor interactions across all spots.
}
