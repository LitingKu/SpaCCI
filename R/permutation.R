###################################################################
#' Perform Global Permutations
#'
#' This function performs global permutations on the spatial transcriptomics data.
#'
#' @param permutationMatrix A matrix containing permutations.
#' @param permut_col A column matrix of permutations.
#' @param cellPropMatrix A matrix of cell type proportions.
#' @param spotGeneMatrix A matrix of gene expressions at spots.
#' @param LigandVectorIndex A vector of ligand indices.
#' @param ReceptorVectorIndex A vector of receptor indices.
#' @param null_expression A matrix of null expression values.
#' @param nBoot Number of bootstrap iterations.
#' @return A matrix with the results of the global permutations.
#' @export
#' @useDynLib SpaCCI, .registration = TRUE
#' @importFrom Rcpp sourceCpp
Global_Permutations <- function(permutationMatrix, permut_col, cellPropMatrix, spotGeneMatrix, LigandVectorIndex, ReceptorVectorIndex, null_expression, nBoot) {
  # This function is implemented in C++.
  # No R code is necessary here.
}

#' Perform Local and Regional Permutations
#'
#' This function performs local and regional permutations on the spatial transcriptomics data.
#'
#' @param permutationMatrix A matrix containing permutations.
#' @param permut_col A column matrix of permutations.
#' @param cellPropMatrix A matrix of cell type proportions.
#' @param spotGeneMatrix A matrix of gene expressions at spots.
#' @param LigandVectorIndex A vector of ligand indices.
#' @param ReceptorVectorIndex A vector of receptor indices.
#' @param null_expression A matrix of null expression values.
#' @param nBoot Number of bootstrap iterations.
#' @return A matrix with the results of the local and regional permutations.
#' @export
#' @useDynLib SpaCCI, .registration = TRUE
#' @importFrom Rcpp sourceCpp
Local_Regional_Permutations <- function(permutationMatrix, permut_col, cellPropMatrix, spotGeneMatrix, LigandVectorIndex, ReceptorVectorIndex, null_expression, nBoot) {
  # This function is implemented in C++.
  # No R code is necessary here.
}
