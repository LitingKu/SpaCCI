###############################################################################################
# modeling
#' Find the spatial neighborhood Spot IDs
#'
#' This function identifies the spatial neighborhood Spot IDs around a given center Spot ID within a specified radius.
#'
#' @param object An object that could be either 1. A Seurat object, or 2. A data frame where the columns are Spot_IDs (i.e., the gene*spot expression matrix).
#' @param spatial_coord A data frame of the spatial coordinates. The column names should include `c("Spot_ID", "imagerow", "imagecol")`, and the row names must be the Spot_ID, which is the same as the row names in the cell type proportion data frame or the column names of the gene*spot expression data frame.
#' @param centerID A vector of length 1, representing a single Spot_ID that serves as the center for the neighborhood.
#' @param radius The radius of the spatial neighborhood, specified as a numeric value.
#' @param enhanced Logical; if `TRUE`, enhances the Seurat object by marking the neighborhood in a special way. Defaults to `FALSE`.
#' @param avern Numeric; the number of samples to average over when determining the unit distance. Defaults to 5.
#'
#' @return A list containing:
#' \describe{
#'   \item{centerID}{The input center Spot_ID.}
#'   \item{closeID}{The Spot_IDs of the neighboring spots within the specified radius.}
#'   \item{unit}{The unit distance used to determine the neighborhood.}
#' }
#'
#' @export
#'
Find_regional_IDs <- function(object,
                              spatial_coord,
                              centerID,
                              enhanced = FALSE,
                              radius  ,
                              avern = 5) {

  stopifnot(exprs = {
    is.character(centerID)
    is.numeric(radius)
    is.numeric(avern)
  })

  #xloc <- object@images$slice1@coordinates$imagecol
  #yloc <- object@images$slice1@coordinates$imagerow
  xloc <- spatial_coord$imagecol
  yloc <- spatial_coord$imagerow

  ID <- sample(seq_len(length(xloc)), size = avern)

  rec_minDist <- rep(NA, avern)

  for(i in seq_len(avern)) {
    rec_minDist[i] <- min(sqrt((xloc[ID[i]] - xloc[-ID[i]])^2 + (yloc[ID[i]] - yloc[-ID[i]])^2))
  }
  unitDist <- mean(rec_minDist)


  if(enhanced) {
    object.enhanced <- object
    allspot <- colnames(object)
    center_x <- xloc[which(allspot == centerID)]
    center_y <- yloc[which(allspot == centerID)]

    alldist <- sqrt((xloc - center_x)^2 + (yloc - center_y)^2)
    closeID <- allspot[which(alldist <= (radius * unitDist+0.1))]

    featureDot <- rep(3, length(allspot))
    featureDot[allspot %in% closeID] <- 2
    featureDot[allspot == centerID] <- (-2)
  } else {
    center_x <- xloc[which(colnames(object) == centerID)]
    center_y <- yloc[which(colnames(object) == centerID)]

    alldist <- sqrt((xloc - center_x)^2 + (yloc - center_y)^2)
    closeID <- colnames(object)[which(alldist <= (radius * unitDist+0.1))]

  }

  return(list(centerID = centerID,
              closeID = closeID,
              unit= unitDist))

}

#' Perform a Deranged Shuffle of Cell Types
#'
#' This function takes a vector of cell types and returns a shuffled version where no element remains in its original position.
#'
#' @param CellType A character vector representing the cell types to be shuffled.
#'
#' @return A character vector of the same length as `CellType`, with elements shuffled such that no element remains in its original position.
#'
#' @examples
#' \donttest{
#' original <- c("B_cell", "T_cell", "NK_cell", "Macrophage")
#' shuffled <- GetShuffledCT(original)
#' print(shuffled)
#' }
#'
#' @export'
#'
GetShuffledCT <- function(CellType) {
  tmprec <- sample(CellType, replace = FALSE)
  idx = (tmprec == CellType)
  while(any(idx)) {
    tmprec[idx] <- sample(tmprec[idx], sum(idx), replace = FALSE)
    idx = (tmprec == CellType)
  }
  return(tmprec)
}


#' Select Closest Spatial IDs to a Center Point: this is used for permutation
#'
#' This function identifies and returns the IDs of the closest spatial points to a specified center point based on Euclidean distance.
#'
#' @param spatial_coord A data frame of the spatial coordinates. The column names should include `c("Spot_ID", "imagerow", "imagecol")`, and the row names must be the Spot_IDs, which is the same as the row names in the cell type proportion data frame or the column names of the gene*spot expression data frame.
#' @param center_id A character string specifying the ID of the center spot from which distances are calculated.
#' @param n_ids An integer specifying the number of closest IDs to select.
#'
#' @return A character vector of the `n_ids` closest IDs to the specified center ID.
#'
#' @examples
#' \donttest{
#' spatial_coord <- data.frame(
#'   imagecol = c(1, 2, 3, 4, 5),
#'   imagerow = c(5, 4, 3, 2, 1),
#'   row.names = c("Spot1", "Spot2", "Spot3", "Spot4", "Spot5")
#' )
#' center_id <- "Spot3"
#' closest_ids <- random_region(spatial_coord, center_id, 3)
#' print(closest_ids)
#'}
#' @importFrom dplyr arrange
#'
#' @export

random_region <- function(spatial_coord, center_id, n_ids){
  center_coords <- spatial_coord[center_id, ]
  distances <- sqrt((spatial_coord$imagecol - center_coords$imagecol)^2 + (spatial_coord$imagerow - center_coords$imagerow)^2)
  distance_df <- data.frame(ID = rownames(spatial_coord), Distance = distances)
  sorted_distance_df <- distance_df %>% arrange(Distance)
  # Select the closest IDs based on the specified number
  selected_ids <- head(sorted_distance_df, n_ids)$ID
  return(selected_ids)
}


#' Infer Cell-Cell Interactions on a Local Scale
#'
#' This function infers cell-cell interactions on a local scale using spatial transcriptomics data. It utilizes permutation testing to identify significant ligand-receptor interactions within specified neighborhoods around randomly selected center spots.
#'
#' @param gene_spot_df A data frame where the rows are genes and the columns are spots (Spot_IDs), representing gene expression levels across spatial spots.
#' @param spot_cell_prop_df A data frame of cell type proportions for each spot. The rows represent spots (Spot_IDs), and the columns represent different cell types.
#' @param spatial_coord A data frame of the spatial coordinates. The column names should include `c("Spot_ID", "imagerow", "imagecol")`, and the row names must be the Spot_IDs, which is the same as the row names in the cell type proportion data frame or the column names of the gene*spot expression data frame.
#' @param prop A numeric value representing the proportion of spots to randomly sample as center spots for local neighborhood analysis.
#' @param radius A numeric value specifying the radius of the spatial neighborhood around each center spot.
#' @param matching_L_R_pairs A data frame containing matching ligand-receptor pairs. Each row corresponds to a ligand-receptor pair, with columns for \code{ligand_vector} and \code{receptor_vector}.
#' @param matching_L_R_pairs_info A data frame providing additional information for each ligand-receptor pair, such as pathway information.
#'
#' @return A list containing:
#' \describe{
#'   \item{dataframelist}{A list of data frames, each representing the inferred interactions for a specific center spot. Each data frame includes information on ligand and receptor cell types, P-values, and adjusted P-values.}
#'   \item{RegionIDs_matrix}{A list of matrices, each containing the IDs of the spots within the specified radius of each center spot.}
#' }
#'
#' @importFrom nnls nnls
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
SpaCCI_local <- function(gene_spot_df,
                         spot_cell_prop_df,
                         spatial_coord,
                         prop,
                         radius,
                         matching_L_R_pairs,
                         matching_L_R_pairs_info){

  dataframelist <- list()
  centerIDs <- sample(c(colnames(gene_spot_df)), size = round(ncol(gene_spot_df)*prop) )
  #centerIDs <- c(colnames(gene_spot_df))
  RegionIDs_matrix <- list()
  message(paste0("Now finding the neighborhoods for ",round(ncol(gene_spot_df)*prop), " spots"))
  for ( i in 1:length(centerIDs)){
    IDs <- Find_regional_IDs(gene_spot_df,spatial_coord, centerIDs[i], enhanced = FALSE, radius = radius, avern = 5)$closeID
    RegionIDs_matrix[[centerIDs[i]]] <- IDs
  }
  RegionIDs_matrix <- RegionIDs_matrix[sapply(RegionIDs_matrix, length) >= 6] # remove those that don't have enough neighborhoods
  centerIDs <- names(RegionIDs_matrix) # update cencter ID
  
  ## prepared for permutation
  cellPropMatrix <- as.matrix(spot_cell_prop_df)
  GeneSpotMatrix <- as.matrix(gene_spot_df)
  pb <- txtProgressBar(min = 0, max =  length(centerIDs), style = 3, file = stderr())

  for (id in 1: length(centerIDs)){
    setTxtProgressBar(pb, id)
    ids <- c(RegionIDs_matrix[[ centerIDs[id] ]])
    ids <- ids[ids %in% colnames(gene_spot_df)]
    nboot <- as.numeric(floor(nrow(spot_cell_prop_df)/500)*100)
    if (nboot == 0){
      nboot <- 1
    }else if(nboot > 1000){
      nboot <- 1000
    }

    premut_center <- sample(colnames(gene_spot_df[,-which(colnames(gene_spot_df) %in% ids )]), size = nboot)
    region_permut_index <- vector("list", length(premut_center))

    # Loop through each center ID
    for (i in seq_along(premut_center)) {
      region_ids <- random_region(spatial_coord, premut_center[i], n_ids = length(ids))
      region_permut_index[[i]] <- as.numeric(match(region_ids, colnames(gene_spot_df) ))
      if(any(is.na( region_permut_index[[i]]))){
        stop("Please make sure the rownames of spatial_coordinates_dataframe match the colnames of gene_spot_expression_dataframe")
      }
      
    }

    # Optionally, convert list to a data frame or matrix if all outputs have the same length
    permutation <- do.call(cbind, lapply(region_permut_index, function(x) as.character(x)))
    permut_col <- replicate(nboot, GetShuffledCT(length(colnames(spot_cell_prop_df))))
    permutationMatrix <- matrix(as.numeric(permutation), nrow = nrow(permutation))


    matching_indices <- numeric(length(ids)) # construct for null indices
    for (ss in seq_along(ids)) {
      matching_indices[ss] <- which(colnames(gene_spot_df) == ids[ss])
    }
    permut_null_regionMatrix <- as.matrix(matching_indices)
    output <- list()

    ############ start ############
    for (i in 1:nrow(matching_L_R_pairs)){
      avg_ligand_list <- vector("list", length = length(matching_L_R_pairs$ligand_vector[[i]]))
      avg_receptor_list <- vector("list", length = length(matching_L_R_pairs$receptor_vector[[i]]))
      #.   Ligand
      for (l in 1: length(c(matching_L_R_pairs$ligand_vector[[i]])) ){

        y1 <- t(gene_spot_df[c(matching_L_R_pairs$ligand_vector[[i]][l]),ids])
        x1 <- cellPropMatrix[ids, ]

        if (length(ids) < 2){
          truth <- rep(0, length(colnames(spot_cell_prop_df)))
          names(truth) <- colnames(spot_cell_prop_df)
        }else{
          mod1 <- nnls(x1, y1)
          truth <- mod1$x
          names(truth) <- colnames(spot_cell_prop_df)

          threshols1 <- mean(y1) / length(colnames(spot_cell_prop_df))
          threshols1.5 <- mean(y1)*(colSums(spot_cell_prop_df)/nrow(spot_cell_prop_df))
          threshols2 <-  1/length(truth > 0 )
          truth[(colMeans(sweep(x1, 2, truth, `*`) ) < threshols1  & colMeans(sweep(x1, 2, truth, `*`) ) < threshols1.5 ) | truth/sum(truth) < threshols2] <- 0


        }


        avg_ligand_list[[l]] <- truth


      }
      # Calculate the product of corresponding elements across all vectors
      product_vector <- Reduce("*", avg_ligand_list)
      # Calculate the nth root of the product, where n is the number of vectors
      n <- length(avg_ligand_list)
      null_avg_ligand <- product_vector^(1/n)



      #.   Receptor
      for (r in 1: length(c(matching_L_R_pairs$receptor_vector[[i]])) ){


        ###True expression from the region ########
        y1 <- t(gene_spot_df[c(matching_L_R_pairs$receptor_vector[[i]][r]),ids])
        x1 <- cellPropMatrix[ids, ]

        if (length(ids) < 2){
          truth <- rep(0, length(colnames(spot_cell_prop_df)))
          names(truth) <- colnames(spot_cell_prop_df)
        }else{
          mod1 <- nnls(x1, y1)
          truth <- mod1$x
          names(truth) <- colnames(spot_cell_prop_df)

          threshols1 <- mean(y1) / length(colnames(spot_cell_prop_df))
          threshols1.5 <- mean(y1)*(colSums(spot_cell_prop_df)/nrow(spot_cell_prop_df))
          threshols2 <-  1/length(truth > 0 )
          truth[ (colMeans(sweep(x1, 2, truth, `*`) ) < threshols1  & colMeans(sweep(x1, 2, truth, `*`) ) < threshols1.5 ) | truth/sum(truth) < threshols2] <- 0

        }


        avg_receptor_list[[r]] <- truth

      }
      # Calculate the product of corresponding elements across all vectors
      product_vector <- Reduce("*", avg_receptor_list)
      # Calculate the nth root of the product, where n is the number of vectors
      n <- length(avg_receptor_list)
      null_avg_receptor <- product_vector^(1/n)
      t <- outer(null_avg_ligand, null_avg_receptor, `*`)


      prop_t <- outer(colSums(spot_cell_prop_df[ids, ])/nrow(spot_cell_prop_df[ids, ]), colSums(spot_cell_prop_df[ids, ])/nrow(spot_cell_prop_df[ids, ]), `*`)
      transformed_t <- t / (0.5 + t)
      null_expression <- as.matrix(transformed_t * prop_t )


      ################# Permutation ##############
      LigandVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$ligand_vector[[i]])))
      ReceptorVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$receptor_vector[[i]])))
      # Call the function

      averageResult <- .Call("_SpaCCI_Local_Regional_Permutations",
                             permutationMatrix,
                             permut_col,
                             cellPropMatrix,
                             GeneSpotMatrix,
                             LigandVectorIndex,
                             ReceptorVectorIndex,
                             null_expression,
                             nboot)

      #averageResult <- Local_Regional_Permutations(permutationMatrix,
      #                                             permut_col,
      #                                             cellPropMatrix,
      #                                             GeneSpotMatrix,
      #                                             LigandVectorIndex,
      #                                             ReceptorVectorIndex,
      #                                             null_expression,
      #                                             nboot)



      colnames(averageResult) <- colnames(spot_cell_prop_df)
      rownames(averageResult) <- colnames(spot_cell_prop_df)

      non_zero_indices <- which(null_expression != 0,arr.ind = TRUE)
      nul <- data.frame(
        Cell_type_Ligand = names(null_avg_ligand)[non_zero_indices[,1]],
        Cell_type_Receptor = names(null_avg_receptor)[non_zero_indices[,2]]
      )
      re <- averageResult[non_zero_indices]
      nul$PValue <- re
      nul$adjusted.PValue <- p.adjust(nul$PValue, method = "BH")
      nul$adjusted.PValue <- nul$PValue

      output[[i]] <- nul

    }

    table <- lapply(seq_along(output), function(i) {
      #df <- get_elements(output[[i]])
      df <- output[[i]]
      interaction_info <- matching_L_R_pairs_info[i, ]
      #interaction_info <- cbind(interaction_info, avgL_mat[i,], avgR_mat[i, ] )
      interaction_expanded <- interaction_info[rep(1, nrow(df)), ]
      rownames(interaction_expanded) <- NULL
      df <- cbind(df, interaction_expanded)
      return(df)
    })

    dataframe <- do.call(rbind, table)


    dataframelist[[centerIDs[id]]] <- dataframe




  }

  close(con = pb)
  message("writing data frame")
  return(list(dataframelist = dataframelist,RegionIDs_matrix = RegionIDs_matrix))



}


#' Infer Cell-Cell Interactions in a Specified Region
#'
#' This function infers cell-cell interactions within a specified region using spatial transcriptomics data. It applies permutation testing to identify significant ligand-receptor interactions in the region.
#'
#' @param gene_spot_df A data frame where the rows are genes and the columns are spots (Spot_IDs), representing gene expression levels across spatial spots.
#' @param spot_cell_prop_df A data frame of cell type proportions for each spot. The rows represent spots (Spot_IDs), and the columns represent different cell types.
#' @param region_spot_IDs A vector of Spot_IDs representing the spots included in the region of interest.
#' @param matching_L_R_pairs A data frame containing matching ligand-receptor pairs. Each row corresponds to a ligand-receptor pair, with columns for \code{ligand_vector} and \code{receptor_vector}.
#' @param matching_L_R_pairs_info A data frame providing additional information for each ligand-receptor pair, such as pathway information.
#'
#' @return A list containing:
#' \describe{
#'   \item{pvalue_df}{A data frame of inferred interactions within the specified region, including information on ligand and receptor cell types, P-values, and adjusted P-values.}
#' }
#'
#' @importFrom nnls nnls
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'

SpaCCI_region <- function(gene_spot_df,
                          spot_cell_prop_df,
                          region_spot_IDs,
                          matching_L_R_pairs,
                          matching_L_R_pairs_info){
  ## prepared for permutation
  cellPropMatrix <- as.matrix(spot_cell_prop_df)
  GeneSpotMatrix <- as.matrix(gene_spot_df)
  # construct permutation IDs
  nboot <- 1000
  indices <- as.numeric(match(colnames(gene_spot_df[,-which(colnames(gene_spot_df) %in% region_spot_IDs )]), colnames(gene_spot_df) ))
  permutation <- replicate(nboot, sample(indices, size = length(region_spot_IDs)))
  permut_col <- replicate(nboot, GetShuffledCT(length(colnames(spot_cell_prop_df))) )


  permutationMatrix <- as.matrix(permutation)
  matching_indices <- numeric(length(region_spot_IDs))
  for (ss in seq_along(region_spot_IDs)) {
    matching_indices[ss] <- which(colnames(gene_spot_df) == region_spot_IDs[ss])
  }
  permut_null_regionMatrix <- as.matrix(matching_indices)

  output <- list()


  pb <- txtProgressBar(min = 0, max =  nrow(matching_L_R_pairs), style = 3, file = stderr())
  ############ start ############
  for (i in 1:nrow(matching_L_R_pairs)){

    setTxtProgressBar(pb, i)
    avg_ligand_list <- vector("list", length = length(matching_L_R_pairs$ligand_vector[[i]]))
    avg_ligand_list_f <- vector("list", length = length(matching_L_R_pairs$ligand_vector[[i]]))
    avg_receptor_list <- vector("list", length = length(matching_L_R_pairs$receptor_vector[[i]]))
    avg_receptor_list_f <- vector("list", length = length(matching_L_R_pairs$receptor_vector[[i]]))

    #.   Ligand
    for (l in 1: length(c(matching_L_R_pairs$ligand_vector[[i]])) ){

      y1 <- t(gene_spot_df[c(matching_L_R_pairs$ligand_vector[[i]][l]),region_spot_IDs])
      x1 <- cellPropMatrix[region_spot_IDs, ]

      if (length(region_spot_IDs) < 2){
        truth <- rep(0, length(colnames(spot_cell_prop_df)))
        names(truth) <- colnames(spot_cell_prop_df)
      }else{
        mod1 <- nnls(x1, y1)
        truth <- mod1$x
        names(truth) <- colnames(spot_cell_prop_df)

        truth_f <-truth
        threshols1 <- mean(y1) / length(colnames(spot_cell_prop_df))
        threshols2 <-  1/sum(truth > 0 )
        truth_f[colMeans(sweep(x1, 2, truth, `*`) ) < threshols1| truth/sum(truth) < threshols2] <- 0


      }


      avg_ligand_list[[l]] <- truth
      avg_ligand_list_f[[l]] <- truth_f

    }
    # Calculate the product of corresponding elements across all vectors
    product_vector <- Reduce("*", avg_ligand_list)
    product_vector_f <- Reduce("*", avg_ligand_list_f)
    # Calculate the nth root of the product, where n is the number of vectors
    n <- length(avg_ligand_list)
    null_avg_ligand <- product_vector^(1/n)
    null_avg_ligand_f <- product_vector_f^(1/n)


    #.   Receptor
    for (r in 1: length(c(matching_L_R_pairs$receptor_vector[[i]])) ){


      ###True expression from the region ########
      y1 <- t(gene_spot_df[c(matching_L_R_pairs$receptor_vector[[i]][r]),region_spot_IDs])
      x1 <- cellPropMatrix[region_spot_IDs, ]

      if (length(region_spot_IDs) < 2){
        truth <- rep(0, length(colnames(spot_cell_prop_df)))
        names(truth) <- colnames(spot_cell_prop_df)
      }else{
        mod1 <- nnls(x1, y1)
        truth <- mod1$x
        names(truth) <- colnames(spot_cell_prop_df)
        truth_f <- truth
        threshols1 <- mean(y1) / length(colnames(spot_cell_prop_df))
        threshols2 <-  1/sum(truth > 0 )
        truth_f[colMeans(sweep(x1, 2, truth, `*`) ) < threshols1 | truth/sum(truth) < threshols2] <- 0

      }


      avg_receptor_list[[r]] <- truth
      avg_receptor_list_f[[r]] <- truth_f

    }
    # Calculate the product of corresponding elements across all vectors
    product_vector <- Reduce("*", avg_receptor_list)
    product_vector_f <- Reduce("*", avg_receptor_list_f)
    # Calculate the nth root of the product, where n is the number of vectors
    n <- length(avg_receptor_list)
    null_avg_receptor <- product_vector^(1/n)
    null_avg_receptor_f <- product_vector_f^(1/n)
    t <- outer(null_avg_ligand, null_avg_receptor, `*`)
    t_f <- outer(null_avg_ligand_f, null_avg_receptor_f, `*`)

    prop <- spot_cell_prop_df[region_spot_IDs, ]
    prop_t <- outer(colSums(prop)/nrow(spot_cell_prop_df[region_spot_IDs, ]), colSums(prop)/nrow(spot_cell_prop_df[region_spot_IDs, ]), `*`)
    transformed_t <- t / (0.5 + t)
    transformed_t_f <- t_f / (0.5 + t_f)
    null_expression <- as.matrix(transformed_t * prop_t )
    null_expression_f <- as.matrix(transformed_t_f * prop_t )

    ################# Permutation ##############
    LigandVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$ligand_vector[[i]])))
    ReceptorVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$receptor_vector[[i]])))
    # Call the function

    averageResult <- .Call("_SpaCCI_Local_Regional_Permutations",
                           permutationMatrix,
                           permut_col,
                           cellPropMatrix,
                           GeneSpotMatrix,
                           LigandVectorIndex,
                           ReceptorVectorIndex,
                           null_expression,
                           nboot)


    #averageResult <- Local_Regional_Permutations(permutationMatrix,
    #                                             permut_col,
    #                                             cellPropMatrix,
    #                                             GeneSpotMatrix,
    #                                             LigandVectorIndex,
    #                                             ReceptorVectorIndex,
    #                                             null_expression,
    #                                             nboot)



    colnames(averageResult) <- colnames(spot_cell_prop_df)
    rownames(averageResult) <- colnames(spot_cell_prop_df)
    #p_vector <- as.vector(averageResult)

    t <- outer(null_avg_ligand, null_avg_receptor, `*`) # row are ligand, col are receptor
    non_zero_indices <- which(null_expression_f != 0,arr.ind = TRUE)
    nul <- data.frame(
      Cell_type_Ligand = names(null_avg_ligand)[non_zero_indices[,1]],
      Cell_type_Receptor = names(null_avg_receptor)[non_zero_indices[,2]]
    )
    re <- averageResult[non_zero_indices]
    nul$PValue <- re
    nul$adjusted.PValue <- p.adjust(nul$PValue, method = "BH")
    nul$adjusted.PValue <- nul$PValue

    output[[i]] <- nul

  }


  close(con = pb)

  message("writing data frame")
  ##############

  table <- lapply(seq_along(output), function(i) {
    df <- output[[i]]
    interaction_info <- matching_L_R_pairs_info[i, ,drop=FALSE]
    interaction_expanded <- interaction_info[rep(1, nrow(df)), ,drop=FALSE]
    rownames(interaction_expanded) <- NULL
    df <- cbind(df, interaction_expanded)
    return(df)
  })
  dataframe <- do.call(rbind, table)




  return(list(pvalue_df = dataframe ))
}


#' Infer Cell-Cell Interactions on a Global Scale
#'
#' This function infers cell-cell interactions on a global scale using spatial transcriptomics data. It applies permutation testing to identify significant ligand-receptor interactions across all spots.
#'
#' @param gene_spot_df A data frame where the rows are genes and the columns are spots (Spot_IDs), representing gene expression levels across spatial spots.
#' @param spot_cell_prop_df A data frame of cell type proportions for each spot. The rows represent spots (Spot_IDs), and the columns represent different cell types.
#' @param matching_L_R_pairs A data frame containing matching ligand-receptor pairs. Each row corresponds to a ligand-receptor pair, with columns for \code{ligand_vector} and \code{receptor_vector}.
#' @param matching_L_R_pairs_info A data frame providing additional information for each ligand-receptor pair, such as pathway information.
#'
#' @return A list containing:
#' \describe{
#'   \item{pvalue_df}{A data frame of inferred interactions across the global scale, including information on ligand and receptor cell types, interaction strength, P-values, and adjusted P-values.}
#' }
#'
#' @importFrom nnls nnls
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
SpaCCI_global <- function(gene_spot_df,
                          spot_cell_prop_df,
                          matching_L_R_pairs,
                          matching_L_R_pairs_info){
  ## prepared for permutation
  cellPropMatrix <- as.matrix(spot_cell_prop_df)
  GeneSpotMatrix <- as.matrix(gene_spot_df)
  region_spot_IDs <- rownames(spot_cell_prop_df)
  # construct permutation IDs
  nboot <- 500
  permutation <- replicate(nboot, sample(ncol(gene_spot_df), size = length(colnames(gene_spot_df))))
  permut_col <- replicate(nboot, GetShuffledCT(length(colnames(spot_cell_prop_df))) )


  permutationMatrix <- as.matrix(permutation)
  matching_indices <- numeric(length(region_spot_IDs))
  for (ss in seq_along(region_spot_IDs)) {
    matching_indices[ss] <- which(colnames(gene_spot_df) == region_spot_IDs[ss])
  }
  permut_null_regionMatrix <- as.matrix(matching_indices)
  output <- list()


  pb <- txtProgressBar(min = 0, max =  nrow(matching_L_R_pairs), style = 3, file = stderr())
  ############ start ############
  for (i in 1:nrow(matching_L_R_pairs)){
    setTxtProgressBar(pb, i)
    avg_ligand_list <- vector("list", length = length(matching_L_R_pairs$ligand_vector[[i]]))
    avg_receptor_list <- vector("list", length = length(matching_L_R_pairs$receptor_vector[[i]]))
    avg_ligand_list_f <- vector("list", length = length(matching_L_R_pairs$ligand_vector[[i]]))
    avg_receptor_list_f <- vector("list", length = length(matching_L_R_pairs$receptor_vector[[i]]))
    #.   Ligand
    for (l in 1: length(c(matching_L_R_pairs$ligand_vector[[i]])) ){

      y1 <- t(gene_spot_df[c(matching_L_R_pairs$ligand_vector[[i]][l]),region_spot_IDs])
      x1 <- cellPropMatrix[region_spot_IDs, ]

      if (length(region_spot_IDs) < 2){
        truth <- rep(0, length(colnames(spot_cell_prop_df)))
        names(truth) <- colnames(spot_cell_prop_df)
      }else{
        mod1 <- nnls(x1, y1)
        truth <- mod1$x
        names(truth) <- colnames(spot_cell_prop_df)
        truth_f <- truth
        threshols2 <-  1/sum(truth > 0 )
        truth_f[ truth/sum(truth) < threshols2] <- 0


      }

      avg_ligand_list[[l]] <- truth
      avg_ligand_list_f[[l]] <- truth_f

    }
    # Calculate the product of corresponding elements across all vectors
    product_vector <- Reduce("*", avg_ligand_list)
    product_vector_f <- Reduce("*", avg_ligand_list_f)
    # Calculate the nth root of the product, where n is the number of vectors
    n <- length(avg_ligand_list)
    null_avg_ligand <- product_vector^(1/n)
    null_avg_ligand_f <- product_vector_f^(1/n)


    #.   Receptor
    for (r in 1: length(c(matching_L_R_pairs$receptor_vector[[i]])) ){


      ###True expression from the region ########
      y1 <- t(gene_spot_df[c(matching_L_R_pairs$receptor_vector[[i]][r]),region_spot_IDs])
      x1 <- cellPropMatrix[region_spot_IDs, ]

      if (length(region_spot_IDs) < 2){
        truth <- rep(0, length(colnames(spot_cell_prop_df)))
        names(truth) <- colnames(spot_cell_prop_df)
      }else{
        mod1 <- nnls(x1, y1)
        truth <- mod1$x
        names(truth) <- colnames(spot_cell_prop_df)
        truth_f <- truth
        threshols2 <-  1/sum(truth > 0 )
        truth_f[ truth/sum(truth) < threshols2] <- 0

      }



      avg_receptor_list[[r]] <- truth
      avg_receptor_list_f[[r]] <- truth_f
    }
    # Calculate the product of corresponding elements across all vectors
    product_vector <- Reduce("*", avg_receptor_list)
    product_vector_f <- Reduce("*", avg_receptor_list_f)
    # Calculate the nth root of the product, where n is the number of vectors
    n <- length(avg_receptor_list)
    null_avg_receptor <- product_vector^(1/n)
    null_avg_receptor_f <- product_vector_f^(1/n)

    t <- outer(null_avg_ligand, null_avg_receptor, `*`)
    t_f <- outer(null_avg_ligand_f, null_avg_receptor_f, `*`)


    #prop <- spot_cell_prop_df
    prop_t <- outer(colSums(spot_cell_prop_df)/nrow(spot_cell_prop_df), colSums(spot_cell_prop_df)/nrow(spot_cell_prop_df), `*`)
    transformed_t <- t / (0.5 + t)
    transformed_t_f <- t_f / (0.5 + t_f)
    null_expression <- as.matrix(transformed_t * prop_t )
    null_expression_f <- as.matrix(transformed_t_f * prop_t )
    #null_expression <- as.matrix(transformed_t  )

    ################# Permutation ##############
    LigandVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$ligand_vector[[i]])))
    ReceptorVectorIndex <- c(which(rownames(gene_spot_df) %in% c(matching_L_R_pairs$receptor_vector[[i]])))
    # Call the function
    ## global1 with null_expression_f works
    averageResult <- .Call("_SpaCCI_Local_Regional_Permutations",
                           permutationMatrix,
                           permut_col,
                           cellPropMatrix,
                           GeneSpotMatrix,
                           LigandVectorIndex,
                           ReceptorVectorIndex,
                           null_expression,
                           nboot)


    #averageResult <- Global_Permutations(permutationMatrix,
    #                                     permut_null_regionMatrix,
    #                                     permut_col,
    #                                     cellPropMatrix,
    #                                     GeneSpotMatrix,
    #                                     LigandVectorIndex,
    #                                     ReceptorVectorIndex,
    #                                     null_expression_f,
    #                                     nboot)



    colnames(averageResult) <- colnames(spot_cell_prop_df)
    rownames(averageResult) <- colnames(spot_cell_prop_df)
    #p_vector <- as.vector(averageResult)

    #t <- outer(null_avg_ligand, null_avg_receptor, `*`) # row are ligand, col are receptor
    non_zero_indices <- which(null_expression != 0,arr.ind = TRUE)
    nul <- data.frame(
      Cell_type_Ligand = names(null_avg_ligand)[non_zero_indices[,1]],
      Cell_type_Receptor = names(null_avg_receptor)[non_zero_indices[,2]]
    )
    re <- averageResult[non_zero_indices]
    strength <- null_expression[non_zero_indices]
    nul$strength <- strength
    nul$PValue <- re
    nul$adjusted.PValue <- p.adjust(nul$PValue, method = "BH")
    nul$adjusted.PValue <- nul$PValue

    output[[i]] <- nul

  }


  close(con = pb)

  message("writing data frame")
  ##############

  table <- lapply(seq_along(output), function(i) {
    df <- output[[i]]
    interaction_info <- matching_L_R_pairs_info[i, ,drop=FALSE]
    interaction_expanded <- interaction_info[rep(1, nrow(df)), ,drop=FALSE]
    rownames(interaction_expanded) <- NULL
    df <- cbind(df, interaction_expanded)
    return(df)
  })
  dataframe <- do.call(rbind, table)


  return(list(pvalue_df = dataframe ))
}




