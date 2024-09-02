###############################################################################################
#' Plot SpaCCI Localized Interaction Results
#'
#' This function provides a unified interface to visualize the localized cell-cell interaction patterns inferred by SpaCCI, either using a Seurat object with a spatial image or a spatial coordinates data frame.
#'
#' @param Seurat_Object Optional. A Seurat object containing spatial data. If provided, the function will plot the interaction patterns on the tissue image.
#' @param spatial_coordinates_dataframe Optional. A data frame containing the spatial coordinates of the spots. The columns should include \code{"Spot_ID"}, \code{"imagerow"}, and \code{"imagecol"}. The row names must be the names of \code{"Spot_ID"}, matching those in the cell type proportion data frame or the gene expression data frame.
#' @param SpaCCI_local_Result_List A list containing the results from a SpaCCI local analysis. This list should include \code{dataframelist} and \code{RegionIDs_matrix}, which are the outputs from \code{run_SpaCCI(..., analysis_scale = "local",...)}.
#' @param Ligand_cell_type The name of the ligand cell type to plot. This should match the cell type names used in the \code{run_SpaCCI} analysis.
#' @param Receptor_cell_type The name of the receptor cell type to plot. This should match the cell type names used in the \code{run_SpaCCI} analysis.
#' @param spot_plot_size A numeric value controlling the size of the spots in the plot.
#' @param specific_LR_pair_name Optional. The name of a specific ligand-receptor pair to plot. If provided, the plot will focus on this interaction. The name should match those in the \code{SpaCCI_local_Result_List$dataframelist}.
#' @param significant_cutoff A numeric value specifying the significance cutoff for the adjusted P-values from the permutation test. Default is \code{0.05}.
#'
#' @return A plot object showing the localized interaction patterns. The plot will be generated using either the Seurat object or the spatial coordinates data frame, depending on the input provided.
#'
#' @examples
#' \dontrun{
#' # Plot localized SpaCCI results using Seurat object
#' Result <- run_SpaCCI(..., analysis_scale = "local",...)
#' plot_SpaCCI_local(Seurat_Object = seurat_object,
#'                   SpaCCI_local_Result_List = Result,
#'                   Ligand_cell_type = "Beta_cells",
#'                   Receptor_cell_type = "T_cells",
#'                   spot_plot_size = 3)
#'
#' # Plot localized SpaCCI results using spatial coordinates
#' plot_SpaCCI_local(spatial_coordinates_dataframe = spatial_coords,
#'                   SpaCCI_local_Result_List = Result,
#'                   Ligand_cell_type = "Beta_cells",
#'                   Receptor_cell_type = "T_cells",
#'                   spot_plot_size = 3)
#' }
#'
#' @export
#''

plot_SpaCCI_local <- function(Seurat_Object = NULL,
                              spatial_coordinates_dataframe = NULL,
                              SpaCCI_local_Result_List,
                              Ligand_cell_type,
                              Receptor_cell_type,
                              spot_plot_size,
                              specific_LR_pair_name = NULL,
                              significant_cutoff = 0.05){

  if (is.null(Seurat_Object) & is.null(spatial_coordinates_dataframe)){
    stop("Please input either a Seurat_object with image or input a spatial_coordinates_dataframe")
  }else if (!is.null(Seurat_Object)){
    print("plotting using Seurat image")
    local_plot <- plot_localized_Seurat(Seurat_object = Seurat_Object,
                                        resultdf_list = SpaCCI_local_Result_List$dataframelist,
                                        RegionIDs_matrix = SpaCCI_local_Result_List$RegionIDs_matrix,
                                        celltype_ligand = Ligand_cell_type,
                                        celltype_receptor = Receptor_cell_type,
                                        plot_size = spot_plot_size ,
                                        L_R_pair_name = specific_LR_pair_name,
                                        alpha = significant_cutoff)

  }else if(!is.null(spatial_coordinates_dataframe)){
    print("plotting using image spatial coordinates")
    local_plot <- plot_localized(spatial_coord = spatial_coordinates_dataframe ,
                                 resultdf_list = SpaCCI_local_Result_List$dataframelist,
                                 RegionIDs_matrix = SpaCCI_local_Result_List$RegionIDs_matrix,
                                 celltype_ligand = Ligand_cell_type,
                                 celltype_receptor = Receptor_cell_type,
                                 plot_size = spot_plot_size ,
                                 L_R_pair_name = specific_LR_pair_name,
                                 alpha = significant_cutoff)

  }

  return(local_plot)

}


# Plotting
#' Plot Localized Hotspot Pattern on Seurat Object
#'
#' Visualize the inferred cell-cell interaction localized pattern on the tissue image with Seurat_object
#' @param Seurat_object A Seurat object
#' @param resultdf_list A result of data frame list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`dataframelist`}
#' @param RegionIDs_matrix A result of matrix list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`RegionIDs_matrix`}
#' @param celltype_ligand Ligand cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.
#' @param celltype_receptor Receptor cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.
#' @param plot_size As this function incorporate with \code{Seurat}'s \code{`SpatialFeaturePlot`}, this parameter could control the plotting size of the each spot.
#' @param L_R_pair_name Initially this is set to \code{NULL}, if one is interested in a specific Ligand-Receptor pair, then one could specify the L_R_pair_name here. Note: the input name should match the L-R pair name exists in the dataframe in the output of SpaCCI_local "dataframelist".
#' @param alpha This is the significant cutoff for the adjusted-p-value of thr permutation test. Initially this is set to \code{0.05}, one could adjust the cutoff.
#'
#' @importFrom Seurat SpatialFeaturePlot
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom dplyr filter select mutate
#'
#' @return The localized plot from the inferred cell-cell interaction on the local scale.
#'
#' #' @examples
#' \dontrun{
#' # Run localized hotspot plot
#' Result <- run_SpaCCI(..., analysis_scale = "local",...)
#' local_plot <- plot_localized_Seurat(Seurat_object = gene_spot_df,
#'                                     resultdf_list = Result$dataframelist,
#'                                     RegionIDs_matrix = Result$RegionIDs_matrix,
#'                                     celltype_ligand = "Beta_cells",
#'                                     celltype_receptor = "T_ells",
#'                                     plot_size = 3)
#' }
#'
#' @export
#'
plot_localized_Seurat <- function(Seurat_object,
                              resultdf_list,
                              RegionIDs_matrix,
                              celltype_ligand,
                              celltype_receptor,
                              plot_size,
                              L_R_pair_name = NULL, alpha = 0.05){

  centerIDs <- names(resultdf_list)
  inter_name <- lapply(resultdf_list, function(df){ g<- unique(df$interaction_name)})
  inter_name <- unique(unlist(inter_name))
  # Assuming 'resultdf_list' is your list of data frames
  if (celltype_ligand == celltype_receptor){
    cleaned_dfs <- lapply(resultdf_list, function(df) {
      df <- df[which(df$Cell_type_Ligand == celltype_ligand & df$Cell_type_Receptor == celltype_receptor),]
      nameL <- paste("Null_avg_L", celltype_ligand)
      nameR <- paste("Null_avg_R", celltype_receptor)

      if (!is.null(L_R_pair_name)){
        df_clean <- df[which( df$interaction_name %in% L_R_pair_name ), ]
      }else {
        df_clean <- df
      }
      # Return the cleaned data frame
      return(df_clean)
    })
  }else{
    cleaned_dfs <- lapply(resultdf_list, function(df) {

      df_clean <- df[which(df$Cell_type_Ligand == celltype_ligand & df$Cell_type_Receptor == celltype_receptor),]
      if (!is.null(L_R_pair_name)){
        df_clean <- df_clean[which( df_clean$interaction_name %in% L_R_pair_name ), ]
      }
      #df_clean <- df_clean[which(df_clean[[nameL]] != 0 & df_clean[[nameR]] != 0), ]
      # Return the cleaned data frame
      return(df_clean)
    })
  }

  # For each data frame, count the rows where 'adjusted.PValue' is less than 0.05
  count_rows <- sapply(cleaned_dfs, function(df) nrow(df[df$adjusted.PValue < alpha, ]))
  g <- as.data.frame(count_rows)
  Seurat_object@meta.data$CCIcount <- 0
  ############################## option 1
  for (row_name in rownames(g)) {
    if (row_name %in% rownames(Seurat_object@meta.data)) {
      Seurat_object@meta.data[row_name, "CCIcount"] <- g[row_name, "count_rows"]
    }
  }
  if (!is.null(L_R_pair_name)){
    max_value <- length(L_R_pair_name)
  }else{max_value <- ceiling(max(Seurat_object@meta.data$CCIcount))}
  min_value <- 0

  #p_center_only <- SpatialFeaturePlot(Seurat_object, features = "CCIcount",pt.size.factor = plot_size) + ggplot2::scale_fill_gradient2(low = "black", high = "yellow", mid = "purple", limits =c(0,max_value), midpoint = max_value/2 , space = "Lab")
  df <- Seurat_object@meta.data
  df <- df[,c("CCIcount"),drop=FALSE]

  reg <- data.frame(nrow = (rownames(Seurat_object@meta.data)), ncol = 0 )
  rownames(reg) <- reg$nrow
  reg$CCIcount <- 0
  for (row_name in rownames(g)) {
    if (row_name %in% rownames(reg)) {
      reg[row_name, "CCIcount"] <- g[row_name, "count_rows"]
    }
  }

  # Add a column in reg for each centerID and initialize it to 0
  for (id in centerIDs) { reg[[paste0(id)]] <- NA}
  # Iterate over centerIDs and update the corresponding column in reg
  for (id in centerIDs) {
    # Get the corresponding IDs from RegionIDs_matrix
    IDs <- RegionIDs_matrix[[id]]
    # Update the corresponding column in reg for rows that match the IDs
    reg[rownames(reg) %in% IDs , id] <- g[id, ]
  }

  row_means <- rowMeans(reg[, 4:ncol(reg)], na.rm = TRUE)
  reg$ncol <- row_means
  reg$ncol[is.na(reg$ncol)] <- 0
  Seurat_object@meta.data$CCIcount <- reg$ncol
  Seurat_object@meta.data$CCI_Probability <- reg$ncol

  min_value <- 0

  max_value <-  ceiling(max(Seurat_object@meta.data$CCIcount))

  p_all_average <- SpatialFeaturePlot(Seurat_object, features = "CCIcount",pt.size.factor = plot_size)  + ggplot2::scale_fill_gradient2(low = "black", high = "yellow", mid = "purple3", limits =c(0,max_value),midpoint = max_value/2, space = "Lab")


  return(p_all_average)

}


#' Plot Localized Hotspot Pattern
#'
#' Visualize the inferred cell-cell interaction localized pattern if NOT using Seurat Object
#'
#' @param spatial_coord A data frame of the spatial coordinates. The columns should include \code{"Spot_ID"}, \code{"imagerow"}, and \code{"imagecol"}. And the row names must be the names of \code{"Spot_ID"}, which is the same as the rownames in cell type proportion data frame or the colnames of the gene* spot expression data frame
#' @param resultdf_list A result of data frame list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`dataframelist`}
#' @param RegionIDs_matrix A result of matrix list from the output of \code{run_SpaCCI(..., analysis_scale = "local",...)} \code{`RegionIDs_matrix`}
#' @param celltype_ligand Ligand cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.
#' @param celltype_receptor Receptor cell type string inputted by user, the name of the cell type should match the names in the \code{`spot_cell_proportion_dataframe`} during the \code{run_SpaCCI} analysis.
#' @param plot_size As this function incorporate with \code{Seurat}'s \code{`SpatialFeaturePlot`}, this parameter could control the plotting size of the each spot.
#' @param L_R_pair_name Initially this is set to \code{NULL}, if one is interested in a specific Ligand-Receptor pair, then one could specify the L_R_pair_name here. Note: the input name should match the L-R pair name exists in the dataframe in the output of SpaCCI_local "dataframelist".
#' @param alpha This is the significant cutoff for the adjusted-p-value of thr permutation test. Initially this is set to \code{0.05}, one could adjust the cutoff.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient2 labs theme
#' @importFrom dplyr filter select mutate
#'
#' @return The localized plot from the inferred cell-cell interaction on the local scale.
#'
#' @examples
#' \dontrun{
#' # Run localized hotspot plot
#' Result <- run_SpaCCI(..., analysis_scale = "local",...)
#' local_plot <- plot_localized(Seurat_object = gene_spot_df,
#'                             resultdf_list = Result$dataframelist,
#'                             RegionIDs_matrix = Result$RegionIDs_matrix,
#'                             celltype_ligand = "Beta_cells",
#'                             celltype_receptor = "T_ells",
#'                             plot_size = 3)

#' }
#'
#' @export
#'

plot_localized <- function(spatial_coord,
                           resultdf_list,
                           RegionIDs_matrix,
                           celltype_ligand,
                           celltype_receptor,
                           plot_size,
                           L_R_pair_name = NULL, alpha = 0.05){

  centerIDs <- names(resultdf_list)
  inter_name <- lapply(resultdf_list, function(df){ g<- unique(df$interaction_name)})
  inter_name <- unique(unlist(inter_name))
  # Assuming 'resultdf_list' is your list of data frames
  if (celltype_ligand == celltype_receptor){
    cleaned_dfs <- lapply(resultdf_list, function(df) {
      df <- df[which(df$Cell_type_Ligand == celltype_ligand & df$Cell_type_Receptor == celltype_receptor),]
      nameL <- paste("Null_avg_L", celltype_ligand)
      nameR <- paste("Null_avg_R", celltype_receptor)

      if (!is.null(L_R_pair_name)){
        df_clean <- df[which( df$interaction_name %in% L_R_pair_name ), ]
      }else {
        df_clean <- df
      }
      # Return the cleaned data frame
      return(df_clean)
    })
  }else{
    cleaned_dfs <- lapply(resultdf_list, function(df) {

      df_clean <- df[which(df$Cell_type_Ligand == celltype_ligand & df$Cell_type_Receptor == celltype_receptor),]
      if (!is.null(L_R_pair_name)){
        df_clean <- df_clean[which( df_clean$interaction_name %in% L_R_pair_name ), ]
      }
      # Return the cleaned data frame
      return(df_clean)
    })
  }

  # For each data frame, count the rows where 'adjusted.PValue' is less than 0.05
  count_rows <- sapply(cleaned_dfs, function(df) nrow(df[df$adjusted.PValue < alpha, ]))
  g <- as.data.frame(count_rows)
  spatial_coord$CCIcount <- 0
  ############################## option 1
  for (row_name in rownames(g)) {
    if (row_name %in% rownames(spatial_coord)) {
      spatial_coord[row_name, "CCIcount"] <- g[row_name, "count_rows"]
    }
  }

  df <- spatial_coord
  df <- df[,c("CCIcount"),drop=FALSE]

  reg <- data.frame(nrow = (rownames(spatial_coord)), ncol = 0 )
  rownames(reg) <- reg$nrow
  reg$CCIcount <- 0
  for (row_name in rownames(g)) {
    if (row_name %in% rownames(reg)) {
      reg[row_name, "CCIcount"] <- g[row_name, "count_rows"]
    }
  }

  for (id in centerIDs) { reg[[paste0(id)]] <- NA}
  # Iterate over centerIDs and update the corresponding column in reg
  for (id in centerIDs) {
    # Get the corresponding IDs from RegionIDs_matrix
    IDs <- RegionIDs_matrix[[id]]
    # Update the corresponding column in reg for rows that match the IDs
    reg[rownames(reg) %in% IDs , id] <- g[id, ]
  }


  row_means <- rowMeans(reg[, 4:ncol(reg)], na.rm = TRUE)
  reg$ncol <- row_means
  reg$ncol[is.na(reg$ncol)] <- 0
  spatial_coord$CCIcount <- reg$ncol

  min_value <- 0

  max_value <-  ceiling(max(spatial_coord$CCIcount ))


  p_all_average <- ggplot(spatial_coord, aes(x = imagecol, y = imagerow, colour = CCIcount)) +
    geom_point(alpha = 0.9, size = plot_size, stroke = 0, shape = 16) +
    #scale_colour_gradient(low = "#F0F0F0", high = '#F29403') +
    scale_color_gradient2(low = "black", mid = "purple", high = "yellow", limits = c(0, max_value), midpoint = max_value / 2, space = "Lab") +
    labs(color = "CCI Count", x = "X Center", y = "Y Center") +
    xlim((floor(min(spatial_coord$imagecol)/25)*25), (ceiling(max(spatial_coord$imagecol)/25)*25) ) + ylim((floor(min(spatial_coord$imagerow)/25)*25), (ceiling(max(spatial_coord$imagerow)/25)*25) ) +
    theme(
      plot.margin = margin(0.35, 0.35, 0.35, 0.35, "cm"),
      panel.background = element_rect(colour = "white", fill = "white"),
      plot.background = element_rect(colour = "white", fill = "white"),
      panel.border = element_rect(colour = "grey39", fill = NA, size = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 4),
      legend.position = 'top',
      legend.key = element_rect(colour = "transparent", fill = "white"),
      legend.key.size = unit(0.5, 'cm')
    )


  return(p_all_average)


}


#' Plot SpaCCI Results on the heatmap
#' Visualize inferred significant cell-cell interactions using a heatmap
#' @param SpaCCI_Result_List A list containing the results from a SpaCCI \code{"regional"} or \code{"global"} analysis. This list should include \code{pvalue_df}, which are the outputs from \code{run_SpaCCI(..., analysis_scale = "regional",...)} or  \code{run_SpaCCI(..., analysis_scale = "global",...)}.
#' @param specific_celltypes A vector of cell types to include in the heatmap, i.e c("Celltype_A","Celltype_B"). NOTE: the cell type names should match the names input in the SpaCCI analysis.
#' @param pathways A vector of pathways to filter the interactions. Initially set to \code{NULL}, if not, then it will aggregate the results of the selected pathways.
#' @param interaction A vector of interactions to filter. Initially set to \code{NULL}, if not, then it will aggregate the results of the selected interactions.
#' @param log1p_transform Logical; whether to apply a log(1 + x) transformation to the count matrix.
#' @param show_rownames Logical; whether to show row names in the heatmap.
#' @param show_colnames Logical; whether to show column names in the heatmap.
#' @param scale Character; whether to scale the data ("row", "column", "none").
#' @param cluster_cols Logical; whether to cluster columns.
#' @param cluster_rows Logical; whether to cluster rows.
#' @param border_color Character; color of the heatmap borders.
#' @param fontsize_row Numeric; font size for row names.
#' @param fontsize_col Numeric; font size for column names.
#' @param family Character; font family for text in the heatmap.
#' @param main Character; title of the heatmap.
#' @param treeheight_col Numeric; height of the column dendrogram.
#' @param treeheight_row Numeric; height of the row dendrogram.
#' @param low_col Character; color for low values in the heatmap.
#' @param mid_col Character; color for mid values in the heatmap.
#' @param high_col Character; color for high values in the heatmap.
#' @param alpha Numeric; significance threshold for p-values, initailly set to \code{0.05}.
#' @param return_tables Logical; whether to return the count matrix and summary tables.
#' @param symmetrical Logical; whether to make the heatmap symmetrical.
#' @param ... Additional arguments passed to `pheatmap`.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr filter group_by summarise left_join
#' @importFrom reshape2 acast
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @importFrom viridis viridis
#'
#' @return If `return_tables` is \code{FALSE} (default), the function returns a heatmap object created by \code{pheatmap}, showing the count of significant cell-cell interactions. If `return_tables` is \code{TRUE}, the function returns a list containing:
#' \describe{
#'   \item{\strong{heatmap}}{The heatmap object showing the significant cell-cell interactions.}
#'   \item{\strong{heatmap_countmatrix}}{The matrix used to generate the heatmap, with cell types as rows and columns, and counts of significant interactions as values.}
#'   \item{\strong{table}}{A data frame summarizing the counts of significant interactions between each ligand and receptor cell type combination.}
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Run global analysis
#' result_global <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
#'                             spot_cell_proportion_dataframe = cell_prop_df,
#'                             spatial_coordinates_dataframe = spatial_coords_df,
#'                             LR_database_list = LR_database("Human", "CellChat", gene_spot_df),
#'                             analysis_scale = "global")
#'
#' # Run the result heatmap for global analysis
#' plot_SpaCCI_heatmap(SpaCCI_Result_List = result_global,
#'                     symmetrical = FALSE, cluster_cols = FALSE, return_tables = FALSE,
#'                     cluster_rows = FALSE, #cellheight = 10, cellwidth = 10,
#'                     specific_celltypes = c(colnames(cell_prop_df)),
#'                     main= "Cell-Cell Interaction Count")
#'
#' # Run regional analysis
#' result_regional <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
#'                               spot_cell_proportion_dataframe = cell_prop_df,
#'                               spatial_coordinates_dataframe = spatial_coords_df,
#'                               LR_database_list = LR_database("Human", "CellChat", gene_spot_df),
#'                               analysis_scale = "regional",
#'                               region_spot_IDs = c("Spot1", "Spot2", "Spot3"))
#'
#' # Run the result heatmap for regional analysis
#' plot_SpaCCI_heatmap(SpaCCI_Result_List = result_regional,
#'                     symmetrical = FALSE, cluster_cols = FALSE, return_tables = FALSE,
#'                     cluster_rows = FALSE, #cellheight = 10, cellwidth = 10,
#'                     specific_celltypes = c(colnames(cell_prop_df)),
#'                     main= "Cell-Cell Interaction Count")
#'
#'}
#'
#' @export
#'

plot_SpaCCI_heatmap <- function(SpaCCI_Result_List , specific_celltypes = NULL, pathways =NULL, interaction=NULL , log1p_transform = FALSE,
                                show_rownames = TRUE, show_colnames = TRUE, scale = "none", cluster_cols = TRUE,
                                cluster_rows = TRUE, border_color = "white", fontsize_row = 11, fontsize_col = 11,
                                family = "Arial", main = "", treeheight_col = 0, treeheight_row = 0, low_col = "dodgerblue4",
                                mid_col = "peachpuff", high_col = "deeppink4", alpha = 0.05, return_tables = FALSE,
                                symmetrical = FALSE, ...) {
  suppressWarnings({
  requireNamespace("reshape2")
  requireNamespace("grDevices")

  all_intr <- SpaCCI_Result_List$pvalue_df
  if (!is.null(pathways) & !is.null(interaction) ){
    all_intr <- all_intr[which(all_intr$pathway_name %in% c(pathways) & all_intr$interaction_name %in% c(interaction) ),]
  }else if(!is.null(pathways) & is.null(interaction) ){
    all_intr <- all_intr[which(all_intr$pathway_name %in% c(pathways) ),]
  }else if (is.null(pathways) & !is.null(interaction) ){
    all_intr <- all_intr[which(all_intr$interaction_name %in% c(interaction) ),]
  }


  if (!is.null(specific_celltypes)){
    all_intr <- all_intr[which(all_intr$Cell_type_Ligand %in% c(specific_celltypes) & all_intr$Cell_type_Receptor %in% c(specific_celltypes) ),]
  }else{
    all_intr <- all_intr
  }
  all_intr <- all_intr[, c("Cell_type_Ligand" ,"Cell_type_Receptor","adjusted.PValue" )] # remain only the data of p-values

  all_combinations <- expand.grid(
    Cell_type_Ligand = unique((specific_celltypes)),
    Cell_type_Receptor = unique((specific_celltypes))
  )
  all_intr$significant <- all_intr$adjusted.PValue < alpha

  all_count <- all_intr %>%
    filter(significant) %>%
    group_by(Cell_type_Ligand, Cell_type_Receptor) %>%
    summarise(COUNT = n(),.groups = 'drop')
  all_count <- left_join(all_combinations, all_count, by = c("Cell_type_Ligand", "Cell_type_Receptor"))
  all_count$COUNT[is.na(all_count$COUNT)] <- 0

  if (any(all_count$COUNT) > 0) {
    count_mat <- reshape2::acast(Cell_type_Ligand ~ Cell_type_Receptor, data = all_count, value.var = "COUNT")
    count_mat[is.na(count_mat)] <- 0
    count_mat <- count_mat[specific_celltypes, specific_celltypes]
    col.heatmap <- (grDevices::colorRampPalette(c(low_col, mid_col, high_col)))(1000)
    if (symmetrical) {
      dcm <- diag(count_mat)
      count_mat <- count_mat + t(count_mat)
      diag(count_mat) <- dcm
    }

    if (log1p_transform == TRUE) {
      count_mat <- log1p(count_mat)
    }

    p <- pheatmap(count_mat, show_rownames = show_rownames, show_colnames = show_colnames,
                  scale = scale, cluster_cols = cluster_cols, border_color = border_color,
                  cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                  main = main, treeheight_row = treeheight_row, family = family, color = col.heatmap,
                  treeheight_col = treeheight_col, ...)
    if (return_tables) {
      return(list(heatmap = p, heatmap_countmatrix = count_mat, table = all_count))
    } else {
      return(p)
    }
  } else {
    stop("There are no significant results using p-value of: ", alpha, call. = FALSE)
  }

  })

}



#' Generate a Color Palette
#'
#' This function generates a color palette. It selects colors from a predefined color space, and if more colors are needed than are available in the predefined set, it generates a palette using color interpolation.
#'
#' @param n An integer specifying the number of colors needed.
#'
#' @return A character vector of colors in hexadecimal format.
#'
#' @examples
#' # Generate a palette with 5 colors
#' palette <- scPalette(5)
#' print(palette)
#'
#' # Generate a palette with 30 colors
#' large_palette <- scPalette(30)
#' print(large_palette)
#'
#' @importFrom grDevices colorRampPalette
#' @export
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#DCF0B9','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}



#' Plot SpaCCI Results on Chord Diagram
#'
#' This function generates a chord diagram to visualize cell-cell interactions based on ligand-receptor pairs. The interactions can be filtered by specific cell types, pathways, or interaction names.
#'
#' @param SpaCCI_Result_List A list containing the results from a SpaCCI \code{"regional"} or \code{"global"} analysis. This list should include \code{pvalue_df}, which are the outputs from \code{run_SpaCCI(..., analysis_scale = "regional",...)} or  \code{run_SpaCCI(..., analysis_scale = "global",...)}.
#' @param specific_celltypes A character vector specifying the cell types to include in the plot, RECOMMEND using colnames of cell type proportion matrix to include all cell types. If \code{NULL}, cell types that involved in significant interactions are included.
#' @param pathway_name A single character string specifying the pathway name to filter the interactions. If \code{NULL}, all pathways are included.
#' @param L_R_pair_name A character vector specifying the ligand-receptor pair names to include in the plot. If \code{NULL}, all interactions are included.
#' @param color A named vector of colors to use for the cell types. If \code{NULL}, a default color palette is used.
#' @param alpha A numeric value specifying the significance threshold for adjusted P-values. Initially, set to \code{0.05}.
#'
#' @return A chord diagram plot visualizing the significant cell-cell interactions.
#'
#' @examples
#' \dontrun{
#'
#' # Run global analysis
#' result_global <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
#'                             spot_cell_proportion_dataframe = cell_prop_df,
#'                             spatial_coordinates_dataframe = spatial_coords_df,
#'                             LR_database_list = LR_database("Human", "CellChat", gene_spot_df),
#'                             analysis_scale = "global")
#'
#' # Run the result heatmap for global analysis
#' plot_SpaCCI_chordDiagram(SpaCCI_Result_List = result_global,
#'                          specific_celltypes = c(colnames(cell_prop_df)),
#'                          L_R_pair_name  = "AREG_EGFR")
#'
#' # Run regional analysis
#' result_regional <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
#'                               spot_cell_proportion_dataframe = cell_prop_df,
#'                               spatial_coordinates_dataframe = spatial_coords_df,
#'                               LR_database_list = LR_database("Human", "CellChat", gene_spot_df),
#'                               analysis_scale = "regional",
#'                               region_spot_IDs = c("Spot1", "Spot2", "Spot3"))
#'
#' # Run the result heatmap for regional analysis
#' plot_SpaCCI_chordDiagram(SpaCCI_Result_List = result_regional,
#'                          specific_celltypes = c(colnames(cell_prop_df)),
#'                          L_R_pair_name  = "AREG_EGFR")
#'
#'}
#'
#' @importFrom circlize chordDiagram
#' @importFrom dplyr filter group_by summarise left_join
#' @importFrom reshape2 acast
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics title
#'
#' @export
plot_SpaCCI_chordDiagram <- function(SpaCCI_Result_List, specific_celltypes = NULL, pathway_name = NULL, L_R_pair_name = NULL, color = NULL, alpha = 0.05) {
  suppressWarnings({
  all_count <- SpaCCI_Result_List$pvalue_df

  if (!is.null(specific_celltypes)) {
    if (!is.null(color)) {
      colors <- color
    } else {
      colors <- scPalette(length(specific_celltypes))
      names(colors) <- specific_celltypes
    }

    all_combinations <- expand.grid(
      Cell_type_Ligand = specific_celltypes,
      Cell_type_Receptor = specific_celltypes
    )
    all_count <- all_count[all_count$Cell_type_Ligand %in% specific_celltypes & all_count$Cell_type_Receptor %in% specific_celltypes, ]
  } else {
    if (!is.null(color)) {
      colors <- color
    } else {
      colors <- scPalette(length(unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor))))
      names(colors) <- unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor))
    }
    all_combinations <- expand.grid(
      Cell_type_Ligand = unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor)),
      Cell_type_Receptor = unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor))
    )
    all_count <- all_count[all_count$Cell_type_Ligand %in% unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor)) & all_count$Cell_type_Receptor %in% unique(c(all_count$Cell_type_Ligand, all_count$Cell_type_Receptor)), ]
  }

  if (!is.null(pathway_name) & is.null(L_R_pair_name)) {
    all_count <- all_count[all_count$pathway_name == pathway_name, ]
    if (nrow(all_count) == 0) {
      stop("Please Enter Correct Name of Pathway")
    }
  } else if (!is.null(pathway_name) & !is.null(L_R_pair_name)) {
    all_count <- all_count[all_count$interaction_name %in% L_R_pair_name, ]
    if (nrow(all_count) == 0) {
      stop("Please Enter Correct Name of Pathway or L-R Pair")
    }
  } else if (is.null(pathway_name) & !is.null(L_R_pair_name)) {
    all_count <- all_count[all_count$interaction_name %in% L_R_pair_name, ]
    if (nrow(all_count) == 0) {
      stop("Please Enter Correct Name of L-R Pair")
    }
  }

  all_count$significant <- all_count$adjusted.PValue < alpha
  all_count <- all_count %>%
    filter(significant) %>%
    group_by(Cell_type_Ligand, Cell_type_Receptor) %>%
    summarise(COUNT = n(), .groups = 'drop')

  all_count <- left_join(all_combinations, all_count, by = c("Cell_type_Ligand", "Cell_type_Receptor"))
  all_count$Cell_type_Ligand <- factor(all_count$Cell_type_Ligand, levels = unique(all_count$Cell_type_Ligand))
  all_count$Cell_type_Receptor <- factor(all_count$Cell_type_Receptor, levels = unique(all_count$Cell_type_Ligand))
  all_count$COUNT[is.na(all_count$COUNT)] <- 0

  if (any(all_count$COUNT) > 0) {
    count_mat <- reshape2::acast(Cell_type_Ligand ~ Cell_type_Receptor, data = all_count, value.var = "COUNT")
    count_mat[is.na(count_mat)] <- 0
  } else {
    stop("There are no significant interactions")
  }

  mat <- count_mat + sum(count_mat) / (ncol(count_mat) * ncol(count_mat))
  df <- data.frame(from = rep(rownames(mat), times = ncol(mat)),
                   to = rep(colnames(mat), each = nrow(mat)),
                   value = as.vector(mat),
                   stringsAsFactors = FALSE)

  chordDiagram(df, annotationTrack = c("name", "grid"),
               grid.col = colors,
               link.arr.type = "big.arrow",
               link.arr.length = 0.08,
               directional = 1,
               self.link = 2,
               diffHeight = -0.02,
               direction.type = c("diffHeight", "arrows"),
               annotationTrackHeight = c(0.14, 0.07),
               link.visible = df[[3]] >= 1, scale = TRUE, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df))))))

  if (!is.null(L_R_pair_name) & is.null(pathway_name)) {
    title(paste("Cell-Cell Interaction of", L_R_pair_name))
  } else if (!is.null(pathway_name)) {
    title(paste("Cell-Cell Interaction of", pathway_name, "pathway"))
  } else {
    title("Cell-Cell Interaction")
  }

  })
}










