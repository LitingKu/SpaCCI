###############################################################################################
# database for deriving ligand-receptor pairs
#' Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' This function identifies possible ligand-receptor (L-R) pairs for cell-cell communication analysis using a selected database. It checks for the presence of all genes involved in each L-R pair within the provided gene expression matrix, filtering based on a specified expression percentage threshold. The function supports multiple databases including CellChat, CellPhoneDB, Cellinker, ICELLNET, and ConnectomeDB.
#'
#' @param species A string specifying the species (\code{"Human"} or \code{"Mouse"}).
#' @param database_name A string specifying the L-R database to use. Options include \code{"CellChat"}, \code{"CellPhoneDB"}, \code{"Cellinker"}, \code{"ICELLNET"}, and \code{"ConnectomeDB"}.
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is \code{10}, meaning the gene express over 10\% of spots.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{possible_LR_pairs}}{A data frame of L-R pairs where all genes are present in the \code{gene_spot_expression_dataframe} and meet the expression threshold. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{\code{possible_LR_pairs_info}}{A data frame with detailed information about the identified L-R pairs, including their original annotations from the selected database.}
#' }
#'
#' @examples
#' library(SpaCCI)
#' #Load the example data
#' data(test_data)
#' gene_spot_df <- test_data$gene_spot_df
#' result <- LR_database(species = "Human",
#'             database_name = "CellChat",
#'             gene_spot_expression_dataframe = gene_spot_df)
#' @export
#'
#'
LR_database <- function(species, database_name, gene_spot_expression_dataframe, percentage = 10){
  if ( database_name %in% c("CellChat","cellchat","Cellchat","cellChat")  ){

    dat <- possible_L_R_pairs_cellchat(species, gene_spot_expression_dataframe, percentage)

  }else if(species %in% c("Mouse","mouse") & database_name %in% c("CellPhoneDB","CellphoneDB","cellphoneDB","cellPhoneDB") ){
    stop("CellPhoneDB only has database for Human, if want to use it for mouse, please use the ortholog human genes for our mouse dataset")
  }else if(species %in% c("Human","human") & database_name %in% c("CellPhoneDB","CellphoneDB","cellphoneDB","cellPhoneDB") ){
    dat <- possible_L_R_pairs_cellphoneDB(gene_spot_expression_dataframe, percentage)

  }else if( database_name %in% c("Cellinker","cellinker")) {
    dat <- possible_L_R_pairs_Cellinker(species,  gene_spot_expression_dataframe, percentage)

  }else if(species %in% c("Mouse","mouse") & database_name %in% c("ICELLNET","icellnet") ){
    stop("ICELLNET only has database for Human, if want to use it for mouse, please use the ortholog human genes for our mouse dataset")
  }else if(species %in% c("Human","human") & database_name %in% c("ICELLNET","icellnet") ){
    dat <- possible_L_R_pairs_ICELLNET(gene_spot_expression_dataframe, percentage)

  }else if(species %in% c("Mouse","mouse") & database_name %in% c("ConnectomeDB","connectomeDB") ){
    stop("ConnectomeDB only has database for Human, if want to use it for mouse, please use the ortholog human genes for our mouse dataset")
  }else if(species %in% c("Human","human") & database_name %in% c("ConnectomeDB","connectomeDB") ){
    dat <- possible_L_R_pairs_connectome(gene_spot_expression_dataframe, percentage)

  }else {
    stop("Invalid species or database name provided.")
  }

  return(list( possible_LR_pairs = dat$possible_L_R_pairs, possible_LR_pairs_info = dat$possible_L_R_pairs_details) )

}



#'
#' CellChat Database: Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' This function identifies possible ligand-receptor (L-R) pairs for cell-cell communication analysis using a subset of the CellChat database. It checks for the presence of all genes involved in each L-R pair within the provided gene expression matrix.
#'
#' @param species A string specifying the species ("Human" or "Mouse"). The function selects the appropriate CellChatDB object, typically 'CellChatDB.human' or 'CellChatDB.mouse', which contains information on ligand-receptor interactions.
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{possible_L_R_pairs}{A data frame of L-R pairs where all genes are present in the 'gene_spot_expression_dataframe'. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{possible_L_R_pairs_details}{A data frame with detailed information about the L-R pairs, including the original annotations from the CellChatDB.}
#' }
#'
#' @examples
#'
#' \donttest{
#' library(SpaCCI)
#' #Load the example data
#' load(system.file("extdata", "Tutorial_example_data.rda", package = "SpaCCI"))
#' Example_Seurat <- NormalizeData(Example_Seurat)
#' gene_spot_df <- as.data.frame(Example_Seurat@assays$Spatial@data)
#' result <- possible_L_R_pairs_cellchat(CellChatDB.human,
#'       gene_spot_expression_dataframe = gene_spot_df)
#' }
#' @importFrom dplyr %>% filter select mutate arrange group_by summarise
#'
#' @export
#'
#'
possible_L_R_pairs_cellchat <- function(species, gene_spot_expression_dataframe, percentage) {
  lr_data_path <- system.file("extdata", "CellChatDB.rda", package = "SpaCCI")
  load(lr_data_path)

  if (species %in% c("Human","human") ){
    CellChatDB <- CellChatDB_human
  }else if (species %in% c("Mouse","mouse") ){
    CellChatDB <- CellChatDB_mouse
  }else{
    stop("Please enter species of either Human or Mouse")
  }

  CellChatDB.use <- CellChatDB[["interaction"]]
  DB_inter <- CellChatDB.use[which(CellChatDB.use$annotation %in% c("Secreted Signaling")),]

  # Function to split and trim whitespace
  split_and_trim <- function(x) {
    sapply(strsplit(x, ",", fixed = TRUE), function(y) trimws(y))
  }
  # Apply the function to ligand and receptor columns
  DB_inter$ligand_vector <- lapply(DB_inter$ligand.symbol, split_and_trim)
  DB_inter$receptor_vector <- lapply(DB_inter$receptor.symbol, split_and_trim)

  # Combine ligand and receptor genes into a single vector for each row
  DB_inter$combined_genes <- mapply(function(l, r) c(l, r),
                                    DB_inter$ligand_vector,
                                    DB_inter$receptor_vector,
                                    SIMPLIFY = FALSE)

  # Function to check if all genes of an LR pair are in the data
  check_genes_presence <- function(gene_vector, gene_data) {
    all(gene_vector %in% rownames(gene_data))
  }
  exp_gene_express <- gene_spot_expression_dataframe[rowSums(gene_spot_expression_dataframe[] )>0,]
  binary_expression_data <- ifelse(exp_gene_express > 0, 1, 0)
  expression_percent <- apply(binary_expression_data, 1, function(row) sum(row > 0) / length(row) * 100)
  expressed_genes <- names(expression_percent)[expression_percent >= percentage ] # 1% 5% 10%
  final_gene_spot_expression_dataframe <- exp_gene_express[expressed_genes,]

  presence_vector <- mapply(check_genes_presence, DB_inter$combined_genes, MoreArgs = list(gene_data = final_gene_spot_expression_dataframe))

  # Subset interactions to get only those L-R pairs where all genes are present in gene_data
  L_R_pairs <- DB_inter[presence_vector, c("ligand_vector", "receptor_vector", "combined_genes")]
  info_inter <- DB_inter[presence_vector, setdiff(names(DB_inter), c("ligand_vector", "receptor_vector", "combined_genes"))]

  return(list(possible_L_R_pairs = L_R_pairs, possible_L_R_pairs_details = info_inter))
}



#' Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' CellPhone Database: This function identifies possible ligand-receptor (L-R) pairs for cell-cell communication analysis using data from a CellPhoneDB dataset. It checks for the presence of all genes involved in each L-R pair within the provided gene expression matrix and filters based on a specified expression percentage threshold.
#'
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{possible_L_R_pairs}{A data frame of L-R pairs where all genes are present in the 'gene_spot_expression_dataframe' and meet the expression threshold. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{possible_L_R_pairs_details}{A data frame with detailed information about the identified L-R pairs, including their original annotations from the CellPhoneDB dataset.}
#' }
#'
#'
#' @importFrom dplyr rename %>% filter select mutate arrange group_by summarise
#' @importFrom utils read.csv
#'
#' @export
#'
possible_L_R_pairs_cellphoneDB <- function(gene_spot_expression_dataframe, percentage ){
  filepath <- system.file("extdata", "cellphoneDB_interactions.csv", package = "SpaCCI")
  cellphone <- read.csv(filepath)
  cpdb_LR <- cellphone[which(cellphone$directionality == "Ligand-Receptor"),]
  cpdb_LR <- cpdb_LR %>% rename(interaction_name = interactors, pathway_name = classification)
  split_interactions <- lapply(strsplit(cpdb_LR$interaction_name, "-"), function(x) {
    sapply(x, function(y) gsub("\\+", ",", y))
  })
  interaction_df <- data.frame(
    Ligand = sapply(split_interactions, function(x) x[1]),
    Receptor = sapply(split_interactions, function(x) x[2])
  )
  split_and_trim <- function(x) {
    sapply(strsplit(x, ",", fixed = TRUE), function(y) trimws(y))
  }
  # Apply the function to ligand and receptor columns
  cpdb_LR$ligand_vector <- lapply(interaction_df$Ligand, split_and_trim)
  cpdb_LR$receptor_vector <- lapply(interaction_df$Receptor, split_and_trim)
  #cpdb_LR$combined_genes <- mapply(c, cpdb_LR$ligand_vector, cpdb_LR$receptor_vector, SIMPLIFY = FALSE)
  cpdb_LR$combined_genes <- mapply(function(l, r) c(l, r),
                                   cpdb_LR$ligand_vector,
                                   cpdb_LR$receptor_vector,
                                   SIMPLIFY = FALSE)
  # function of checking if all genes of LR are in the data
  check_genes_presence <- function(gene_vector, gene_data) { all(gene_vector %in% rownames(gene_data))}

  exp_gene_express <- gene_spot_expression_dataframe[rowSums(gene_spot_expression_dataframe[] )>0,]
  binary_expression_data <- ifelse(exp_gene_express > 0, 1, 0)
  expression_percent <- apply(binary_expression_data, 1, function(row) sum(row > 0) / length(row) * 100)
  expressed_genes <- names(expression_percent)[expression_percent >= percentage ] # 1% 5% 10%
  final_gene_spot_expression_dataframe <- exp_gene_express[expressed_genes,]

  presence_vector <- mapply(check_genes_presence, cpdb_LR$combined_genes, MoreArgs = list(gene_data = final_gene_spot_expression_dataframe))
  # Subset inter to get only those L-R pairs where all genes are present in gene_data
  L_R_pairs <- cpdb_LR[presence_vector, c("ligand_vector", "receptor_vector","combined_genes" )]
  info_inter <- cpdb_LR[presence_vector, setdiff(names(cpdb_LR), c("ligand_vector", "receptor_vector", "combined_genes"))]
  return(list( possible_L_R_pairs = L_R_pairs, possible_L_R_pairs_details = info_inter))
}



#' Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' Cellinker Database: This function identifies possible ligand-receptor (L-R) pairs for cell-cell communication analysis using data from the Cellinker database. It checks for the presence of all genes involved in each L-R pair within the provided gene expression matrix, filtering based on a specified expression percentage threshold.
#'
#' @param species A string specifying the species ("Human" or "Mouse"). The function selects the appropriate Cellinker interaction file based on this input.
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{possible_L_R_pairs}{A data frame of L-R pairs where all genes are present in the 'gene_spot_expression_dataframe' and meet the expression threshold. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{possible_L_R_pairs_details}{A data frame with detailed information about the identified L-R pairs, including their original annotations from the Cellinker dataset.}
#' }
#'
#' @importFrom dplyr rename %>% filter select mutate arrange group_by summarise
#' @importFrom utils read.delim
#'
#' @export
#'
possible_L_R_pairs_Cellinker <- function(species,  gene_spot_expression_dataframe, percentage){
  if (species %in% c("Human","human") ){
    filepath <- system.file("extdata", "Cellinker_Human_interactions.txt", package = "SpaCCI")
  }else if (species %in% c("Mouse","mouse") ){
    filepath <- system.file("extdata", "Cellinker_Mouse_interactions.txt", package = "SpaCCI")
  }else{
    stop("Please enter species of either Human or Mouse")
  }

  db <- read.delim(filepath)
  db <- db[which(db$Ligand_location == "Secreted" & db$Receptor.location =="Membrane"),]
  split_and_trim <- function(x) { sapply(strsplit(x, ";", fixed = TRUE), function(y) trimws(y))}

  # Replace comma with plus
  ligand_string <- gsub(";", " + ", db$Ligand_symbol)
  receptor_string <- gsub(";", " + ", db$Receptor_symbol)
  inter <- paste0(ligand_string,"_",receptor_string)
  db$interaction_name <- inter
  rownames(db) <- db$interaction_name
  db <- db %>% rename(pathway_name = KEGG.pathway)
  # Apply the function to ligand and receptor columns
  db$ligand_vector <- lapply(db$Ligand_symbol, split_and_trim)
  db$receptor_vector <- lapply(db$Receptor_symbol, split_and_trim)
  db$combined_genes <- mapply(function(l, r) c(l, r),
                                   db$ligand_vector,
                                   db$receptor_vector,
                                   SIMPLIFY = FALSE)
  # function of checking if all genes of LR are in the data
  check_genes_presence <- function(gene_vector, gene_data) { all(gene_vector %in% rownames(gene_data))}
  exp_gene_express <- gene_spot_expression_dataframe[rowSums(gene_spot_expression_dataframe[] )>0,]
  binary_expression_data <- ifelse(exp_gene_express > 0, 1, 0)
  expression_percent <- apply(binary_expression_data, 1, function(row) sum(row > 0) / length(row) * 100)
  expressed_genes <- names(expression_percent)[expression_percent >= percentage ] # 1% 5% 10%
  final_gene_spot_expression_dataframe <- exp_gene_express[expressed_genes,]

  presence_vector <- mapply(check_genes_presence, db$combined_genes, MoreArgs = list(gene_data = final_gene_spot_expression_dataframe))
  # Subset inter to get only those L-R pairs where all genes are present in gene_data
  L_R_pairs <- db[presence_vector, c("ligand_vector", "receptor_vector","combined_genes" )]
  info_inter <- db[presence_vector, setdiff(names(db), c("ligand_vector", "receptor_vector", "combined_genes"))]
  return(list( possible_L_R_pairs = L_R_pairs, possible_L_R_pairs_details = info_inter))
}


#' Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' ICELLENT Database: This function identifies possible ligand-receptor (L-R) pairs for cell-cell communication analysis using data from the ICELLNET database. It checks for the presence of all genes involved in each L-R pair within the provided gene expression matrix, filtering based on a specified expression percentage threshold.
#'
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{possible_L_R_pairs}{A data frame of L-R pairs where all genes are present in the 'gene_spot_expression_dataframe' and meet the expression threshold. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{possible_L_R_pairs_details}{A data frame with detailed information about the identified L-R pairs, including their original annotations from the ICELLNET dataset.}
#' }
#'
#' @importFrom utils read.csv
#'
#' @export
#'
possible_L_R_pairs_ICELLNET <- function(gene_spot_expression_dataframe, percentage  ){
  filepath <- system.file("extdata", "ICELLNET_Human_interactions.csv", package = "SpaCCI")
  db <- read.csv(filepath, sep=";")

  db$Ligand <- apply(db[, grepl("^Ligand", names(db))], 1, function(row) {
    cleaned_row <- row[!is.na(row) & row != ""]
    paste(cleaned_row, collapse = ",")
  })
  db$Receptor <- apply(db[, grepl("^Receptor", names(db))], 1, function(row) {
    cleaned_row <- row[!is.na(row) & row != ""]
    paste(cleaned_row, collapse = ",")
  })
  db$pathway_name <- paste0(db$Family,"_",db$Subfamily)

  split_and_trim <- function(x) {
    sapply(strsplit(x, ",", fixed = TRUE), function(y) trimws(y))
  }
  # Replace comma with plus
  ligand_string <- gsub(",", " + ", db$Ligand)
  receptor_string <- gsub(",", " + ", db$Receptor)
  inter <- paste0(ligand_string,"_",receptor_string)
  db$interaction_name <- inter
  rownames(db) <- db$interaction_name
  # Apply the function to ligand and receptor columns
  db$ligand_vector <- lapply(db$Ligand, split_and_trim)
  db$receptor_vector <- lapply(db$Receptor, split_and_trim)
  db$combined_genes <- mapply(function(l, r) c(l, r),
                              db$ligand_vector,
                              db$receptor_vector,
                              SIMPLIFY = FALSE)
  # function of checking if all genes of LR are in the data
  check_genes_presence <- function(gene_vector, gene_data) { all(gene_vector %in% rownames(gene_data))}
  exp_gene_express <- gene_spot_expression_dataframe[rowSums(gene_spot_expression_dataframe[] )>0,]
  binary_expression_data <- ifelse(exp_gene_express > 0, 1, 0)
  expression_percent <- apply(binary_expression_data, 1, function(row) sum(row > 0) / length(row) * 100)
  expressed_genes <- names(expression_percent)[expression_percent >= percentage ] # 1, 5 ,10
  final_gene_spot_expression_dataframe <- exp_gene_express[expressed_genes,]

  presence_vector <- mapply(check_genes_presence, db$combined_genes, MoreArgs = list(gene_data = final_gene_spot_expression_dataframe))
  # Subset inter to get only those L-R pairs where all genes are present in gene_data
  L_R_pairs <- db[presence_vector, c("ligand_vector", "receptor_vector","combined_genes" )]
  info_inter <- db[presence_vector, setdiff(names(db), c("ligand_vector", "receptor_vector", "combined_genes"))]
  return(list( possible_L_R_pairs = L_R_pairs, possible_L_R_pairs_details = info_inter))
}


#' Identify Possible Ligand-Receptor Pairs for Cell-Cell Communication
#'
#' ConnectomeDB 2020 Database: This function identifies possible ligand-receptor (L-R) pairs based on gene expression data.
#'
#' @param gene_spot_expression_dataframe A gene expression data frame with genes as row names and Spot IDs as column names. This data frame is used to verify the presence of all genes involved in the L-R pairs.
#' @param percentage A numeric value specifying the minimum percentage of spots in which a gene must be expressed to be considered. The default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{possible_L_R_pairs}{A data frame of L-R pairs where all genes are present in the 'gene_spot_expression_dataframe' and meet the expression threshold. The data frame includes the ligand and receptor vectors, and the combined gene vectors.}
#'   \item{possible_L_R_pairs_details}{A data frame with detailed information about the identified L-R pairs, including their original annotations from the ConnectomeDB 2020 dataset.}
#' }
#'
#' @export

possible_L_R_pairs_connectome <- function(gene_spot_expression_dataframe, percentage ){
  lr_data_path <- system.file("extdata", "ConnectomeDB_2020_Human_interactions.RData", package = "SpaCCI")
  load(lr_data_path)
  db <- LR_pairs_ConnectomeDB_2020
  db$interaction_name <- paste0(db$ligand,"_",db$receptor)
  rownames(db) <- db$interaction_name
  db <- db[which(db$Ligand.location %in% c("secreted", "plasma membrane; secreted") ),]
  split_and_trim <- function(x) {
    sapply(strsplit(x, ",", fixed = TRUE), function(y) trimws(y))
  }
  db$pathway_name <- "Not_Annotated_in_this_database"
  # Apply the function to ligand and receptor columns
  db$ligand_vector <- lapply(db$ligand, split_and_trim)
  db$receptor_vector <- lapply(db$receptor, split_and_trim)
  db$combined_genes <- mapply(function(l, r) c(l, r),
                              db$ligand_vector,
                              db$receptor_vector,
                              SIMPLIFY = FALSE)
  # function of checking if all genes of LR are in the data
  check_genes_presence <- function(gene_vector, gene_data) { all(gene_vector %in% rownames(gene_data))}
  exp_gene_express <- gene_spot_expression_dataframe[rowSums(gene_spot_expression_dataframe[] )>0,]
  binary_expression_data <- ifelse(exp_gene_express > 0, 1, 0)
  expression_percent <- apply(binary_expression_data, 1, function(row) sum(row > 0) / length(row) * 100)
  expressed_genes <- names(expression_percent)[expression_percent >= percentage ] # 1, 5 ,10
  final_gene_spot_expression_dataframe <- exp_gene_express[expressed_genes,]

  presence_vector <- mapply(check_genes_presence, db$combined_genes, MoreArgs = list(gene_data = final_gene_spot_expression_dataframe ))
  # Subset inter to get only those L-R pairs where all genes are present in gene_data
  L_R_pairs <- db[presence_vector, c("ligand_vector", "receptor_vector","combined_genes" )]
  info_inter <- db[presence_vector, setdiff(names(db), c("ligand_vector", "receptor_vector", "combined_genes"))]
  return(list( possible_L_R_pairs = L_R_pairs, possible_L_R_pairs_details = info_inter))
}






