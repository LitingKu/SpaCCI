---
subtitle: "Download Example Data"
layout: single
permalink: /datas/data_download/
toc: true
toc_sticky: true
sidebar:
  nav: datas
---


The following are the example dataset for the <span style="font-family: 'ABeeZee', sans-serif; color: #94C47D;">SpaCCI</span> tutorial.

## 1. Seurat Data for Tutorial
For the Seurat data used in the tutorial, please first install package and load the data as the following

    library(SpaCCI)
    # Load the tutorial data
    load(system.file("extdata", "Tutorial_example_data.rda", package = "SpaCCI"))

After loading Tutorial_example_data.rda, the following objects should be available:
- Example_Seurat: A Seurat object containing 10X Visium data with an associated image.
- cell_prop_df: A data frame containing spot-level cell type proportions.
- interest_region_Spot_IDs: Spot IDs for conducting "regional" analysis.
    
## 2. Non Seurat Data for Tutorial

- **Gene Expression Data**: A spot-level data frame with genes as row
    names and spot IDs as column names.
- **Cell Type Proportion Data**: A spot-level data frame with spot IDs
    as row names and cell type labels as column names.
- **Spatial Coordinates** :  

