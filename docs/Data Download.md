---
title: "Download Example Data"
layout: single
permalink: /datas/data_download/
toc: true
toc_sticky: true
sidebar:
  nav: datas
---


The following are the example dataset for the <span style="font-family: 'ABeeZee', sans-serif; color: #94C47D;">SpaCCI</span> tutorial.

## 1. Seurat Data for Tutorial
For the Seurat object used in the tutorial, please first install package and load the data as the following

    library(SpaCCI)
    # Load the tutorial data
    load(system.file("extdata", "Tutorial_example_data.rda", package = "SpaCCI"))

After loading `Tutorial_example_data.rda`, the following objects should be available:
- `Example_Seurat`: A Seurat object containing 10X Visium data with an associated image.
- `cell_prop_df`: A data frame containing spot-level cell type proportions.
- `interest_region_Spot_IDs`: Spot IDs for conducting "regional" analysis.
    
## 2. Non Seurat Data for Tutorial
If NOT using the Seurat object, please also first download the following data file. 

- **<a href="https://github.com/LitingKu/SpaCCI/blob/main/SpaCCI_tutorial_files/data/normalized_gene_spot_df.csv" style="color: #38a3a5;"> Gene Expression Data</a>**: A spot-level data frame with genes as row
    names and spot IDs as column names.
- **<a href="https://github.com/LitingKu/SpaCCI/blob/main/SpaCCI_tutorial_files/data/cell_prop_df.csv" style="color: #38a3a5;">Cell Type Proportion Data</a>**: A spot-level data frame with spot IDs
    as row names and cell type labels as column names.
- **<a href="https://github.com/LitingKu/SpaCCI/blob/main/SpaCCI_tutorial_files/data/spatial_coords_df.csv" style="color: #38a3a5;">Spatial Coordinates</a>** : A spot-level data frame with spot IDs as row names and `c("imagerow","imagecol")` as column names.

For regional analysis, please also load the data as following to obtain the `interest_region_Spot_IDs`:

    library(SpaCCI)
    # Load the tutorial data
    load(system.file("extdata", "Tutorial_example_data.rda", package = "SpaCCI"))


## 3. Real 10X Visium Data

<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233254"> GSE233254</a>: Contains the 10X Visium data.

<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233293"> GSE233293</a>: Contains the RCTD(Robust Cell Type Deconvolution) of the 10X Visium data.

from <a href="https://doi.org/10.1158/2159-8290.CD-22-1200"> Sans M, Makino Y, Min J, et al. Spatial Transcriptomics of Intraductal Papillary Mucinous Neoplasms of the Pancreas Identifies NKX6-2 as a Driver of Gastric Differentiation and Indolent Biological Potential. Cancer Discov. 2023;13(8):1844-1861. doi:10.1158/2159-8290.CD-22-1200</a>




    

