---
layout: splash
title: "Installation"
permalink: /installation/
subtitle: "Spatially Aware Cell-Cell Interaction Analysis"
header:
  overlay_color: "#333"
  overlay_image: cover.jpg  # Path to your cover image
  overlay_filter: 0.5  # Adjust the darkness of the image
---

<span style="font-family: 'ABeeZee', sans-serif; color: #94C47D;">SpaCCI</span> R package can be installed from GitHub using the following steps:

## Dependencies

- R version (>= 3.5.0)
- R packages: Rcpp (>= 1.0.13), RcppArmadillo, nnls, ggrepel, pheatmap, circlize (>= 0.4.12), Matrix, Seurat (>= 4.0.0), dplyr, patchwork, grDevices, reshape2, graphics, ggplot2


### 1. Install devtools

    if (!requireNamespace("devtools", quietly = TRUE)) {
         install.packages("devtools")
    }



### 2. Install <span style="font-family: 'ABeeZee', sans-serif; color: #94C47D;">SpaCCI</span> 

    devtools::install_github("LitingKu/SpaCCI")
    

### 3. Load <span style="font-family: 'ABeeZee', sans-serif; color: #94C47D;">SpaCCI</span>    

    library(SpaCCI)


The package has currently been test on:
- MAC: OS Sonoma (14.6.1)







