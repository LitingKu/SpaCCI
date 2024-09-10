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

Finally we run the local analysis on the whole slide with specifying
`specific_LR_pair`.

    # local analysis
    result_local <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
                               spot_cell_proportion_dataframe = cell_prop_df,
                               spatial_coordinates_dataframe = spatial_coords_df,
                               LR_database_list = LRDB,
                               specific_LR_pair = "EDN2_EDNRA",
                               analysis_scale = "local",
                               local_scale_proportion = 1,
                               neighborhood_radius = 2.5)

    ## [1] "Now analyzing localized detection using 100% of spots in the whole slide, with a radius of 2.5."

Then we plot the localized plot to access the local signalling hotspot.

    # if you have the Seurat Object with image 
    plot_SpaCCI_local(Seurat_Object = Example_Seurat,
                      SpaCCI_local_Result_List = result_local,
                      Ligand_cell_type = "ductal",
                      Receptor_cell_type = "activated_stellate",
                      spot_plot_size = 6)

    ## [1] "plotting using Seurat image"

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](SpaCCI_tutorial_files/figure-markdown_strict/unnamed-chunk-12-1.png)

    # if you DON'T have the Seurat Object with image, Please use the spatial_coords_df
    plot_SpaCCI_local(spatial_coordinates_dataframe = spatial_coords_df,
                      SpaCCI_local_Result_List = result_local,
                      Ligand_cell_type = "ductal",
                      Receptor_cell_type = "activated_stellate",
                      spot_plot_size = 6)

    ## [1] "plotting using image spatial coordinates"

![](SpaCCI_tutorial_files/figure-markdown_strict/unnamed-chunk-12-2.png)
