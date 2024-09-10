


    ## [1] "writing data frame"

    # plot the result heatmap
    plot_SpaCCI_heatmap(SpaCCI_Result_List = result_regional,
                        symmetrical = FALSE, cluster_cols = FALSE, return_tables = FALSE,
                        cluster_rows = FALSE, cellheight = 15, cellwidth = 15,
                        specific_celltypes = c(colnames(cell_prop_df)), alpha = 0.05,
                        main= "Significant Cell-Cell Interaction Count in the Interested Region")

![](SpaCCI_tutorial_files/figure-markdown_strict/unnamed-chunk-10-1.png)

#### (C) Local analysis

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
