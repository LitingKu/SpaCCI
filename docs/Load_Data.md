---
title: "Load Data"
layout: single
permalink: /docs/Load_Data/
toc: true
toc_sticky: true
sidebar:
  nav: docs
---

This tutorial demonstrates how to load data and the input files that SpaCCI need.

### Required input data for SpaCCI analysis

SpaCCI requires two main input datasets:

1.  **Gene Expression Data**: A spot-level data frame with genes as row
    names and spot IDs as column names.
2.  **Cell Type Proportion Data**: A spot-level data frame with spot IDs
    as row names and cell type labels as column names.

With one correspond:

3. **Spatial Coordinates** : A spot-level data frame with spot IDs as row
   names and c("imagerow","imagecol") as column names


