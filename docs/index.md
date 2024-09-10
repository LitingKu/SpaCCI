---
layout: splash
homepage: true
title: "Welcome to SpaCCI"
header:
  overlay_color: "#333"
  overlay_image: cover.jpg  # Path to your cover image
  overlay_filter: 0.5  # Adjust the darkness of the image
  caption: "SpaCCI: Spatially Aware Cell-Cell Interaction Analysis"
---

# Welcome to SpaCCI

<p align="center">
  <img width="1000" src="SPACCI Method Diagram.png" alt="SpaCCI Method Diagram">
</p>

## SpaCCI key features
---
- **Integration of Multiple Ligand-Receptor Databases**: SpaCCI leverages widely used ligand-receptor interaction databases, including CellChat, CellPhoneDB, Cellinker, ICELLNET, and ConnectomeDB, for both human and mouse species. This ensures comprehensive coverage and accuracy in identifying potential cell-cell interactions.
- **Spatial Awareness**: Unlike traditional methods, SpaCCI incorporates spatial information from ST platforms like 10x Visium. This allows for more accurate detection of cell-cell interactions by considering the physical proximity of cells and the spatial distribution of cell types across tissue samples.
- **Cell Type Mixture Modeling**: SpaCCI acknowledges the mixture of cell types within each spot of the tissue sample. It uses non-negative least square regression combined with cell type proportion data to estimate cell type-specific gene expression, offering a more accurate quantification of interactions.
- **Global, Regional, and Local Interaction Quantification**: The framework is designed to analyze cell-cell interactions at local, regional, and global scales. By borrowing information from neighboring spots, SpaCCI provides a detailed view of interaction patterns within and across different regions of the tissue.

[about](about.md) 


