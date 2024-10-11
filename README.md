<p align="center">
  <img width="1000"  src="SPACCI Method Diagram.png">
</p>


# SpaCCI: Spatially Aware Cell-Cell Interaction Analysis

---

## SpaCCI key features

- **Integration of Multiple Ligand-Receptor Databases**: SpaCCI leverages widely used ligand-receptor interaction databases, including CellChat, CellPhoneDB, Cellinker, ICELLNET, and ConnectomeDB, for both human and mouse species. This ensures comprehensive coverage and accuracy in identifying potential cell-cell interactions.
- **Spatial Awareness**: Unlike traditional methods, SpaCCI incorporates spatial information from ST platforms like 10x Visium. This allows for more accurate detection of cell-cell interactions by considering the physical proximity of cells and the spatial distribution of cell types across tissue samples.
- **Cell Type Mixture Modeling**: SpaCCI acknowledges the mixture of cell types within each spot of the tissue sample. It uses non-negative least square regression combined with cell type proportion data to estimate cell type-specific gene expression, offering a more accurate quantification of interactions.
- **Global, Regional and Local Interaction Quantification**: The framework is designed to analyze cell-cell interactions at among local, regional and global scales. By borrowing information from neighboring spots, SpaCCI provides a detailed view of interaction patterns within and across different regions of the tissue.

---

## Installation

SpaCCI R package can be easily installed from 
1. r cran:

```
install.packages("SpaCCI")
```

2. Github using devtools:

```
install.packages("devtools")
devtools::install_github("LitingKu/SpaCCI")
```
---

## SpaCCI Tutorials

- [Simple Tutorials](https://litingku.github.io/SpaCCI/docs/tutorial1/)
---

## Help, Suggestion and Contribution

---


