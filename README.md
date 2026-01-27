# DenMark
DenMark (Density-dependent Marked Point process framework) is a model-based statistical framework to quantify how gene expression varies with local cell density and to identify density-correlated genes (DCGs). It is designed for single-cell resolution spatial transcriptomics data sich as MERFISH, Xenium and SeqFISH, where cell location and gene expression at the single-cell resolution is provided.

------
## DenMark Workflow 
![Overview the DenMark workflow](Images/denmark-workflow.png)

By aligning an “unknown” spectrum (or dataset) to a reference while forcing matched spectra to share the same m/z values, DenMark enables downstream analyses such as:

- jointly quantify the spatial heterogeneity of the cell locations and a typical gene expression (candidate gene);
- quantify the correlation between cell density and gene expression
- identify the DCGs in the provided single-cell resolution spatial transcriptomics dataset 

This repository contains the reference R implementation used in the DenMark manuscript.


------
## The files description 





------


