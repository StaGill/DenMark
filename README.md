# DenMark
**DenMark** (<ins>**Den**</ins>sity-dependent <ins>**Mark**</ins>ed Point process framework) is a model-based statistical framework to quantify how gene expression varies with local cell density and to identify density-correlated genes (DCGs). It is designed for single-cell resolution spatial transcriptomics data such as MERFISH, Xenium and SeqFISH, where cell location and gene expression at the single-cell resolution is provided.

------
## DenMark Workflow 
![Overview the DenMark workflow](Images/denmark-workflow.png)

Implemented with a density-dependent marked point process, as well as comparing to the one with independent marked point process, DenMark enables downstream analyses such as:

- jointly quantify the spatial heterogeneity of the cell locations and a typical gene expression (candidate gene);
- quantify the correlation between cell density and gene expression
- identify the DCGs in the provided single-cell resolution spatial transcriptomics dataset 

This repository contains the reference R implementation used in the DenMark manuscript.


------
## Repository layout

- `DenMark/`  
  Core R implementation of DenMark (grid discretization and main modeling functions).

- `CodeInPaper/`  
  Scripts used to generate the figures and results in the manuscript
  (simulation study, MERFISH mouse brain data and Xenium breast cancer data).

- `Images/`  
  Images in the paper and this repo.

- `Tutorial_DenMark.rmd`  
  An RMarkdown tutorial that walks through two single-cell resolution datasets
  (MERFISH mouse brain and Xenium breast cancer data).

- `LICENSE`  
  License for using and modifying this code.


-----
## A Quick Start 

Please refer to the tutorial file Tutorial_DenMark.rmd. 

------
## Reproducing results from the manuscript

The scripts in CodeInPaper/ (to be documented) correspond to the main analyses:

- Simulation study  
Evaluates two approximation performance (grid-based approach vs. the actual marked point process; HSGP vs. exact GP).

-  MERFISH mouse brain data 
Quanfity the spatial heterogeneity in cell locations and gene expression, and candidate gene expression correlation to cell density. Identification of DCGs is also provided. 

- Xenium breast cancer data
Quanfity the spatial heterogeneity in cell locations and gene expression, and candidate gene expression correlation to cell density. Identification of DCGs is also provided. 

Data sources:

- Mouse brain MERFISH:  https://console.cloud.google.com/storage/browser/public-datasets-vizgen-merfish ;

- Human breast cancer Xenium: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast .




