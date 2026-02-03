##################################
#    MERFISH dataset analysis    #
##################################


# Note: Users should install stan, cmdstar first. 
#==================================================================
# Instructions about the tuning parameters: 
# 
# grid_size: the resolution of grids (default is 50*50)
# Basis: basis function approximation in HSGP (default is 25*25)  
#
# Framework of DenMark is:
# Input:  data of (1) cell locations and (2) gene expression 
# Model:  discretization with HSGP
# Output:  (1) posterior summary and (2) latent field decomposition 
#
# Overview of the approach:
# [Step 1]. Decide the study window. In thi example we take the whole brain tissue as the study window 
# [Step 2]. Input. Grid resolution  
# [Step 3]. DenMark with HSGP 
# [Step 4]. Output 
#==================================================================


# load the pkgs 
library('rstan')
library('MASS')
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
library(loo)
# for viz:
library('ggplot2')
library(cowplot)
library(viridis)
library(rethinking)
library(reshape2)
library(ggridges)
library(dplyr)
library(ggrepel)

# quantiles:
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


# for results showing usage 
options(cmdstanr_max_rows = 20)



# [Step 0]: load the spatial information matrix and gene expression matrix  ----------
# as well as transforming the whole coordinates systems to around origin 
# as required by HSGP approximation 
# 1. Translation 
# 2. scaling 
# 3. decide the hyperparameter to be used in the approximation 


# the spatial information & the gene expression information: 
load('data/cell-location.Rdata')
load('data/cell-gene-expression.Rdata')



# [Step 1]: Decide the window for analysis and the ROI normalization -----------
# short-range vs. large-range analysis 
# if short-range, decide the meaning window you are interested in 
# if large-range, decide the whole cell occurrence as window 
# Window format is [-L1, L1]*[L2, L2] and needs transformation of coordinates 

center_slice_cell_location<-  center_slice_cell_location/1000 # scale the micrometers into millimeters  
GeneIDls<- c(387,  403,  73, 435, 421 ) # These are the gene indeces of the five candidate genes shown in Figure 3. 
rownames(center_slice_gene_location) # check the gene name 




# User-defined hyperparameters: 
grid_size<- 50          #<- This is the grid size. 
basis_choice <- 25      #<- This is the number of basis functions in each of the xy directions    
boundary_choice <- 1.5  #<- This is the scaling factors 
















