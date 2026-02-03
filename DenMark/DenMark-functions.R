#################################
#     Functions in  DenMark     #
#################################

#****** This file contains following functions in developing DenMark: *************************************************************************
# 1. buildgridpp      : a function to discretize the single-cell resolution ST data (cell locations and gene expression per cell) by grids. 
# 2. grid_pts_centroid: a function to calculate and summaize the grid centroids, grid area, and the total number of grids 
# 3. DenMark          : a Bayesian framework to get the posterior 
#**********************************************************************************************************************************************
# Note: before running DenMark, you should have installed \texttt{stan} and \texttt{cmdstanr} for the Monte Carlo sampling, \texttt{fields} for generating the grids, \texttt{loo} for calculating WAIC  



#-------- 0. Load the Packages --------#

# Pkg for the algorithm framework 
library('rstan')
library(cmdstanr)

# Pkg for vizualization:
library('ggplot2')
library(cowplot)
library(viridis)
library(rethinking)
library(reshape2)
library(ggridges)
library(dplyr)
library(ggrepel)

# other Pkgs
library('MASS')
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(loo) 




#-------- 1. Function buildgridpp --------#

# Input: matrix with locations, gene expression at each point 
# Output: two vectors (Y) and (M); 2 matrics of Y and M for the analysis 
buildgridpp<- function(grid_size=grid_size, 
                       dataset=dataset){
  #========input===================
  # dataset: a dataframe with locations (x, y) and 1 gene expression as the colnames 
  # grid_size: grid resolution at each of the x- and y-axis
  #================================
  # the min and max locations
  colnames(dataset)<- c('x', 'y', 'gene')
  loc<- as.data.frame(cbind(dataset$x, dataset$y))
  colnames(loc)<- c('x', 'y')
  # sum up the grids at each grid 
  loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
  marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
  
  # initial values 
  init_x <- min(loc$x)
  init_y <- min(loc$y)
  
  # the delta: 
  step_x <- (max(loc$x)-min(loc$x))/grid_size
  step_y <- (max(loc$y)-min(loc$y))/grid_size
  
  # fill in the two matrix: 
  for (i in 1:grid_size){
    for (j in 1:grid_size){
      x_bound_1 <- init_x + (i-1) * step_x
      x_bound_2 <- init_x + i * step_x
      y_bound_1 <- init_y + (j-1) * step_y
      y_bound_2 <- init_y + j * step_y
      # index for each grid
      pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                       loc$x <= x_bound_2),
                               which(loc$y >= y_bound_1 & 
                                       loc$y <=  y_bound_2))
      # sum up the gene expression 
      marks_100[i,j]<- sum(dataset$gene[pts_index])
      loc_100[i,j] <- length(loc[pts_index,][,1])
    }
  }
  
  #======= return vectors ==== 
  # 1st element: 
  # the first column is the vectorized points 
  # the second column is the vectorized marks 
  # 2nd element: matrix for the points 
  # 3rd element: matrix for the marks 
  #===========================
  loc_mark_mat<- cbind(as.vector(loc_100), 
                       as.vector(marks_100))
  
  colnames(loc_mark_mat)<- c('point', 'mark')
  
  # summarize the results 
  loc_mark_res<- list(loc_mark_mat, 
                      loc_100, 
                      marks_100)
  
  return(loc_mark_res)
}



#---------- Function grid_pts_centroid ----------#

# Build another function for the locations of the centroids of those grids: 
# This is the function to sum up the characterisics we can get from the grids: 
# grid locations, grid_area, and the total number of grids 
grid_pts_centroid<- function(grid_size=grid_size, 
                             dataset=dataset){
  #========input===================
  # dataset: a dataframe with locations (x, y) and 1 gene expression as the colnames 
  # grid_size: grid resolution at each of the x- and y-axis
  #================================
  
  # the min and max locations
  loc<- as.data.frame(cbind(dataset$x, dataset$y))
  colnames(loc)<- c('x', 'y')
  # some key parameter for the bbox method
  max_coords_x <- max(loc$x)
  min_coords_x <- min(loc$x)
  max_coords_y <- max(loc$y)
  min_coords_y <- min(loc$y)
  
  # length of grid in two directions 
  len.grid.x <- ( max_coords_x - min_coords_x )/grid_size
  len.grid.y <- ( max_coords_y - min_coords_y )/grid_size
  
  # === The elements we care about ===  
  grid_area <- (len.grid.x * len.grid.y)
  
  # === The centroids of each grid === 
  grid_points <- expand.grid(x = seq(  min_coords_x-0.001,
                                       max_coords_x+0.001,
                                       length.out=grid_size), 
                             y = seq( min_coords_y-0.001,
                                      max_coords_y+0.001,
                                      length.out=grid_size))
  
  # === The number of total grid points 
  N <- nrow(grid_points)  # Total number of points
  
  
  #======= return list ==== 
  # the 1st element is the locations of the grids
  # the 2nd element is the grid area 
  # the 3rd element is the total number of the grid 
  #===========================
  grid_res_summary<- list(grid_points, 
                          grid_area, 
                          N)
  
  return(grid_res_summary)
}




#---------- quantile functions ----------#
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)



# 




#---------- DenMark framework ----------#
# input: 
# (1) grid_size
# (2) ROI region: first define and transform the ROI into a form [-L1, L1]*[-L2, L2], then provide the half-length in x and y as L1 and L2 
# (3) the gridded cell counts and gene expression per grid, calculated from the listed functions 
# (4) basis function at each direction (HSGP hyperparameters)
# (5) boundary factor (HSGP hyperparameters)
# 
# Expected output from the function: 
# (1) the sampled from the posterior distribution 
# (2) the Rhat statistics for convergenece check 
# (3) the WAIC criteria



DenMarkmain <- function(grid_size = grid_size, 
                        dataset = dataset, 
                        basis = basis, 
                        boundfactor=boundfactor){
# grid_size: total number of grids 
# basis:  number of basis function per direction (xy two dimensions)

  # packages
  require(cmdstanr)
  require(fields)
  require(posterior)
  require(loo)
 # ==================================================
  # Step 1: Grid geometry
  # ==================================================

  grid_info <- grid_pts_centroid(
    grid_size = grid_size,
    dataset   = dataset
  )

  coords      <- grid_info$coords
  grid_area   <- grid_info$grid_area
  log_grid_area <- log(grid_area)
  N           <- grid_info$N

  # Normalize coordinates
  L1 <- max(abs(coords[,1] - mean(coords[,1])))
  L2 <- max(abs(coords[,2] - mean(coords[,2])))

  coords[,1] <- coords[,1] - mean(coords[,1])
  coords[,2] <- coords[,2] - mean(coords[,2])

  # Distance-based prior scale
  distMat <- fields::rdist(coords)
  mrange  <- max(distMat) / (2 * 2.75)
# ==================================================
  # Step 2: Build gridded Y and M
  # ==================================================

  grid_data <- buildgridpp(
    grid_size = grid_size,
    dataset   = dataset
  )

  Y <- grid_data$vec[, "Y"]
  M <- grid_data$vec[, "M"]

  # ==================================================
  # Step 3: HSGP construction
  # ==================================================

  m1 <- basis
  m2 <- basis
  mstar <- m1 * m2

  Lstar <- c(L1, L2)
  L <- boundfactor * Lstar

  indices <- as.matrix(
    expand.grid(S2 = 1:m1, S1 = 1:m2)[, 2:1]
  )

  stan_data <- list(
    n = N,
    y1 = Y,
    y2 = M,
    d = 2,
    mstar = mstar,
    coords = coords,
    log_grid_area = log_grid_area,
    mrange = mrange,
    indices = indices,
    L = L,
    is_centerted_PHI = 0
  )

  # ==================================================
  # Step 4: Fit Stan
  # ==================================================

  init_fun <- function() list(
    Astar = matrix(c(1, 0, 0.5, 1), 2, 2),
    sigma = c(1, 1),
    beta0 = 0,
    beta1 = -1,
    ell   = c(1, 1),
    betab1 = rep(0, mstar),
    betab2 = rep(0, mstar)
  )

  mod <- cmdstan_model("../DenMark-M2.stan")

  fit <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    init = init_fun
  )

  # ==================================================
  # Step 5: Diagnostics & WAIC
  # ==================================================

  log_lik <- posterior::as_draws_matrix(
    fit$draws("log_lik")
  )

  list(
    fit = fit,
    rhat = fit$summary()$rhat,
    waic = loo::waic(log_lik),
    Ymat = grid_data$Ymat,
    Mmat = grid_data$Mmat
  )
  
} 











