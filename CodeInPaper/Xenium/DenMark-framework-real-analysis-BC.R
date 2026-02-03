###############################
# DenMark general framework 
# to be used in the real analysis for MERFISH and BC cells data 
# 
###############################

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
# [Step 1]. Decide the study window. Short-range or large-range?
# [Step 2]. Input. Grid resolution and  
# [Step 3]. DenMark with HSGP 
# [Step 4]. Output 
#==================================================================

# load packages: 
library('ggplot2')
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
library(cowplot)
library(viridis)
library(rethinking)
library(reshape2)
library(posterior)

# quantiles:
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


# for results showing usage 
options(cmdstanr_max_rows = 20)

# [Step 0]: load the spatial information matrix and gene expression matrix, ------
# as well as transforming the whole coordinates systems to around origin 
# as required by HSGP approximation 
# 1. Translation 
# 2. scaling 
# 3. decide the hyperparameter to be used in the approximation 

# the spatial information & the gene expression information: 
cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/cell_spatial.csv')
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000

# check the genes in the dataset:
colnames(cell_gene)



# 0C. hyperparameters specification -----
# grid_size
# hypyparameters in HSGP
grid_size<- 50
basis_choice <- 25
boundary_choice <- 1.5


# QC: filter out the unrelated genes: -----







# [Step 1]: Decide the window for analysis-----------
# short-range vs. large-range analysis 
# if short-range, decide the meaning window you are interested in 
# if large-range, decide the whole cell occurrence as window 
# Window format is [-L1, L1]*[L2, L2] and needs transformation of coordinates 


## 1A. decide to use one cell type/ one region of interest/... ----

# Please Select from: 
# Choice 1 (large range, but restrict on 1 cell type)
# Choice 2 (short range, restrict on one region)


### Choice 1: extract the  proliferative invasive tumor cell type ---
cell_spatial_tumor<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
plot(cell_spatial_tumor$x, 
     cell_spatial_tumor$y)
gene_spatial_tumor<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]

# df is a data frame containing x-coords, y-coords, and gene expression 
df<- as.data.frame(cbind( cell_spatial_tumor$x, 
                          cell_spatial_tumor$y ,
                          gene_spatial_tumor$CCND1 ))
colnames(df)<- c('x','y','Gene')

# center_slice_cell_location is a data frame containing spatial locations 
center_slice_cell_location<- as.data.frame(cbind(cell_spatial_tumor$x, 
                                                 cell_spatial_tumor$y ))

colnames(center_slice_cell_location)<- c('x','y')

# summarize: 
loc<- center_slice_cell_location
eg_gene_express<- df$Gene

### Choice 2: extract a small window containing all cells (short-range analysis) ---


#cell_spatial$ERBB2<- cell_gene$ERBB2
#cell_spatial$ESR1<- cell_gene$ESR1
#cell_spatial$PGR<- cell_gene$PGR


## 1B. decide the window of interest -----

center_slice_cell_location_DCIS_I<- cell_spatial[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5 & cell_spatial$y>= 1.8 & cell_spatial$y <= 2.8),]
center_slice_cell_location_DCIS_II<- cell_spatial[which(cell_spatial$x>=5 & cell_spatial$x<= 6 & cell_spatial$y>= 0.75 & cell_spatial$y <= 1.75),]
center_slice_cell_location_invasive<- cell_spatial[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y <= 3),]


gene_dcis1<- cell_gene[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5 & cell_spatial$y>= 1.8 & cell_spatial$y <= 2.8),]
gene_dcis2<- cell_gene[which(cell_spatial$x>=5 & cell_spatial$x<= 6 & cell_spatial$y>= 0.75 & cell_spatial$y <= 1.75),]
gene_invasive<- cell_gene[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y <= 3),]



# Quality control: 
# filter out the genes with (1) the lowest 5% , (2) zeros with more than 5% locations
gene_invasive<-gene_invasive[,-1]

# Convert to matrix (if not already)
mat <- as.matrix(gene_invasive)

# -------------------------------
# (1) Compute variance per gene
# -------------------------------
gene_var <- apply(mat, 2, var)

# Keep top 95% high-variance genes
var_threshold <- quantile(gene_var, 0.05)
genes_high_var <- names(gene_var[gene_var >= var_threshold])

# -------------------------------
# (2) Non-zero filter (>=5%)
# -------------------------------
nonzero_fraction <- colSums(mat > 0) / nrow(mat)
genes_nonzero <- names(nonzero_fraction[nonzero_fraction >= 0.05])

# -------------------------------
# Combine both filters
# -------------------------------
#genes_keep <- intersect(genes_high_var, genes_nonzero)

# Filtered dataset
#gene_invasive_filtered <- mat[, genes_keep, drop = FALSE]

# save the keeped gene names:
#write.csv(genes_keep, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-I-after-QC.csv')
#write.csv(gene_dcis1_filtered, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-I-mat.csv')
#write.csv(gene_dcis2_filtered, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-II-mat.csv')
#write.csv(gene_invasive_filtered, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/Invasive-mat.csv')




# for the large range analysis: ----

center_slice_cell_location_DCIS<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
center_slice_cell_location_tumor<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
center_slice_cell_location_immune<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),]


gene_dcis<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
gene_tumor<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]




# Quality control: 
# filter out the genes with (1) the lowest 5% , (2) zeros with more than 5% locations
gene_dcis<-gene_dcis[,-1]

# Convert to matrix (if not already)
mat <- as.matrix(gene_dcis)

# -------------------------------
# (1) Compute variance per gene
# -------------------------------
gene_var <- apply(mat, 2, var)

# Keep top 95% high-variance genes
var_threshold <- quantile(gene_var, 0.05)
genes_high_var <- names(gene_var[gene_var >= var_threshold])

# -------------------------------
# (2) Non-zero filter (>=5%)
# -------------------------------
nonzero_fraction <- colSums(mat > 0) / nrow(mat)
genes_nonzero <- names(nonzero_fraction[nonzero_fraction >= 0.05])

# -------------------------------
# Combine both filters
# -------------------------------
genes_keep <- intersect(genes_high_var, genes_nonzero)

# Filtered dataset
gene_dcis_filtered <- mat[, genes_keep, drop = FALSE]

write.csv(genes_keep, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/DCIS-after-QC.csv')
write.csv(gene_dcis_filtered, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/DCIS-mat.csv')




# Then we can calculate the window in the long-range analysis ----

#center_slice_cell_location<- cell_spatial[which(cell_spatial$x>=4 & cell_spatial$x<= 4.5 & cell_spatial$y>= 4.7 & cell_spatial$y <= 5.2),]

# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))

# the centroids of each grid: 
coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))

# distance matrix:
distMat <- fields::rdist(coords)
# the empirical range, used for prior specification 
mrange <- max(distMat)/(2*2.75)
# the area of the grid, used for adjustment of offset 
grid_area<- (2*L1*2*L2)/(grid_size^2)
log_grid_area <- log(grid_area) 




# [Step 2]: Preprocessing for the input Y and M ------


## 2A. Filter out uninformative genes-----
# (just for whole gene analysis usage)



## 2B. Combine cell-level N and M into a data frame. (colnames: x-coordinates, y-coordinates, and gene expression)-----




# To sum up for the 
#length(loc$x)==length(eg_gene_express) # check the length are the same 

## 2C. Create grid-level N and M for the gridded model ---- 

# create grids:
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of grids


# the location 
loc<- center_slice_cell_location
eg_gene_express<- center_slice_cell_location$PGR

# prepare the grids: 
init_x <- min(loc$x)
max(loc$x)
init_y <- min(loc$y)
max(loc$y)

# the delta: 
step_x <- (max(loc$x)-min(loc$x))/grid_size
step_y <- (max(loc$y)-min(loc$y))/grid_size
#grid_area<- step_x*step_y

# sum up the grids
loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
for (i in 1:grid_size){
  for (j in 1:grid_size){
    x_bound_1 <- init_x + (i-1) * step_x
    x_bound_2 <- init_x + i * step_x
    y_bound_1 <- init_y + (j-1) * step_y
    y_bound_2 <- init_y + j * step_y
    
    pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                     loc$x <= x_bound_2),
                             which(loc$y >= y_bound_1 & 
                                     loc$y <=  y_bound_2))
    marks_100[i,j]<- sum(eg_gene_express[pts_index])
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}

# To sum up for the gridded data:
Y<- as.vector(loc_100)
M<- as.vector(marks_100)


# Viz for the gridded data (can be hidden): 

# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# ---- Prepare data frames for ggplot ----
df_loc <- melt(loc_mat)
df_marks <- melt(marks_mat)
# Ensure same shape before dividing
stopifnot(all(dim(loc_mat) == dim(marks_mat)))

# ---- Compute expression per cell ----
expr_per_cell_mat <- marks_mat / (loc_mat+0.05)
df_expr <- melt(expr_per_cell_mat)


# But we are only interested in one gene to do analysis at a time



# [Step 3]: DenMark framework with HSGP -----


## 3A. The prepration for the HSGP-----
################################################
# HSGP approximation to general model -----
################################################

xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# 22*22
m1 <- basis_choice; m2 <- basis_choice; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(boundary_choice, boundary_choice)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)


## 3B. the input and the method -----

##########################################################
## Order the input
##########################################################

# adding the grid_area term,
input <- list(n = N,
              y1 = Y, # outcome1: points 
              y2 = M, # outcome2: marks 
              d = 2, # dim=2
              mstar = mstar, 
              coords = coords,
              log_grid_area=log_grid_area,
              mrange=mrange,
              indices = S,
              L=L,
              is_centerted_PHI = 0)

# define initial values 
init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 1), 2, 2),
  sigma=c(1, 1),
  beta0=0, 
  beta1=-1, 
  ell=c(1, 1),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)


str(input)

stan_file <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/M2-latent-field-mod.stan'
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_HS <- mod$sample(data = input, 
                             chains = 4,
                             parallel_chains = 4,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             adapt_delta = 0.95,
                             init=init_fun, 
                             max_treedepth = 10,
                             step_size = 1)
elapsed_time <- cmdstan_fit_HS$time()
elapsed_time
elapsed_time$total/60


q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


fit_summary_HS <- cmdstan_fit_HS$summary(variables = c("beta0", "beta1",   "ell[1]", 'ell[2]', 
                                                       'Astar[2,1]', 'sigma[1]','sigma[2]'), 
                                         c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS


# save the draws for the results: 
# Extract draws for log_lik
log_lik_array <- cmdstan_fit_HS$draws("log_lik") # 3D array: iterations × chains × n

# Convert to matrix for loo
log_lik_mat <- posterior::as_draws_matrix(log_lik_array) # iterations × n

# save all draws from the model: 
mod_draws<- cmdstan_fit_HS$draws()

# save the draws: 
#saveRDS(mod_draws, file = "C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/draws-for-BC-data/cmdstan_fit_draws_ERBB2.rds")
#saveRDS(mod_draws, file = "C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/draws-for-BC-data/cmdstan_fit_draws_ESR1.rds")
#saveRDS(mod_draws, file = "C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/draws-for-BC-data/cmdstan_fit_draws_PGR.rds")

cmdstan_fit_HS <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/dataset-analysis-genes/analysis_Aqp4.rds')



# results: 
#> fit_summary_HS
## A tibble: 7 × 11
#variable     mean     sd  `2.5%`   `25%`  `50%`  `75%` `97.5%`  rhat ess_bulk
#<chr>       <dbl>  <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl> <dbl>    <dbl>
#  1 beta0       5.63  0.665   4.05    5.26    5.73   6.09    6.68   1.00    1283.
#2 beta1      -0.472 0.522  -1.64   -0.781  -0.416 -0.116   0.405  1.00    1321.
#3 ell[1]      0.113 0.0278  0.0687  0.0931  0.110  0.130   0.178  1.00    1167.
#4 ell[2]      0.408 0.855   0.0171  0.0465  0.143  0.506   2.02   1.00     754.
#5 Astar[2,1]  0.915 0.106   0.601   0.892   0.958  0.983   0.998  1.00    2023.
#6 sigma[1]    1.40  0.348   0.858   1.15    1.35   1.60    2.18   1.00    1362.
#7 sigma[2]    1.19  0.325   0.714   0.960   1.14   1.37    1.97   1.00    1371.




# [Step 4]: Output Analysis-----
# (1) Posterior summary  
# (2) Latent field decomposition


# the latent field decomposition: 
w1_summary <- cmdstan_fit_HS$summary(variables = paste0("w1[",1:grid_size^2,"]"), 
                                     c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))

w2_summary <- cmdstan_fit_HS$summary(variables = paste0("w2[",1:grid_size^2,"]"), 
                                     c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))


# summary of the posterior 
samples_beta0_HS<- draws_df_HS[['beta0']]
samples_beta1_HS<- draws_df_HS[['beta1']]
samples_corr_HS<- draws_df_HS[['Astar[2,1]']]
samples_rho1_HS<- draws_df_HS[['ell[1]']]
samples_rho2_HS<- draws_df_HS[['ell[2]']]
samples_sigma1_HS<- draws_df_HS[['sigma[1]']]
samples_sigma2_HS<- draws_df_HS[['sigma[2]']]
# A matrix 
samples_a11_HS<- samples_sigma1_HS
samples_a21_HS<- samples_sigma2_HS*samples_corr_HS
samples_a22_HS<- samples_sigma2_HS*sqrt(1-samples_corr_HS^2)



# a11*w1= sigma1*w1
cal_a11w1<- function(i){
  w1<-cmdstan_fit_HS$draws(variables = paste0("w1[",i, "]"))
  samples_a11_HS*w1
}

# a21*w1 = a21_star*sigma2*w1
cal_a21w1<- function(i){
  w1<-cmdstan_fit_HS$draws(variables = paste0("w1[",i, "]"))
  samples_a21_HS*w1
}

# a22*w2= sqrt(1-a21.star^2)*sigma2 * w2
cal_a22w2<- function(i){
  w2<-cmdstan_fit_HS$draws(variables = paste0("w2[",i, "]"))
  samples_a22_HS*w2
}

# eta1 = a11*w1
cal_eta1<- cal_a11w1


# eta2 = a21*w1+a22*w2
cal_eta2<- function(i){
  w1<- cmdstan_fit_HS$draws(variables = paste0("w1[",i, "]"))
  w2<- cmdstan_fit_HS$draws(variables = paste0("w2[",i, "]"))
  (samples_a21_HS*w1) + (samples_a22_HS*w2)
}


# w1 
cal_w1<- function(i){
  cmdstan_fit_HS$draws(variables = paste0("w1[",i, "]"))
}
samples_w1<- sapply(1:(grid_size^2), cal_w1)

# w2 
cal_w2<- function(i){
  cmdstan_fit_HS$draws(variables = paste0("w2[",i, "]"))
}
samples_w2<- sapply(1:(grid_size^2), cal_w2)


# samples for the decomposition 
samples_a11_w1<- sapply(1:(grid_size^2), cal_a11w1)
samples_a21_w1<- sapply(1:(grid_size^2), cal_a21w1)
samples_a22_w2<- sapply(1:(grid_size^2), cal_a22w2)
samples_eta1<- samples_a11_w1
samples_eta2<- sapply(1:(grid_size^2), cal_eta2)



# add legend bar to each of the plot: 
library(viridis)
library(fields)

# --- helper for rotation ---
rotate_reflect <- function(m) {
  t(apply(m, 2, rev))
}

# --- compute matrices ---
mat1 <- rotate_reflect(apply(matrix(apply(samples_beta0_HS + samples_a11_w1, 2, mean),
                                    nrow = grid_size, byrow = TRUE), 2, rev))
mat2 <- rotate_reflect(apply(matrix(apply(samples_a11_w1, 2, mean),
                                    nrow = grid_size, byrow = TRUE), 2, rev))
mat3 <- rotate_reflect(apply(matrix(apply(samples_beta1_HS + samples_eta2, 2, mean),
                                    nrow = grid_size, byrow = TRUE), 2, rev))
mat4 <- rotate_reflect(apply(matrix(apply(samples_a21_w1, 2, mean),
                                    nrow = grid_size, byrow = TRUE), 2, rev))
mat5 <- rotate_reflect(apply(matrix(apply(samples_a22_w2, 2, mean),
                                    nrow = grid_size, byrow = TRUE), 2, rev))

# --- define row-specific scales ---
vals_row1 <- c(mat1, mat2)
breaks_row1 <- seq(min(vals_row1), max(vals_row1), length.out = 16)
cols_row1   <- viridis(length(breaks_row1) - 1, option = "magma")

vals_row2 <- c(mat3,mat4, mat5)
breaks_row2 <- seq(min(vals_row2), max(vals_row2), length.out = 16)
cols_row2   <- viridis(length(breaks_row2) - 1, option = "magma")



# the plot for the lambda1 and lambda2:
layout(matrix(c(1,2,3,4, 5,6,7,8), nrow=2, byrow=TRUE))
par(mar = c(3, 3, 2, 1), oma = c(0, 0, 0, 0))

# Row 1 (shared scale for row 1)
image.plot(exp(mat1), col = cols_row1, breaks = seq(min(exp(mat1)),
                                                    max(exp(mat1)),
                                                    length.out=16),
           main = expression(lambda[1](g)),
           xlab = '', ylab = '', axes = FALSE,
           legend.mar = 4)

image.plot(mat1, col = cols_row1, breaks = breaks_row1,
           main = expression(log(lambda[1](g))),
           xlab = '', ylab = '', axes = FALSE,
           legend.mar = 4)

image.plot(mat2, col = cols_row1, breaks = breaks_row1,
           main = expression(a[11]*w[1](g)),
           xlab = '', ylab = '', axes = FALSE,
           legend.mar = 4)

plot.new()  # empty slot

# Row 2 (shared scale for row 2)
image.plot(exp(mat3), col = cols_row2, breaks = seq(min(exp(mat3)),
                                                    max(exp(mat3)),
                                                    length.out=16),
           main = expression(lambda[2](g)),
           xlab = '', ylab = '', axes = FALSE,  legend.mar = 4)

image.plot(mat3, col = cols_row2, breaks = breaks_row2,
           main = expression(log(lambda[2](g))),
           xlab = '', ylab = '', axes = FALSE,  legend.mar = 4)

image.plot(mat4, col = cols_row2, breaks = breaks_row2,
           main = expression(a[21]*w[1](g)),
           xlab = '', ylab = '', axes = FALSE,
           legend.mar = 4)

image.plot(mat5, col = cols_row2, breaks = breaks_row2,
           main = expression(a[22]*w[2](g)),
           xlab = '', ylab = '', axes = FALSE,
           legend.mar = 4)







# the correlation Corr(Y, M/Y):

lambda1_g <- vector("list", grid_size^2)
lambda2_g <- vector("list", grid_size^2)
sim.Y <- vector("list", grid_size^2)
sim.M <- vector("list", grid_size^2)
sim.Z <- vector("list", grid_size^2)

for (g in 1:grid_size^2){
  w1_col <- paste0("w1[", g, "]")
  w2_col <- paste0("w2[", g, "]")
  # A matrix: 
  a11 <- draws_df_HS$`sigma[1]`
  a21 <- draws_df_HS$`sigma[2]`*draws_df_HS$`Astar[2,1]`
  a22<-  sqrt((draws_df_HS$`sigma[2]`)^2-a21^2)
  # log(lambda1(g)) = beta0 + a11*w1(g)
  lambda1_g[[g]] <- draws_df_HS$beta0 + a11 * draws_df_HS[[w1_col]]
  # simulate the points: 
  sim.Y[[g]] <- rpois(length(lambda1_g[[g]]),lambda=exp(log_grid_area +lambda1_g[[g]]))
  # log(lambda2(g)) = beta1 + log(sim.Y1 + 0.05) + a21*w1 + a22*w2
  lambda2_g[[g]] <- draws_df_HS$beta1 + log(sim.Y[[g]] + 0.05) + a21 * draws_df_HS[[w1_col]] + a22 * draws_df_HS[[w2_col]] 
  # simulate the marks:
  sim.M[[g]] <- rpois(length(lambda1_g[[g]]), lambda=exp(lambda2_g[[g]]))
  # calculate the scaled gene expression: 
  sim.Z[[g]] <- sim.M[[g]]/(sim.Y[[g]]+0.05)
}

# convert the Y and Z into matrix 
Y_mat <- do.call(cbind, sim.Y)  # [n_draws x grid_size^2]
Z_mat <- do.call(cbind, sim.Z)  # same

# the correlation 
corr_dist_1 <- sapply(1:nrow(Y_mat), function(i) {
  cor(Y_mat[i, ], Z_mat[i, ])
})

par(mfrow=c(1,1))
hist(corr_dist_1)

quantile(corr_dist_1, 0.025)
quantile(corr_dist_1, 0.975)


plot(density(corr_dist_1), col='red', lwd=2, xlab='Correlation', ylab='Density',main='' )
plot(density(corr_dist_1), col='red', lwd=2, xlab='Correlation', ylab='Density',main='' )
abline(v= cor(Y, M/(Y+0.05)), col='black', lwd=5)



par(mar=c(5,5,4,2), mfrow=c(1,1))
dens <- density(corr_dist_1)
plot(dens, col='red', lwd=5, 
     xlab='Correlation', ylab='Density', main='',
     ylim=c(0, max(dens$y) * 1.2), bty="l",
     xlim=c(0,0.5),
     cex.lab=2,
     cex.axis=2)  # "l" removes top/right box


cor(Y, M/(Y+0.05))






dens <- density(samples_corr_HS)
par(mar=c(5,5,4,2), mfrow=c(1,1))
# Extend y-axis so the CI segment has space
plot(dens, col='red', lwd=5, 
     xlab=expression('Correlation'~rho), ylab='Density', main='',
     ylim=c(0, max(dens$y) * 1.2), bty="l",
     xlim=c(0,1), # add the x-axis limit for all genes 
     cex.lab=2,
     cex.axis=2)  # "l" removes top/right box

# Add prior line
lines(density(rlkjcorr(100000, 2, eta=1.2)[,1,2]), col='blue', 
      lty=2, lwd=5)

# Compute CI bounds
ci_lower <- quantile(samples_corr_HS, 0.025)
ci_upper <- quantile(samples_corr_HS, 0.975)
med_ci<- median(samples_corr_HS)

# Add horizontal line segment at the top
y_max <- max(dens$y)
y_seg <- y_max * 1.1
segments(ci_lower, y_seg, ci_upper, y_seg, col='black', lwd=8)

# Add dots at ends of the segment
points(c(ci_lower, ci_upper), c(y_seg, y_seg), pch=16, cex=2)





### the correlation decay function -----
distMat <- fields::rdist(coords)
max_dist<- max(distMat)

# 1.2 Define Matérn covariance function for nu = 3/2
matern32 <- function(d, sigma, rho) {
  sigma^2*(1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}



# corr decay plot -----
dist_choice<- seq(0, max_dist, by=0.1)

samples_rho1<- summary_gene_22$`ell[1]`
samples_rho2<- summary_gene_22$`ell[2]`
samples_a21_HS<- summary_gene_22$`sigma[2]`*summary_gene_22$`Astar[2,1]`
samples_a22_HS<- summary_gene_22$`sigma[2]`*sqrt(1- summary_gene_22$`Astar[2,1]`^2)


post_rho1<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1))
post_rho2<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho2))

post_mean_rho1<- rep(NA, length(dist_choice))
post_upper_rho1<- rep(NA, length(dist_choice))
post_lower_rho1<- rep(NA, length(dist_choice))

post_mean_rho2<- rep(NA, length(dist_choice))
post_upper_rho2<- rep(NA, length(dist_choice))
post_lower_rho2<- rep(NA, length(dist_choice))


for(k in 1:length(dist_choice)){
  post_rho1[k,]<- matern32(d=dist_choice[k], sigma=1, rho= samples_rho1)
  post_rho2[k,]<- (samples_a21_HS^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1)+
                     samples_a22_HS^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2))/(samples_a21_HS^2+samples_a22_HS^2)
}


for (k in 1:length(dist_choice)){
  post_mean_rho1[k]<-mean(post_rho1[k,])
  post_upper_rho1[k]<-quantile(post_rho1[k,], 0.975)
  post_lower_rho1[k]<-quantile(post_rho1[k,], 0.025)
  
  post_mean_rho2[k]<-mean(post_rho2[k,])
  post_upper_rho2[k]<-quantile(post_rho2[k,], 0.975)
  post_lower_rho2[k]<-quantile(post_rho2[k,], 0.025)
}




df_rho1<- cbind( dist_choice ,post_mean_rho1, post_upper_rho1, post_lower_rho1)
colnames(df_rho1)<- c('dist_choice',  'post_mean_rho1', 'post_upper_rho1', 'post_lower_rho1')
est_range_1<-dist_choice[which(abs(df_rho1[,2]-0.05)==min(abs(df_rho1[,2]-0.05)))]


eta1_corr<-ggplot(data=df_rho1, aes(x=dist_choice))+
  geom_ribbon(aes(ymin=post_lower_rho1, ymax= post_upper_rho1), fill = "grey70")+
  geom_line(aes(y=post_mean_rho1), color='black', linewidth=1.05)+
  geom_hline(yintercept=0.05, linetype='dashed', col='#ED1B2F', linewidth=1.5)+
  labs(x='Distance (Millimeter)', y=expression( 'Posterior distribution of '~ r[1]) )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -0, hjust = 0))+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = "none")+
  annotate('text', 
           x = max_dist * 0.7,  
           y = 0.2,                    
           label = paste0("Practical Range = ",est_range_1 ), 
           color = "#ED1B2F", 
           size = 8)


df_rho2<- cbind( dist_choice ,post_mean_rho2, post_upper_rho2, post_lower_rho2)
colnames(df_rho2)<- c('dist_choice',  'post_mean_rho2', 'post_upper_rho2', 'post_lower_rho2')
est_range_2<-dist_choice[which(abs(df_rho2[,2]-0.05)==min(abs(df_rho2[,2]-0.05)))]


eta2_corr<-ggplot(data=df_rho2, aes(x=dist_choice))+
  geom_ribbon(aes(ymin=post_lower_rho2, ymax= post_upper_rho2), fill = "grey70")+
  geom_line(aes(y=post_mean_rho2), color='black', linewidth=1.05)+
  geom_hline(yintercept=0.05, linetype='dashed', col='#ED1B2F', linewidth=1.5)+
  labs(x='Distance (Millimeter)', y=expression( 'Posterior distribution of '~ r[2]))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -0, hjust = 0))+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = "none"
  )+
  annotate('text', 
           x = max_dist * 0.7,  
           y = 0.2,                     
           label = paste0("Practical Range = ",est_range_2 ), 
           color = "#ED1B2F", 
           size = 8)

library(cowplot)
plot_grid(eta1_corr, eta2_corr)

est_range_1
est_range_2




# II. the point and gene expression plot ----
p0<- ggplot(cell_spatial_tumor, aes(x=x, y=y))+
  geom_point(size=0.9)+
  labs(title='',x='',y='')+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )+
  theme_void()
#theme_classic(base_size = 14)


df<- as.data.frame(cbind(cell_spatial_tumor$x, 
                         cell_spatial_tumor$y,
                         gene_spatial_tumor))

colnames(df)[1]<- 'x'
colnames(df)[2]<- 'y'


p1<- ggplot(df, aes(x = x, y = y, color = CCND1)) +
  geom_point(size = 0.5, alpha = 1) +
  scale_color_gradient(low = "snow1", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme_void()+
  theme(
    text = element_text(size = 40),
    legend.title = element_text(size = 50, margin = margin(b = 30)), # more space
    legend.text  = element_text(size = 40),
    legend.margin = margin(t = 10, b = 10),   # add breathing room
    legend.spacing.y = unit(0.5, "cm"),       # vertical spacing between title and labels
    strip.text.x = element_text(size = 40),
    strip.text.y = element_text(size = 40)
  )+
  guides(color = guide_colorbar(
    barheight = unit(10, "cm"),  # adjust vertical length
    barwidth  = unit(1.5, "cm")
  ))
#theme_classic(base_size = 14)



# Girdded area: 

# ---- Consistent color breaks ----
breaks <- seq(
  min(c(df_loc$value, df_expr$value), na.rm = TRUE),
  max(c(df_loc$value, df_expr$value), na.rm = TRUE),
  length.out = 20
)

# ---- Define base theme ----
base_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Plot 1: Gridded Cell Counts ----
p1 <- ggplot(df_loc, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(
    option = "magma",
    breaks = breaks,
    labels = function(x) round(x, 0),  # integer labels
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 5,  # longer color bar
      title.position = "top"
    )
  ) +
  labs(
    title = "Gridded Cell Counts",
    x = "X Grid Index",
    y = "Y Grid Index",
    fill = "Count"
  ) +
  coord_fixed() +
  base_theme

# ---- Plot 2: Expression per Cell ----
p2 <- ggplot(df_expr, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(
    option = "magma",
    breaks = breaks,
    labels = function(x) round(x, 2),  # rounded expression
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 15,  # longer legend color bar
      title.position = "top"
    )
  ) +
  labs(
    title = paste0("Expression per Cell, ", colnames(center_slice_cell_location)[7]),
    x = "X Grid Index",
    y = "Y Grid Index",
    fill = "Avg Exprssion"
  ) +
  coord_fixed() +
  base_theme

# ---- Combine side-by-side ----
plot_grid(p1, p2, nrow=1, align='v')





# Plot: the correlation decay for each of the genes at three regions ------
# Figure 1: At DCIS region -----
cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_spatial.csv')
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000
grid_size<- 50
basis_choice <- 25
boundary_choice <- 1.5

center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]


gene_spatial_DCIS<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]

center_slice_cell_location$Gene<- gene_spatial_DCIS$CDH1



# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))

# the centroids of each grid: 
coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))




draws_df_HS_cdh1 <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis_CDH1.rds')
draws_df_HS_postn <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis_POSTN.rds')
draws_df_HS_trac <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis_TRAC.rds')




distMat <- fields::rdist(coords)
max_dist<- max(distMat)

# 1.2 Define Matérn covariance function for nu = 3/2
matern32 <- function(d, sigma, rho) {
  sigma^2*(1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}

summary_gene_cdh1<- as_draws_df(draws_df_HS_cdh1)
summary_gene_postn<- as_draws_df(draws_df_HS_postn)
summary_gene_trac<- as_draws_df(draws_df_HS_trac)





# corr decay plot -----
dist_choice<- seq(0, max_dist, by=0.1)

samples_rho1_cdh1<- summary_gene_cdh1$`ell[1]`
samples_rho2_cdh1<- summary_gene_cdh1$`ell[2]`
samples_a21_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*summary_gene_cdh1$`Astar[2,1]`
samples_a22_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*sqrt(1- summary_gene_cdh1$`Astar[2,1]`^2)

samples_rho1_postn<- summary_gene_postn$`ell[1]`
samples_rho2_postn<- summary_gene_postn$`ell[2]`
samples_a21_HS_postn<- summary_gene_postn$`sigma[2]`*summary_gene_postn$`Astar[2,1]`
samples_a22_HS_postn<- summary_gene_postn$`sigma[2]`*sqrt(1- summary_gene_postn$`Astar[2,1]`^2)

samples_rho1_trac<- summary_gene_trac$`ell[1]`
samples_rho2_trac<- summary_gene_trac$`ell[2]`
samples_a21_HS_trac<- summary_gene_trac$`sigma[2]`*summary_gene_trac$`Astar[2,1]`
samples_a22_HS_trac<- summary_gene_trac$`sigma[2]`*sqrt(1- summary_gene_trac$`Astar[2,1]`^2)



post_rho2_cdh1<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_cdh1))
post_rho2_postn<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_cdh1))
post_rho2_trac<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_cdh1))


post_mean_rho2_cdh1<- rep(NA, length(dist_choice))
post_upper_rho2_cdh1<- rep(NA, length(dist_choice))
post_lower_rho2_cdh1<- rep(NA, length(dist_choice))

post_mean_rho2_postn<- rep(NA, length(dist_choice))
post_upper_rho2_postn<- rep(NA, length(dist_choice))
post_lower_rho2_postn<- rep(NA, length(dist_choice))

post_mean_rho2_trac<- rep(NA, length(dist_choice))
post_upper_rho2_trac<- rep(NA, length(dist_choice))
post_lower_rho2_trac<- rep(NA, length(dist_choice))


for(k in 1:length(dist_choice)){
  # cdh1
  post_rho2_cdh1[k,]<- (samples_a21_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_cdh1)+
                     samples_a22_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_cdh1))/(samples_a21_HS_cdh1^2+samples_a22_HS_cdh1^2)
  # postn
  post_rho2_postn[k,]<- (samples_a21_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_postn)+
                          samples_a22_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_postn))/(samples_a21_HS_postn^2+samples_a22_HS_postn^2)
  # trac
  post_rho2_trac[k,]<- (samples_a21_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_trac)+
                          samples_a22_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_trac))/(samples_a21_HS_trac^2+samples_a22_HS_trac^2)
  
  }


for (k in 1:length(dist_choice)){
  # the cdh gene: 
  post_mean_rho2_cdh1[k]<-mean(post_rho2_cdh1[k,])
  post_upper_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.975)
  post_lower_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.025)
  # the postn gene: 
  post_mean_rho2_postn[k]<-mean(post_rho2_postn[k,])
  post_upper_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.975)
  post_lower_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.025)
  # the trac gene: 
  post_mean_rho2_trac[k]<-mean(post_rho2_trac[k,])
  post_upper_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.975)
  post_lower_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.025)
}



df_rho2<- cbind( dist_choice ,
                 post_mean_rho2_cdh1, 
                 post_upper_rho2_cdh1, 
                 post_lower_rho2_cdh1,
                 post_mean_rho2_postn, 
                 post_upper_rho2_postn, 
                 post_lower_rho2_postn,
                 post_mean_rho2_trac, 
                 post_upper_rho2_trac, 
                 post_lower_rho2_trac)
df_rho2<- as.data.frame(df_rho2)

colnames(df_rho2)<- c('dist_choice', 
                      'post_mean_rho2_cdh1', 
                      'post_upper_rho2_cdh1', 
                      'post_lower_rho2_cdh1',
                      'post_mean_rho2_postn', 
                      'post_upper_rho2_postn', 
                      'post_lower_rho2_postn',
                      'post_mean_rho2_trac', 
                      'post_upper_rho2_trac', 
                      'post_lower_rho2_trac')
est_range_2_cdh1<-dist_choice[which(abs(df_rho2[,2]-0.05)==min(abs(df_rho2[,2]-0.05)))]
est_range_2_postn<-dist_choice[which(abs(df_rho2[,5]-0.05)==min(abs(df_rho2[,5]-0.05)))]
est_range_2_trac<-dist_choice[which(abs(df_rho2[,8]-0.05)==min(abs(df_rho2[,8]-0.05)))]


eta2_corr_dcis<-ggplot(data = df_rho2, aes(x = dist_choice)) +
  geom_ribbon(aes(ymin = post_lower_rho2_cdh1, ymax = post_upper_rho2_cdh1, fill = "CDH1"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_cdh1, color = "CDH1"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_postn, ymax = post_upper_rho2_postn, fill = "POSTN"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_postn, color = "POSTN"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_trac, ymax = post_upper_rho2_trac, fill = "TRAC"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_trac, color = "TRAC"), linewidth = 1.05) +
  
  xlim(0, 2.5) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', col = '#ED1B2F', linewidth = 1.5) +
  
  # Annotation text for ranges - centered at x=1.8, increased spacing
  annotate("text", x = 1.5, y = 0.2, 
           label = paste0("CDH1: Range = ", round(est_range_2_cdh1, 2)), 
           color = "#8B6914", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.3, 
           label = paste0("POSTN: Range = ", round(est_range_2_postn, 2)), 
           color = "#00B050", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.4, 
           label = paste0("TRAC: Range = ", round(est_range_2_trac, 2)), 
           color = "#8B3A6E", size = 10, hjust = 0.5, fontface = "bold") +
  
  scale_color_manual(
    name = "Gene",
    values = c("CDH1" = "#8B6914", "POSTN" = "#00B050", "TRAC" = "#8B3A6E")
  ) +
  scale_fill_manual(
    name = "Gene",
    values = c("CDH1" = "#C9A855", "POSTN" = "#66D991", "TRAC" = "#B87BA2")
  ) +
  
  labs(x = '', y = expression('Posterior of ' ~ r[2])) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 30),
    legend.position = c(0.8, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = NA),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 25)
  )



eta2_corr_dcis


# 2. For Invasive tumor region-----

cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_spatial.csv')
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000
grid_size<- 50
basis_choice <- 25
boundary_choice <- 1.5

center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_spatial_invasive<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]


# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))

# the centroids of each grid: 
coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))


draws_df_HS_cdh1 <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/invasive_CDH1.rds')
draws_df_HS_postn <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/invasive_POSTN.rds')
draws_df_HS_trac <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/invasive_TRAC.rds')




distMat <- fields::rdist(coords)
max_dist<- max(distMat)

# 1.2 Define Matérn covariance function for nu = 3/2
matern32 <- function(d, sigma, rho) {
  sigma^2*(1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}

summary_gene_cdh1<- as_draws_df(draws_df_HS_cdh1)
summary_gene_postn<- as_draws_df(draws_df_HS_postn)
summary_gene_trac<- as_draws_df(draws_df_HS_trac)





# corr decay plot -----
dist_choice<- seq(0, max_dist, by=0.1)

samples_rho1_cdh1<- summary_gene_cdh1$`ell[1]`
samples_rho2_cdh1<- summary_gene_cdh1$`ell[2]`
samples_a21_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*summary_gene_cdh1$`Astar[2,1]`
samples_a22_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*sqrt(1- summary_gene_cdh1$`Astar[2,1]`^2)

samples_rho1_postn<- summary_gene_postn$`ell[1]`
samples_rho2_postn<- summary_gene_postn$`ell[2]`
samples_a21_HS_postn<- summary_gene_postn$`sigma[2]`*summary_gene_postn$`Astar[2,1]`
samples_a22_HS_postn<- summary_gene_postn$`sigma[2]`*sqrt(1- summary_gene_postn$`Astar[2,1]`^2)

samples_rho1_trac<- summary_gene_trac$`ell[1]`
samples_rho2_trac<- summary_gene_trac$`ell[2]`
samples_a21_HS_trac<- summary_gene_trac$`sigma[2]`*summary_gene_trac$`Astar[2,1]`
samples_a22_HS_trac<- summary_gene_trac$`sigma[2]`*sqrt(1- summary_gene_trac$`Astar[2,1]`^2)



post_rho2_cdh1<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_postn ))
post_rho2_postn<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho2_postn))
post_rho2_trac<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho2_postn))


post_mean_rho2_cdh1<- rep(NA, length(dist_choice))
post_upper_rho2_cdh1<- rep(NA, length(dist_choice))
post_lower_rho2_cdh1<- rep(NA, length(dist_choice))

post_mean_rho2_postn<- rep(NA, length(dist_choice))
post_upper_rho2_postn<- rep(NA, length(dist_choice))
post_lower_rho2_postn<- rep(NA, length(dist_choice))

post_mean_rho2_trac<- rep(NA, length(dist_choice))
post_upper_rho2_trac<- rep(NA, length(dist_choice))
post_lower_rho2_trac<- rep(NA, length(dist_choice))


for(k in 1:length(dist_choice)){
  # cdh1
  post_rho2_cdh1[k,]<- (samples_a21_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_cdh1)+
                          samples_a22_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_cdh1))/(samples_a21_HS_cdh1^2+samples_a22_HS_cdh1^2)
  # postn
  post_rho2_postn[k,]<- (samples_a21_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_postn)+
                           samples_a22_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_postn))/(samples_a21_HS_postn^2+samples_a22_HS_postn^2)
  # trac
  post_rho2_trac[k,]<- (samples_a21_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_trac)+
                          samples_a22_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_trac))/(samples_a21_HS_trac^2+samples_a22_HS_trac^2)
  
}


for (k in 1:length(dist_choice)){
  # the cdh gene: 
  post_mean_rho2_cdh1[k]<-mean(post_rho2_cdh1[k,])
  post_upper_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.975)
  post_lower_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.025)
  # the postn gene: 
  post_mean_rho2_postn[k]<-mean(post_rho2_postn[k,])
  post_upper_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.975)
  post_lower_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.025)
  # the trac gene: 
  post_mean_rho2_trac[k]<-mean(post_rho2_trac[k,])
  post_upper_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.975)
  post_lower_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.025)
}



df_rho2<- cbind( dist_choice ,
                 post_mean_rho2_cdh1, 
                 post_upper_rho2_cdh1, 
                 post_lower_rho2_cdh1,
                 post_mean_rho2_postn, 
                 post_upper_rho2_postn, 
                 post_lower_rho2_postn,
                 post_mean_rho2_trac, 
                 post_upper_rho2_trac, 
                 post_lower_rho2_trac)
df_rho2<- as.data.frame(df_rho2)

colnames(df_rho2)<- c('dist_choice', 
                      'post_mean_rho2_cdh1', 
                      'post_upper_rho2_cdh1', 
                      'post_lower_rho2_cdh1',
                      'post_mean_rho2_postn', 
                      'post_upper_rho2_postn', 
                      'post_lower_rho2_postn',
                      'post_mean_rho2_trac', 
                      'post_upper_rho2_trac', 
                      'post_lower_rho2_trac')
est_range_2_cdh1<-dist_choice[which(abs(df_rho2[,2]-0.05)==min(abs(df_rho2[,2]-0.05)))]
est_range_2_postn<-dist_choice[which(abs(df_rho2[,5]-0.05)==min(abs(df_rho2[,5]-0.05)))]
est_range_2_trac<-dist_choice[which(abs(df_rho2[,8]-0.05)==min(abs(df_rho2[,8]-0.05)))]


eta2_corr_invasive<-ggplot(data = df_rho2, aes(x = dist_choice)) +
  geom_ribbon(aes(ymin = post_lower_rho2_cdh1, ymax = post_upper_rho2_cdh1, fill = "CDH1"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_cdh1, color = "CDH1"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_postn, ymax = post_upper_rho2_postn, fill = "POSTN"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_postn, color = "POSTN"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_trac, ymax = post_upper_rho2_trac, fill = "TRAC"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_trac, color = "TRAC"), linewidth = 1.05) +
  
  xlim(0, 2.5) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', col = '#ED1B2F', linewidth = 1.5) +
  
  # Annotation text for ranges - centered at x=1.8, increased spacing
  annotate("text", x = 1.5, y = 0.2, 
           label = paste0("CDH1: Range = ", round(est_range_2_cdh1, 2)), 
           color = "#8B6914", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.3, 
           label = paste0("POSTN: Range = ", round(est_range_2_postn, 2)), 
           color = "#00B050", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.4, 
           label = paste0("TRAC: Range = ", round(est_range_2_trac, 2)), 
           color = "#8B3A6E", size = 10, hjust = 0.5, fontface = "bold") +
  
  scale_color_manual(
    name = "Gene",
    values = c("CDH1" = "#8B6914", "POSTN" = "#00B050", "TRAC" = "#8B3A6E")
  ) +
  scale_fill_manual(
    name = "Gene",
    values = c("CDH1" = "#C9A855", "POSTN" = "#66D991", "TRAC" = "#B87BA2")
  ) +
  
  labs(x = 'Distance (Millimeter)', y = expression('Posterior of ' ~ r[2])) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 30),
    legend.position = 'none',
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 25)
  )



eta2_corr_invasive





# 3. For Immune region-----

cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/cell_spatial.csv')
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000
grid_size<- 50
basis_choice <- 25
boundary_choice <- 1.5
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),]
gene_spatial_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]


# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))

# the centroids of each grid: 
coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))


draws_df_HS_cdh1 <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/immune_CDH1.rds')
draws_df_HS_postn <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/immune_POSTN.rds')
draws_df_HS_trac <- readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/immune_TRAC.rds')




distMat <- fields::rdist(coords)
max_dist<- max(distMat)

# 1.2 Define Matérn covariance function for nu = 3/2
matern32 <- function(d, sigma, rho) {
  sigma^2*(1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}

summary_gene_cdh1<- as_draws_df(draws_df_HS_cdh1)
summary_gene_postn<- as_draws_df(draws_df_HS_postn)
summary_gene_trac<- as_draws_df(draws_df_HS_trac)





# corr decay plot -----
dist_choice<- seq(0, max_dist, by=0.1)

samples_rho1_cdh1<- summary_gene_cdh1$`ell[1]`
samples_rho2_cdh1<- summary_gene_cdh1$`ell[2]`
samples_a21_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*summary_gene_cdh1$`Astar[2,1]`
samples_a22_HS_cdh1<- summary_gene_cdh1$`sigma[2]`*sqrt(1- summary_gene_cdh1$`Astar[2,1]`^2)

samples_rho1_postn<- summary_gene_postn$`ell[1]`
samples_rho2_postn<- summary_gene_postn$`ell[2]`
samples_a21_HS_postn<- summary_gene_postn$`sigma[2]`*summary_gene_postn$`Astar[2,1]`
samples_a22_HS_postn<- summary_gene_postn$`sigma[2]`*sqrt(1- summary_gene_postn$`Astar[2,1]`^2)

samples_rho1_trac<- summary_gene_trac$`ell[1]`
samples_rho2_trac<- summary_gene_trac$`ell[2]`
samples_a21_HS_trac<- summary_gene_trac$`sigma[2]`*summary_gene_trac$`Astar[2,1]`
samples_a22_HS_trac<- summary_gene_trac$`sigma[2]`*sqrt(1- summary_gene_trac$`Astar[2,1]`^2)



post_rho2_cdh1<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_postn))
post_rho2_postn<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_postn))
post_rho2_trac<- matrix(NA, nrow=length(dist_choice), ncol=length(samples_rho1_postn))


post_mean_rho2_cdh1<- rep(NA, length(dist_choice))
post_upper_rho2_cdh1<- rep(NA, length(dist_choice))
post_lower_rho2_cdh1<- rep(NA, length(dist_choice))

post_mean_rho2_postn<- rep(NA, length(dist_choice))
post_upper_rho2_postn<- rep(NA, length(dist_choice))
post_lower_rho2_postn<- rep(NA, length(dist_choice))

post_mean_rho2_trac<- rep(NA, length(dist_choice))
post_upper_rho2_trac<- rep(NA, length(dist_choice))
post_lower_rho2_trac<- rep(NA, length(dist_choice))


for(k in 1:length(dist_choice)){
  # cdh1
  post_rho2_cdh1[k,]<- (samples_a21_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_cdh1)+
                          samples_a22_HS_cdh1^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_cdh1))/(samples_a21_HS_cdh1^2+samples_a22_HS_cdh1^2)
  # postn
  post_rho2_postn[k,]<- (samples_a21_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_postn)+
                           samples_a22_HS_postn^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_postn))/(samples_a21_HS_postn^2+samples_a22_HS_postn^2)
  # trac
  post_rho2_trac[k,]<- (samples_a21_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho1_trac)+
                          samples_a22_HS_trac^2*matern32(d=dist_choice[k], sigma=1, rho= samples_rho2_trac))/(samples_a21_HS_trac^2+samples_a22_HS_trac^2)
  
}


for (k in 1:length(dist_choice)){
  # the cdh gene: 
  post_mean_rho2_cdh1[k]<-mean(post_rho2_cdh1[k,])
  post_upper_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.975)
  post_lower_rho2_cdh1[k]<-quantile(post_rho2_cdh1[k,], 0.025)
  # the postn gene: 
  post_mean_rho2_postn[k]<-mean(post_rho2_postn[k,])
  post_upper_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.975)
  post_lower_rho2_postn[k]<-quantile(post_rho2_postn[k,], 0.025)
  # the trac gene: 
  post_mean_rho2_trac[k]<-mean(post_rho2_trac[k,])
  post_upper_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.975)
  post_lower_rho2_trac[k]<-quantile(post_rho2_trac[k,], 0.025)
}



df_rho2<- cbind( dist_choice ,
                 post_mean_rho2_cdh1, 
                 post_upper_rho2_cdh1, 
                 post_lower_rho2_cdh1,
                 post_mean_rho2_postn, 
                 post_upper_rho2_postn, 
                 post_lower_rho2_postn,
                 post_mean_rho2_trac, 
                 post_upper_rho2_trac, 
                 post_lower_rho2_trac)
df_rho2<- as.data.frame(df_rho2)

colnames(df_rho2)<- c('dist_choice', 
                      'post_mean_rho2_cdh1', 
                      'post_upper_rho2_cdh1', 
                      'post_lower_rho2_cdh1',
                      'post_mean_rho2_postn', 
                      'post_upper_rho2_postn', 
                      'post_lower_rho2_postn',
                      'post_mean_rho2_trac', 
                      'post_upper_rho2_trac', 
                      'post_lower_rho2_trac')
est_range_2_cdh1<-dist_choice[which(abs(df_rho2[,2]-0.05)==min(abs(df_rho2[,2]-0.05)))]
est_range_2_postn<-dist_choice[which(abs(df_rho2[,5]-0.05)==min(abs(df_rho2[,5]-0.05)))]
est_range_2_trac<-dist_choice[which(abs(df_rho2[,8]-0.05)==min(abs(df_rho2[,8]-0.05)))]


eta2_corr_immune<-ggplot(data = df_rho2, aes(x = dist_choice)) +
  geom_ribbon(aes(ymin = post_lower_rho2_cdh1, ymax = post_upper_rho2_cdh1, fill = "CDH1"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_cdh1, color = "CDH1"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_postn, ymax = post_upper_rho2_postn, fill = "POSTN"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_postn, color = "POSTN"), linewidth = 1.05) +
  
  geom_ribbon(aes(ymin = post_lower_rho2_trac, ymax = post_upper_rho2_trac, fill = "TRAC"), alpha = 0.3) +
  geom_line(aes(y = post_mean_rho2_trac, color = "TRAC"), linewidth = 1.05) +
  
  xlim(0, 2.5) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', col = '#ED1B2F', linewidth = 1.5) +
  
  # Annotation text for ranges - centered at x=1.8, increased spacing
  annotate("text", x = 1.5, y = 0.2, 
           label = paste0("CDH1: Range = ", round(est_range_2_cdh1, 2)), 
           color = "#8B6914", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.3, 
           label = paste0("POSTN: Range = ", round(est_range_2_postn, 2)), 
           color = "#00B050", size = 10, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.4, 
           label = paste0("TRAC: Range = ", round(est_range_2_trac, 2)), 
           color = "#8B3A6E", size = 10, hjust = 0.5, fontface = "bold") +
  
  scale_color_manual(
    name = "Gene",
    values = c("CDH1" = "#8B6914", "POSTN" = "#00B050", "TRAC" = "#8B3A6E")
  ) +
  scale_fill_manual(
    name = "Gene",
    values = c("CDH1" = "#C9A855", "POSTN" = "#66D991", "TRAC" = "#B87BA2")
  ) +
  
  #labs(x = 'Distance (Millimeter)', y = expression('Posterior of ' ~ r[2])) +
  labs(x = '', y = expression('Posterior of ' ~ r[2])) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = 'none',
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 25)
  )



eta2_corr_immune

plot_grid(eta2_corr_dcis, 
          eta2_corr_invasive, 
          eta2_corr_immune, ncol=1, align='v')




# Sankey plot for large range analysis: ------
library(ggplot2)
library(ggalluvial)
# load the data 
dcis_neg<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/dcis-neg.csv')[,2]
dcis_pos<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/dcis-pos.csv')[,2]

immune_neg<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/immune-neg.csv')[,2]
immune_pos<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/immune-pos.csv')[,2]

tumor_neg<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/tumor-neg.csv')[,2]
tumor_pos<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/tumor-pos.csv')[,2]

groups <- list(
  DCIS_pos   = dcis_pos,
  DCIS_neg   = dcis_neg,
  Tumor_pos  = tumor_pos,
  Tumor_neg  = tumor_neg
)

immune_groups <- list(
  Immune_pos = immune_pos,
  Immune_neg = immune_neg
)

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggalluvial)
library(stringr)
library(cowplot)

# ==================== PLOT 1: DCIS/Tumor to Immune ====================
flows <- map_df(names(groups), function(g_left) {
  map_df(names(immune_groups), function(g_right) {
    genes_left  <- groups[[g_left]]
    genes_right <- immune_groups[[g_right]]
    
    data.frame(
      From = g_left,
      To   = g_right,
      Flow = length(intersect(genes_left, genes_right))
    )
  })
})

# Recolor the path and set ordering
flows <- flows %>%
  mutate(
    sign_pair = case_when(
      grepl("_pos$", From) & grepl("_neg$", To) ~ "mixed",
      grepl("_neg$", From) & grepl("_pos$", To) ~ "mixed",
      grepl("_pos$", From) & grepl("_pos$", To) ~ "same",
      grepl("_neg$", From) & grepl("_neg$", To) ~ "same",
      TRUE ~ NA_character_
    )
  ) %>%
  # Manually set the order: mixed, same, mixed, same for each group
  mutate(
    flow_sequence = case_when(
      From == "DCIS_pos" & To == "Immune_neg" ~ 1,   # mixed (red)
      From == "DCIS_pos" & To == "Immune_pos" ~ 2,   # same (green)
      From == "DCIS_neg" & To == "Immune_pos" ~ 3,   # mixed (red)
      From == "DCIS_neg" & To == "Immune_neg" ~ 4,   # same (green)
      From == "Tumor_pos" & To == "Immune_neg" ~ 5,  # mixed (red)
      From == "Tumor_pos" & To == "Immune_pos" ~ 6,  # same (green)
      From == "Tumor_neg" & To == "Immune_pos" ~ 7,  # mixed (red)
      From == "Tumor_neg" & To == "Immune_neg" ~ 8,  # same (green)
      TRUE ~ NA_real_
    ),
    From = factor(From, levels = c("DCIS_pos", "DCIS_neg", "Tumor_pos", "Tumor_neg")),
    To = factor(To, levels = c("Immune_pos", "Immune_neg"))
  ) %>%
  arrange(flow_sequence)

flows$Flow <- as.numeric(flows$Flow)

pf1 <- ggplot(flows,
              aes(axis1 = From, axis2 = To, y = Flow)) +
  
  # Alluvial flows with ordering
  geom_alluvium(
    aes(fill = sign_pair, order = flow_sequence),
    width = 1/12,
    alpha = 0.85
  ) +
  
  # Grey strata WITHOUT text
  geom_stratum(
    width = 0.25,
    fill = "grey95",
    color = "grey40"
  ) +
  
  # Flow numbers
  geom_text(
    stat = "alluvium",
    aes(label = Flow, order = flow_sequence),
    size = 5,
    color = "black"
  ) +
  
  # Top label: DCIS / Tumor
  annotate(
    "text",
    x = 1,
    y = Inf,
    label = "DCIS / Tumor",
    vjust = 1.5,
    size = 4.5,
    fontface = "bold"
  ) +
  
  # Top label: Immune
  annotate(
    "text",
    x = 2,
    y = Inf,
    label = "Immune",
    vjust = 1.5,
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "mixed" = "#FFB3A5",
      "same"  = "#7ABD91"
    ),
    name = "Sign consistency",
    labels = c(
      "mixed" = "Discordant",
      "same"  = "Concordant"
    )
  ) +
  
  scale_x_discrete(
    limits = c("DCIS / Tumor", "Immune"),
    expand = c(.25, .25)
  ) +
  
  theme_void(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(t = 20, r = 2, b = 10, l = 15),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  
  labs(
    y = "Number of shared genes"
  )


# ==================== PLOT 2: DCIS to Tumor ====================
dcis_groups <- list(
  DCIS_pos = dcis_pos,
  DCIS_neg = dcis_neg
)

tumor_groups <- list(
  Tumor_pos = tumor_pos,
  Tumor_neg = tumor_neg
)

# Create flows between DCIS and Tumor
flows2 <- map_df(names(dcis_groups), function(g_left) {
  map_df(names(tumor_groups), function(g_right) {
    genes_left  <- dcis_groups[[g_left]]
    genes_right <- tumor_groups[[g_right]]
    
    data.frame(
      From = g_left,
      To   = g_right,
      Flow = length(intersect(genes_left, genes_right))
    )
  })
})

# Recolor the path and set ordering
flows2 <- flows2 %>%
  mutate(
    sign_pair = case_when(
      grepl("_pos$", From) & grepl("_neg$", To) ~ "mixed",
      grepl("_neg$", From) & grepl("_pos$", To) ~ "mixed",
      grepl("_pos$", From) & grepl("_pos$", To) ~ "same",
      grepl("_neg$", From) & grepl("_neg$", To) ~ "same",
      TRUE ~ NA_character_
    )
  ) %>%
  # Manually set the desired order: mixed, same, mixed, same
  mutate(
    flow_sequence = case_when(
      From == "DCIS_pos" & To == "Tumor_neg" ~ 1,  # mixed (red)
      From == "DCIS_pos" & To == "Tumor_pos" ~ 2,  # same (green)
      From == "DCIS_neg" & To == "Tumor_pos" ~ 3,  # mixed (red)
      From == "DCIS_neg" & To == "Tumor_neg" ~ 4,  # same (green)
      TRUE ~ NA_real_
    ),
    From = factor(From, levels = c("DCIS_pos", "DCIS_neg")),
    To = factor(To, levels = c("Tumor_pos", "Tumor_neg"))
  ) %>%
  arrange(flow_sequence)

flows2$Flow <- as.numeric(flows2$Flow)

pf2 <- ggplot(flows2,
              aes(axis1 = From, axis2 = To, y = Flow)) +
  
  # Alluvial flows with ordering
  geom_alluvium(
    aes(fill = sign_pair, order = flow_sequence),
    width = 1/12,
    alpha = 0.85
  ) +
  
  # Grey strata WITHOUT text
  geom_stratum(
    width = 0.25,
    fill = "grey95",
    color = "grey40"
  ) +
  
  # Flow numbers
  geom_text(
    stat = "alluvium",
    aes(label = Flow, order = flow_sequence),
    size = 4.5,
    color = "black"
  ) +
  
  # Top label: DCIS
  annotate(
    "text",
    x = 1,
    y = Inf,
    label = "DCIS",
    vjust = 1.5,
    size = 4.5,
    fontface = "bold"
  ) +
  
  # Top label: Tumor
  annotate(
    "text",
    x = 2,
    y = Inf,
    label = "Tumor",
    vjust = 1.5,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "mixed" = "#FFB3A5",
      "same"  = "#7ABD91"
    ),
    name = "Sign consistency",
    labels = c(
      "mixed" = "Discordant",
      "same"  = "Concordant"
    )
  ) +
  
  scale_x_discrete(
    limits = c("DCIS", "Tumor"),
    expand = c(.25, .25)
  ) +
  
  theme_void(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(t = 20, r = 15, b = 10, l = 2),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  
  labs(
    y = "Number of shared genes"
  )


# ==================== COMBINE PLOTS ====================
# Create plots WITHOUT legends and with reduced margins
pf1_no_legend <- pf1 + 
  theme(
    legend.position = "none",
    plot.margin = margin(t = 20, r = 5, b = 10, l = 15)
  )

pf2_no_legend <- pf2 + 
  theme(
    legend.position = "none",
    plot.margin = margin(t = 20, r = 15, b = 10, l = 5)
  )

# Extract the legend from one of the plots
legend <- get_legend(
  pf1 + 
    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

# Combine: legend on top, plots below with reduced spacing
plot_grid(
  legend,
  plot_grid(pf1_no_legend, pf2_no_legend, ncol = 2, rel_widths = c(1, 1)),
  ncol = 1,
  rel_heights = c(0.05, 1)
)
