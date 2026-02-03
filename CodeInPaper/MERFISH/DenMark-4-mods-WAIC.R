####################################################
# Test for WAIC in the 4 mod: if they work 
####################################################

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

load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-location.Rdata')
load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-gene-expression.Rdata')


center_slice_cell_location<-  center_slice_cell_location/1000
GeneIDls<- c(387, 403,  73, 435, 421 )

# model 2:----
GeneID<- GeneIDls[5]
# the gene name
rownames(center_slice_gene_location)[GeneID]


#  scale the center slice 
df<- cbind( center_slice_cell_location ,center_slice_gene_location[GeneID,])
colnames(df)<- c('x','y','Gene')


set.seed(1)
options(cmdstanr_max_rows = 20)

## 1.1 define grid and distance matrix----
# Define grid size and coordinates 

grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


length(loc$x)==length(eg_gene_express) # check the length are the same 


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



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


Y<- loc_100_v
M<- marks_100_v


# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))


coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))



distMat <- fields::rdist(coords)
mrange<-max(distMat)/(2*2.75)

# the area 
grid_area<- (2*L1*2*L2)/(grid_size^2)
log_grid_area <- log(grid_area) 

################################################
# HSGP approximation to general model -----
################################################

xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# 22*22
m1 <- 25; m2 <- 25; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)



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
  sigma=c(1.3, 1.1),
  beta0=0, 
  beta1=-1, 
  ell=c(0.5, 0.5),
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


# Extract draws for log_lik
log_lik_array <- cmdstan_fit_HS$draws("log_lik") # 3D array: iterations × chains × n

# Convert to matrix for loo
log_lik_mat <- posterior::as_draws_matrix(log_lik_array) # iterations × n


# Compute WAIC
waic_result <- loo::waic(log_lik_mat)


save_df<- as.data.frame(cbind(fit_summary_HS,waic_result$waic ))


# write the csv files 
write.csv(save_df, 
          paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/',
                 rownames(center_slice_gene_location)[GeneID], '_M2_final.csv'))

# also save the samples: 

# save the Aqp4 as the file: 
saveRDS(cmdstan_fit_HS, 
        file = paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/cmdstan_fit_M2_',
                      rownames(center_slice_gene_location)[GeneID], '.rds'))
          
 
          
          
          
# load the saved data: 
#cmdstan_fit_HS <- readRDS("C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/cmdstan_fit_M2_Aqp4.rds")





# Model 1: -----
GeneID<- GeneIDls[1]
# the gene name
rownames(center_slice_gene_location)[GeneID]


#  scale the center slice 
df<- cbind( center_slice_cell_location ,center_slice_gene_location[GeneID,])
colnames(df)<- c('x','y','Gene')


set.seed(1)
options(cmdstanr_max_rows = 20)

## 1.1 define grid and distance matrix----
# Define grid size and coordinates 

grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


length(loc$x)==length(eg_gene_express) # check the length are the same 


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



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


Y<- loc_100_v
M<- marks_100_v


# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))


coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))



distMat <- fields::rdist(coords)

# the area 
grid_area<- (2*L1*2*L2)/(grid_size^2)
log_grid_area <- log(grid_area) 

# adding the grid_area term,
input <- list(n = N,
              y1 = Y, # outcome1: points 
              y2 = M, # outcome2: marks 
              log_grid_area=log_grid_area)

# define initial values 
init_fun <- function()list(
  beta0=0, 
  beta1=-1, 
  beta_y=1)


str(input)

stan_file <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/M1-fixed-effect-mod.stan'
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


fit_summary_HS <- cmdstan_fit_HS$summary(variables = c("beta0", "beta1",  'beta_y'), 
                                         c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS


# Extract draws for log_lik
log_lik_array <- cmdstan_fit_HS$draws("log_lik") # 3D array: iterations × chains × n

# Convert to matrix for loo
log_lik_mat <- posterior::as_draws_matrix(log_lik_array) # iterations × n


# Compute WAIC
waic_result <- loo::waic(log_lik_mat)


save_df<- as.data.frame(cbind(fit_summary_HS,waic_result$waic ))


# write the csv files 
write.csv(save_df, 
          paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/',
                 rownames(center_slice_gene_location)[GeneID], '_M1.csv'))






# Model 1-revised (BenchMark)-----------

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

load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-location.Rdata')
load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-gene-expression.Rdata')


center_slice_cell_location<-  center_slice_cell_location/1000
GeneIDls<- c(387, 403,  73, 435, 421 )

# changing from this 

GeneID<- GeneIDls[5]
# the gene name
rownames(center_slice_gene_location)[GeneID]


#  scale the center slice 
df<- cbind( center_slice_cell_location ,center_slice_gene_location[GeneID,])
colnames(df)<- c('x','y','Gene')


set.seed(1)
options(cmdstanr_max_rows = 20)

## 1.1 define grid and distance matrix----
# Define grid size and coordinates 

grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


length(loc$x)==length(eg_gene_express) # check the length are the same 


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



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


Y<- loc_100_v
M<- marks_100_v

# normalization of the locations:
L1<- max(abs(center_slice_cell_location$x -(min(center_slice_cell_location$x)+max(center_slice_cell_location$x))/2 ))
L2<- max(abs(center_slice_cell_location$y -(min(center_slice_cell_location$y)+max(center_slice_cell_location$y))/2 ))


coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size), 
                                       y = seq(-L2, L2, length.out = grid_size))))



distMat <- fields::rdist(coords)
mrange<-max(distMat)/(2*2.75)

# the area 
grid_area<- (2*L1*2*L2)/(grid_size^2)
log_grid_area <- log(grid_area) 


################################################
# HSGP approximation to general model -----
################################################

xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# 22*22
m1 <- 25; m2 <- 25; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)



# adding the grid_area term,
input <- list(n = N,
              y1 = Y, # outcome1: points 
              y2 = M, # outcome2: marks 
              log_grid_area=log_grid_area,
              d=2,
              mstar = mstar, 
              coords = coords,
              mrange=mrange,
              indices = S,
              L=L,
              is_centerted_PHI = 0)




# define initial values 
init_fun <- function()list(
  a11=1,
  a22=0.8,
  sigma=c(1.3, 1.1),
  beta0=0, 
  beta1=-1, 
  ell=c(2, 2),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)



str(input)

stan_file <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/M1-RE-basic-mod.stan'
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



fit_summary_HS <- cmdstan_fit_HS$summary(variables = c("beta0", "beta1",  "ell[1]", 'ell[2]', 
                                                       'a11','a22'), 
                                         c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS


# Extract draws for log_lik
log_lik_array <- cmdstan_fit_HS$draws("log_lik") # 3D array: iterations × chains × n

# Convert to matrix for loo
log_lik_mat <- posterior::as_draws_matrix(log_lik_array) # iterations × n


# Compute WAIC
waic_result <- loo::waic(log_lik_mat)
print(waic_result)



save_df<- as.data.frame(cbind(fit_summary_HS,waic_result$waic ))


# write the csv files 
write.csv(save_df, 
          paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/',
                 rownames(center_slice_gene_location)[GeneID], '_M1-rev.csv'))













