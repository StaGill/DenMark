#######################################################################
# Compare the point-level and grid-level data 
# The general model is run on the M2, which is 
# log(lambda1(g)) = beta0 + a11*w1(g)
# log(lambda2(g)) = beta1 + log(Y) + a21*w1(g) + a22*w2(g)
#######################################################################

# 0. input the parameter 
args=commandArgs(trailingOnly = TRUE)

# the packages
library('rstan')
library('MASS')
library(fields)
library(tidybayes)
library(coda)
library(cmdstanr)
library(RandomFields)
library(spatstat)
library(spatstat.random)
library(bayesplot)
library(ggplot2)
library(cowplot)


# 
if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}


# the iter parameters
sim_num <- as.numeric(args[1])
iter<- sim_num
print(iter)

set.seed(iter)



#set.seed(1)
options(cmdstanr_max_rows = 20)

# build a function to discretize the MPP with grids:
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


# the functions 
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


# From the true parameter to simulate the MPP -----

#================ parameters (gene 18) ========
# increase the windows 
# the beta_y is fixed at 1 only, and check the fitting criteria 

beta_pts<- -1 # 
beta_marks <- -1  
lscale <-  c(1.5, 1) # the range parameter 
sigma<- c(2, 1)
# the correlation matrix 
R <- cbind(c(1,0.8),
           c(0.8,1))
# elements of A: 
a11 <- sigma[1]
a21 <- sigma[2]*R[2,1]
a22 <- sigma[2]*sqrt(1-(R[2,1])^2) 
# window: 
L1 <- 4.44
L2 <- 3.56
#============================================


eigen(R)$values


# Simulate the MPP ----

# define the window in the real-dataset
win_length_x<- 2*L1 # 2*half-range-x 
win_length_y<- 2*L2 # 2*half-range-y

win<- owin(c(-L1,L1), 
           c(-L2,L2))

# A. simulate points ----
large_res<- 100 # can change the pixel resolution to generate a reasonable number of cells
spatstat.options(npixel = large_res)
x_coords <- seq(-L1, L1, length.out = large_res)
y_coords <- seq(-L2, L2, length.out = large_res)
grid <- expand.grid(x = x_coords, y = y_coords)

#spatstat.options(npixel = large_res)
#exp(beta_pts) * diff(range(win$x)) * diff(range(win$y)) # expected points 

# set seeds 
params <- list(var = a11^2, scale = lscale, nu = 3/2)
# use rLGCP to simulate a LGCP 
sim.pts<-rLGCP(model='matern', 
               mu=beta_pts,
               param= params,
               win=win,
               saveLambda = TRUE )


coords_pp <- unmark(sim.pts) 

# extract the simulated lambda 
sim.lambda1 <- attr(sim.pts,  "Lambda")
sim.log.lambda1<- log(sim.lambda1)
gp_at_points<- sim.log.lambda1[coords_pp]

simulate_pp<- cbind(sim.pts$x,  
                    sim.pts$y)

sim.pts$n # number of simulated points, here is 2,817 points


#plot(sim.pts$x, sim.pts$y)

# check the simulate gaussian process fro the points: w1(s)
pts.df<- cbind(sim.pts$x, sim.pts$y,  (gp_at_points-beta_pts)/a11 )
colnames(pts.df)<- c('x','y','omega1')


# save the log(lambda1)
log_lambda1 <- gp_at_points

## simulate the marks:
RFoptions(spConform=F, seed=iter+1)

matern.model<- RMmatern(nu=1.5, 
                        var=1, 
                        scale=lscale[2])

sim.omega2<- RFsimulate(matern.model, 
                        x=simulate_pp[,1],
                        y=simulate_pp[,2])


# simulate the intensity of marks: 
log_lambda2 <- beta_marks + a21*((gp_at_points-beta_pts)/a11) +a22*sim.omega2

# the exponential of the latent field 
sim.lambda2<- exp(log_lambda2)



# III. Then simulate marks: gene expression---
sim.marks.v<- rpois(length(sim.lambda2), 
                    sim.lambda2)

# attach marks to points 
coords_pp_marks <- cbind(pts.df[,1], pts.df[,2], sim.marks.v)
colnames(coords_pp_marks)<- c('x','y','mark')


# To sum up: 
# the generated mpp dataframe: coords_pp_marks
coords_pp_marks <- as.data.frame(coords_pp_marks)
colnames(coords_pp_marks)<- c('x', 'y', 'mark')





# B/ Fit with different grid sizes -----


# Build grids for all 4 resolutions 
# (50, 40, 30, 20) 
# 1. analyze with 50*50 grids ----
grid_size_large <- 50 

grid_summary_large <- grid_pts_centroid(grid_size=grid_size_large, 
                                        dataset=coords_pp_marks)

grid_data_large <- buildgridpp(grid_size = grid_size_large, 
                               dataset = coords_pp_marks)

#plot(grid_summary_large[[1]][,1], grid_summary_large[[1]][,2], cex= grid_data_large[[1]][,1])

# 2. analyze with 40*40 grids ----
grid_size_mid <- 40 

grid_summary_mid <- grid_pts_centroid(grid_size=grid_size_mid, 
                                      dataset=coords_pp_marks)

grid_data_mid <- buildgridpp(grid_size = grid_size_mid, 
                             dataset = coords_pp_marks)

# 3. analyze with 30*30 grids ----
grid_size_mid2 <- 30 

grid_summary_mid2 <- grid_pts_centroid(grid_size=grid_size_mid2, 
                                       dataset=coords_pp_marks)

grid_data_mid2 <- buildgridpp(grid_size = grid_size_mid2, 
                              dataset = coords_pp_marks)

# 4. analyze with 20*20 grids ----
grid_size_low <- 20

grid_summary_low <- grid_pts_centroid(grid_size=grid_size_low, 
                                      dataset=coords_pp_marks)

grid_data_low <- buildgridpp(grid_size = grid_size_low, 
                             dataset = coords_pp_marks)


# fit with different grid resolutions ----

# To sum up, we will have: 
# The original res data, high and low resolution data : 
# 50*50, 40*40, 30*30, 20*20, as vectors 
Y_high <-  grid_data_large[[1]][,1]
Y_mid <- grid_data_mid[[1]][,1]
Y_mid2<- grid_data_mid2[[1]][,1]
Y_low<- grid_data_low[[1]][,1]

M_high <-  grid_data_large[[1]][,2]
M_mid <- grid_data_mid[[1]][,2]
M_mid2<- grid_data_mid2[[1]][,2]
M_low<- grid_data_low[[1]][,2]


# the number of grids 
n<-  grid_summary_large[[3]]
n1<- grid_summary_mid[[3]]
n2<- grid_summary_mid2[[3]]
n3<- grid_summary_low[[3]]


# and the grid areas for the three conditions 
grid_area_high <- grid_summary_large[[2]]
grid_area_mid  <- grid_summary_mid[[2]]
grid_area_mid2 <- grid_summary_mid2[[2]]
grid_area_low  <- grid_summary_low[[2]]

# then fit with different resolution cases ----
#stan_file_path<- '/lustre09/project/6003552/mingchi/git/M2-latent-field-mod-S2.stan'
stan_file_path<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation2/general-model/M2-latent-field-mod-S2-v2.stan'


# The following analysis is without adjusting the optimal choice of basis functions 

# A. Analyze with HSGP, highest resolution, -----

# change the window length if possible 
xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# adjust the basis function to find proper number/ 22
#m1 <- initial_m[1]; m2 <- initial_m[2]; mstar <- m1*m2
m1<- 20; m2<- 20; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
#c<- rep(initial_c, 2)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)

mlscale<- max(dist(grid_summary_large[[1]]))/(2*2.75)

# adding the grid_area term,
input <- list(n = n,
              y1 = Y_high, # outcome1: points
              y2 = M_high,
              d = 2, # dim=2
              mstar = mstar, 
              coords =  grid_summary_large[[1]],
              indices = S,
              #mrange=mlscale,
              L=L,
              log_grid_area=log(grid_area_high),
              is_centerted_PHI = 0)

# define initial values 
init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 1), 2, 2),
  beta0=2, 
  beta1=1, 
  sigma=c(1, 1),
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)

#str(input)
stan_file <- stan_file_path
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_high <- mod$sample(data = input, 
                               chains = 4,
                               parallel_chains = 4,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               adapt_delta = 0.95,
                               init=init_fun, 
                               max_treedepth = 10,
                               step_size = 1)
elapsed_time <- cmdstan_fit_high$time()
elapsed_time
elapsed_time$total/60


fit_summary_HS_high <- cmdstan_fit_high$summary(variables = c("beta0", 'beta1', "ell[1]",'ell[2]', 'Astar[2,1]'), 
                                                c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_high




# B. Analyze with HSGP, middle resolution -----

# change the window length if possible 
xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# adjust the basis function to find proper number/ 22
#m1 <- initial_m[1]; m2 <- initial_m[2]; mstar <- m1*m2
m1<- 20; m2<- 20; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
#c<- rep(initial_c, 2)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)

mlscale<- max(dist(grid_summary_mid[[1]]))/(2*2.75)

# adding the grid_area term,
input <- list(n = n1,
              y1 = Y_mid, # outcome1: points
              y2 = M_mid,
              d = 2, # dim=2
              mstar = mstar, 
              coords =  grid_summary_mid[[1]],
              indices = S,
              L=L,
              log_grid_area=log(grid_area_mid),
              is_centerted_PHI = 0)

# define initial values 
init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 1), 2, 2),
  beta0=2, 
  beta1=1, 
  sigma=c(1, 1),
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)

#str(input)
stan_file <- stan_file_path
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_mid <- mod$sample(data = input, 
                              chains = 4,
                              parallel_chains = 4,
                              iter_warmup = 1000,
                              iter_sampling = 1000,
                              adapt_delta = 0.95,
                              init=init_fun, 
                              max_treedepth = 10,
                              step_size = 1)
elapsed_time <- cmdstan_fit_mid$time()
elapsed_time
elapsed_time$total/60


fit_summary_HS_mid <- cmdstan_fit_mid$summary(variables = c("beta0", "beta1",  "ell[1]", 'ell[2]', 
                                                            'Astar[2,1]'), 
                                              c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_mid





# C. Analyze with HSGP, second middle resolution -----

# change the window length if possible 
xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# adjust the basis function to find proper number/ 22
#m1 <- initial_m[1]; m2 <- initial_m[2]; mstar <- m1*m2
m1<- 20; m2<- 20; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
#c<- rep(initial_c, 2)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)

mlscale<- max(dist(grid_summary_mid2[[1]]))/(2*2.75)

# adding the grid_area term,
input <- list(n = n2,
              y1 = Y_mid2, # outcome1: points
              y2 = M_mid2,
              d = 2, # dim=2
              mstar = mstar, 
              coords =  grid_summary_mid2[[1]],
              indices = S,
              L=L,
              log_grid_area=log(grid_area_mid2),
              is_centerted_PHI = 0)

# define initial values 
init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 1), 2, 2),
  beta0=2, 
  beta1=1, 
  sigma=c(1, 1),
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)

#str(input)
stan_file <- stan_file_path
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_mid2 <- mod$sample(data = input, 
                               chains = 4,
                               parallel_chains = 4,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               adapt_delta = 0.95,
                               init=init_fun, 
                               max_treedepth = 10,
                               step_size = 1)
elapsed_time <- cmdstan_fit_mid2$time()
elapsed_time
elapsed_time$total/60


fit_summary_HS_mid2 <- cmdstan_fit_mid2$summary(variables = c("beta0", "beta1",  "ell[1]", 'ell[2]', 
                                                              'Astar[2,1]'), 
                                                c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_mid2




# D. Analyze with HSGP, lowest resolution -----

# change the window length if possible 
xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

# adjust the basis function to find proper number/ 22
#m1 <- initial_m[1]; m2 <- initial_m[2]; mstar <- m1*m2
m1<- 20; m2<- 20; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5, 1.5)
#c<- rep(initial_c, 2)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, 
                                  S1 = 1:m2)[,2:1]))
str(S)

mlscale<- max(dist(grid_summary_low[[1]]))/(2*2.75)

# adding the grid_area term,
input <- list(n = n3,
              y1 = Y_low, # outcome1: points
              y2 = M_low,
              d = 2, # dim=2
              mstar = mstar, 
              coords =  grid_summary_low[[1]],
              indices = S,
              L=L,
              log_grid_area=log(grid_area_low),
              is_centerted_PHI = 0)

# define initial values 
init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 1), 2, 2),
  beta0=2, 
  beta1=1, 
  sigma=c(1, 1),
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)

#str(input)
stan_file <- stan_file_path
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_low <- mod$sample(data = input, 
                              chains = 4,
                              parallel_chains = 4,
                              iter_warmup = 1000,
                              iter_sampling = 1000,
                              adapt_delta = 0.95,
                              init=init_fun, 
                              max_treedepth = 10,
                              step_size = 1)
elapsed_time <- cmdstan_fit_low$time()
elapsed_time
elapsed_time$total/60


fit_summary_HS_low <- cmdstan_fit_low$summary(variables = c("beta0", "beta1", "ell[1]", 'ell[2]', 
                                                            'Astar[2,1]'), 
                                              c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_low



# the summary of the fitting 
fit_summary_HS_high
fit_summary_HS_mid
fit_summary_HS_mid2
fit_summary_HS_low


# compare and save the criteria 
draws_df_high <- cmdstan_fit_high$draws(format = "df")
draws_df_mid  <- cmdstan_fit_mid$draws(format = "df")
draws_df_mid2 <- cmdstan_fit_mid2$draws(format = "df")
draws_df_low  <- cmdstan_fit_low$draws(format = "df")


# summary of the posterior: high resolution  
samples_beta0_HS<- draws_df_high[['beta0']]
samples_beta1_HS<- draws_df_high[['beta1']]
samples_rho1_HS<- draws_df_high[['ell[1]']]
samples_rho2_HS<- draws_df_high[['ell[2]']]
samples_sigma1_HS<- draws_df_high[['sigma[1]']]
samples_sigma2_HS<- draws_df_high[['sigma[2]']]
samples_corr_HS<- draws_df_high[['Astar[2,1]']]
# a11, a21 and a22
samples_a11_HS <- samples_sigma1_HS
samples_a21_HS <- samples_sigma2_HS*samples_corr_HS
samples_a22_HS <- samples_sigma2_HS*sqrt(1-samples_corr_HS^2)


# summary of the posterior: mid resolution 
samples_beta0_mid<- draws_df_mid[['beta0']]
samples_beta1_mid<- draws_df_mid[['beta1']]
samples_rho1_mid<- draws_df_mid[['ell[1]']]
samples_rho2_mid<- draws_df_mid[['ell[2]']]
samples_sigma1_mid<- draws_df_mid[['sigma[1]']]
samples_sigma2_mid<- draws_df_mid[['sigma[2]']]
samples_corr_mid<- draws_df_mid[['Astar[2,1]']]
# a11, a21 and a22
samples_a11_mid <- samples_sigma1_mid
samples_a21_mid <- samples_sigma2_mid*samples_corr_mid
samples_a22_mid <- samples_sigma2_mid*sqrt(1-samples_corr_mid^2)


# summary of the posterior: mid2 resolution 
samples_beta0_mid2<- draws_df_mid2[['beta0']]
samples_beta1_mid2<- draws_df_mid2[['beta1']]
samples_rho1_mid2<- draws_df_mid2[['ell[1]']]
samples_rho2_mid2<- draws_df_mid2[['ell[2]']]
samples_sigma1_mid2<- draws_df_mid2[['sigma[1]']]
samples_sigma2_mid2<- draws_df_mid2[['sigma[2]']]
samples_corr_mid2<- draws_df_mid2[['Astar[2,1]']]
# a11 a21 and a22
samples_a11_mid2 <- samples_sigma1_mid2
samples_a21_mid2 <- samples_sigma2_mid2*samples_corr_mid2
samples_a22_mid2 <- samples_sigma2_mid2*sqrt(1-samples_corr_mid2^2)


# summary of the posterior: low resolution 
samples_beta0_low<- draws_df_low[['beta0']]
samples_beta1_low<- draws_df_low[['beta1']]
samples_rho1_low<- draws_df_low[['ell[1]']]
samples_rho2_low<- draws_df_low[['ell[2]']]
samples_sigma1_low<- draws_df_low[['sigma[1]']]
samples_sigma2_low<- draws_df_low[['sigma[2]']]
samples_corr_low<- draws_df_low[['Astar[2,1]']]
# a11 a21 and a22
samples_a11_low <- samples_sigma1_low
samples_a21_low <- samples_sigma2_low*samples_corr_low
samples_a22_low <- samples_sigma2_low*sqrt(1-samples_corr_low^2)


# the extracted mean values 
df_pts_lambda<-as.data.frame(cbind(sim.pts$x, sim.pts$y, log_lambda1, log_lambda2))
colnames(df_pts_lambda)<- c('x','y','log_lambda1', 'log_lambda2')
# the differences between grid and point-level data 

# create vectors for the grids of w1
vectors_w1_high <- vector("list", n)  # create empty list
vectors_w1_mid <- vector("list", n1)  # create empty list
vectors_w1_mid2 <- vector("list", n2)  # create empty list
vectors_w1_low <- vector("list", n3)  # create empty list

vectors_w2_high <- vector("list", n)  # create empty list
vectors_w2_mid <- vector("list", n1)  # create empty list
vectors_w2_mid2 <- vector("list", n2)  # create empty list
vectors_w2_low <- vector("list", n3)  # create empty list


# create vectors for the fitted grids 
lambda1_fitted_high<- rep(NA, n)
lambda1_fitted_mid<- rep(NA, n1)
lambda1_fitted_mid2 <- rep(NA, n2)
lambda1_fitted_low <- rep(NA, n3)

log_lambda1_fitted_high<- rep(NA, n)
log_lambda1_fitted_mid<- rep(NA, n1)
log_lambda1_fitted_mid2 <- rep(NA, n2)
log_lambda1_fitted_low <- rep(NA, n3)


lambda2_fitted_high<- rep(NA, n)
lambda2_fitted_mid<- rep(NA, n1)
lambda2_fitted_mid2 <- rep(NA, n2)
lambda2_fitted_low <- rep(NA, n3)

log_lambda2_fitted_high<- rep(NA, n)
log_lambda2_fitted_mid<- rep(NA, n1)
log_lambda2_fitted_mid2 <- rep(NA, n2)
log_lambda2_fitted_low <- rep(NA, n3)




# high res
for (i in 1:n) {
  vectors_w1_high[[i]] <- draws_df_high[[paste0("w1[", i, "]")]]  # extract by name
  vectors_w2_high[[i]] <- draws_df_high[[paste0("w2[", i, "]")]]  # extract by name
  # fitted lambda1
  lambda1_fitted_high[i] <- mean( exp(samples_beta0_HS + samples_a11_HS*vectors_w1_high[[i]]) )
  # fitted lambda2
  lambda2_fitted_high[i] <-  mean( exp(samples_beta1_HS + samples_a21_HS*vectors_w1_high[[i]] + samples_a22_HS*vectors_w2_high[[i]] )*(Y_high[i]+0.05) )
  # fitted log_lambda1
  log_lambda1_fitted_high[i] <- mean( samples_beta0_HS + samples_a11_HS*vectors_w1_high[[i]]) 
  # fitted log_lambda2
  log_lambda2_fitted_high[i] <- mean( samples_beta1_HS + samples_a21_HS*vectors_w1_high[[i]] + samples_a22_HS*vectors_w2_high[[i]] + log(Y_high[i]+0.05) )
}



# mid res 
for (i in 1:n1) {
  vectors_w1_mid[[i]] <- draws_df_mid[[paste0("w1[", i, "]")]]  # extract by name
  vectors_w2_mid[[i]] <- draws_df_mid[[paste0("w2[", i, "]")]]  # extract by name
  
  # fitted lambda1
  lambda1_fitted_mid[i] <- mean(exp(samples_beta0_mid +  samples_a11_mid*vectors_w1_mid[[i]]))
  # fitted lambda2
  lambda2_fitted_mid[i] <-  mean( exp(samples_beta1_mid + samples_a21_mid*vectors_w1_mid[[i]] + samples_a22_mid*vectors_w2_mid[[i]] )*(Y_mid[i]+0.05) )
  
  # fitted log_lambda1
  log_lambda1_fitted_mid[i] <- mean( samples_beta0_mid +  samples_a11_mid*vectors_w1_mid[[i]]) 
  
  # fitted log_lambda2
  log_lambda2_fitted_mid[i] <- mean( samples_beta1_mid + samples_a21_mid*vectors_w1_mid[[i]] + samples_a22_mid*vectors_w2_mid[[i]] + log(Y_mid[i]+0.05) )
}




#  mid2 res 
for (i in 1:n2) {
  vectors_w1_mid2[[i]] <- draws_df_mid2[[paste0("w1[", i, "]") ]]  # extract by name
  vectors_w2_mid2[[i]] <- draws_df_mid2[[paste0("w2[", i, "]")]]  # extract by name
  
  # fitted lambda1
  lambda1_fitted_mid2[i] <- mean(exp(samples_beta0_mid2 + samples_a11_mid2*vectors_w1_mid2[[i]]))
  # fitted lambda2
  lambda2_fitted_mid2[i] <-  mean( exp(samples_beta1_mid2 + samples_a21_mid2*vectors_w1_mid2[[i]] + samples_a22_mid2*vectors_w2_mid2[[i]] )*(Y_mid2[i]+0.05) )
  
  # fitted log_lambda1
  log_lambda1_fitted_mid2[i] <- mean( samples_beta0_mid2 + samples_a11_mid2*vectors_w1_mid2[[i]]) 
  
  # fitted log_lambda2
  log_lambda2_fitted_mid2[i] <- mean( samples_beta1_mid2 + samples_a21_mid2*vectors_w1_mid2[[i]] + samples_a22_mid2*vectors_w2_mid2[[i]] + log(Y_mid2[i]+0.05) )
}



# low res 
for (i in 1:n3) {
  vectors_w1_low[[i]] <- draws_df_low[[paste0("w1[", i, "]")]]  # extract by name
  vectors_w2_low[[i]] <- draws_df_low[[paste0("w2[", i, "]")]]  # extract by name
  
  # fitted lambda1
  lambda1_fitted_low[i] <- mean(exp(samples_beta0_low + samples_a11_low*vectors_w1_low[[i]]))
  # fitted lambda2
  lambda2_fitted_low[i] <-  mean( exp(samples_beta1_low + samples_a21_low*vectors_w1_low[[i]] + samples_a22_low*vectors_w2_low[[i]] )*(Y_low[i]+0.05) )
  
  # fitted log_lambda1
  log_lambda1_fitted_low[i] <- mean( samples_beta0_low + samples_a11_low*vectors_w1_low[[i]]) 
  
  # fitted log_lambda2
  log_lambda2_fitted_low[i] <- mean( samples_beta1_low + samples_a21_low*vectors_w1_low[[i]] + samples_a22_low*vectors_w2_low[[i]] + log(Y_low[i]+0.05) )
}




# 3.1. Differences

# calculate the criteria 
mpp.df2<- data.frame(cbind(sim.pts$x, 
                           sim.pts$y,  
                           exp(log_lambda1), 
                           exp(log_lambda2)))
colnames(mpp.df2)<- c('x','y','lambda1', 'lambda2')


mpp.df2.log<- data.frame(cbind(sim.pts$x, 
                               sim.pts$y,  
                               log_lambda1, 
                               log_lambda2))
colnames(mpp.df2.log)<- c('x','y','loglambda1', 'loglambda2')



# 3. the fitted intensity on those points 
# fitted intensity and the grids: 
# high
fitted_res_high<- data.frame(cbind(grid_summary_large[[1]], 
                                   lambda1_fitted_high, 
                                   lambda2_fitted_high))
colnames(fitted_res_high)<- c('x', 'y', 'lambda1', 'lambda2')


fitted_res_log_high<- data.frame(cbind(grid_summary_large[[1]], 
                                       log_lambda1_fitted_high, 
                                       log_lambda2_fitted_high))
colnames(fitted_res_log_high)<- c('x', 'y', 'loglambda1', 'loglambda2')




# mid
fitted_res_mid<- data.frame(cbind(grid_summary_mid[[1]], 
                                  lambda1_fitted_mid, 
                                  lambda2_fitted_mid))
colnames(fitted_res_mid)<- c('x', 'y', 'lambda1', 'lambda2')


fitted_res_log_mid<- data.frame(cbind(grid_summary_mid[[1]], 
                                      log_lambda1_fitted_mid, 
                                      log_lambda2_fitted_mid))
colnames(fitted_res_log_mid)<- c('x', 'y', 'loglambda1', 'loglambda2')


# mid2
fitted_res_mid2<- data.frame(cbind(grid_summary_mid2[[1]], 
                                   lambda1_fitted_mid2, 
                                   lambda2_fitted_mid2))
colnames(fitted_res_mid2)<- c('x', 'y', 'lambda1', 'lambda2')

fitted_res_log_mid2<- data.frame(cbind(grid_summary_mid2[[1]], 
                                       log_lambda1_fitted_mid2, 
                                       log_lambda2_fitted_mid2))
colnames(fitted_res_log_mid2)<- c('x', 'y', 'loglambda1', 'loglambda2')




# low
fitted_res_low<- data.frame(cbind(grid_summary_low[[1]], 
                                  lambda1_fitted_low, 
                                  lambda2_fitted_low))
colnames(fitted_res_low)<- c('x', 'y', 'lambda1', 'lambda2')


fitted_res_log_low<- data.frame(cbind(grid_summary_low[[1]], 
                                      log_lambda1_fitted_low, 
                                      log_lambda2_fitted_low))
colnames(fitted_res_log_low)<- c('x', 'y', 'loglambda1', 'loglambda2')




# the breaks at different resolution 
x_breaks_h <- seq(-L1, L1, length.out= grid_size_large+1 )
y_breaks_h <- seq(-L2, L2, length.out= grid_size_large+1 )

x_breaks_m <- seq(-L1, L1, length.out= grid_size_mid+1 )
y_breaks_m <- seq(-L2, L2, length.out= grid_size_mid+1 )

x_breaks_m2 <- seq(-L1, L1, length.out= grid_size_mid2+1 )
y_breaks_m2 <- seq(-L2, L2, length.out= grid_size_mid2+1 )

x_breaks_l <- seq(-L1, L1, length.out= grid_size_low+1 )
y_breaks_l <- seq(-L2, L2, length.out= grid_size_low+1 )



cell_x_h <- findInterval(mpp.df2$x, x_breaks_h, rightmost.closed = TRUE)
cell_y_h <- findInterval(mpp.df2$y, y_breaks_h, rightmost.closed = TRUE)

cell_x_m <- findInterval(mpp.df2$x, x_breaks_m, rightmost.closed = TRUE)
cell_y_m <- findInterval(mpp.df2$y, y_breaks_m, rightmost.closed = TRUE)

cell_x_m2 <- findInterval(mpp.df2$x, x_breaks_m2, rightmost.closed = TRUE)
cell_y_m2 <- findInterval(mpp.df2$y, y_breaks_m2, rightmost.closed = TRUE)

cell_x_l <- findInterval(mpp.df2$x, x_breaks_l, rightmost.closed = TRUE)
cell_y_l <- findInterval(mpp.df2$y, y_breaks_l, rightmost.closed = TRUE)


cell_id_h  <- (cell_y_h - 1) * grid_size_large + cell_x_h
cell_id_m  <- (cell_y_m - 1) * grid_size_mid + cell_x_m
cell_id_m2 <- (cell_y_m2 - 1) * grid_size_mid2 + cell_x_m2
cell_id_l  <- (cell_y_l - 1) * grid_size_low + cell_x_l



# the formula: 
# build a function to compute Bias^2 and Variance for fitted values
compute_bias_var <- function(true_values, 
                             posterior_matrix, 
                             cell_ids) {
  # posterior_matrix: rows = posterior draws, cols = locations
  posterior_means <- colMeans(posterior_matrix)
  posterior_vars  <- apply(posterior_matrix, 2, var)
  
  
  # assign fitted means and variances to each point via its cell id
  fitted_means <- posterior_means[cell_ids]
  fitted_vars  <- posterior_vars[cell_ids]
  
  # compute Bias^2 and Variance averaged over points
  #bias2    <- mean((fitted_means - true_values)^2)
  #bias2    <- mean((fitted_means - true_values))^2
  bias2    <- (mean( (fitted_means - true_values) ))^2
  variance <- mean(fitted_vars)
  
  return(list(bias2 = bias2,
              variance = variance, 
              mse = bias2 + variance))
}

# Number of grid cells for each resolution
ncell_h <- max(cell_id_h)
ncell_m <- max(cell_id_m)
ncell_m2<- max(cell_id_m2)
ncell_l <- max(cell_id_l)

# 4 matrics: 
# Posterior matrices: log(lambda1)
posterior_mat_loglambda1_h <- sapply(1:ncell_h, function(j) {
  samples_beta0_HS + samples_a11_HS * draws_df_high[[paste0("w1[", j, "]")]]
}, simplify = "matrix")

posterior_mat_loglambda1_m <- sapply(1:ncell_m, function(j) {
  samples_beta0_mid + samples_a11_mid * draws_df_mid[[paste0("w1[", j, "]")]]
}, simplify = "matrix")

posterior_mat_loglambda1_m2 <- sapply(1:ncell_m2, function(j) {
  samples_beta0_mid2 + samples_a11_mid2 * draws_df_mid2[[paste0("w1[", j, "]")]]
}, simplify = "matrix")

posterior_mat_loglambda1_l <- sapply(1:ncell_l, function(j) {
  samples_beta0_low + samples_a11_low * draws_df_low[[paste0("w1[", j, "]")]]
}, simplify = "matrix")

# Posterior matrices: log(lambda2)
posterior_mat_loglambda2_h <- sapply(1:ncell_h, function(j) {
  samples_beta1_HS +
    samples_a21_HS * draws_df_high[[paste0("w1[", j, "]")]] +
    samples_a22_HS * draws_df_high[[paste0("w2[", j, "]")]] +
    log(Y_high[j] + 0.05)
}, simplify = "matrix")

posterior_mat_loglambda2_m <- sapply(1:ncell_m, function(j) {
  samples_beta1_mid +
    samples_a21_mid * draws_df_mid[[paste0("w1[", j, "]")]] +
    samples_a22_mid * draws_df_mid[[paste0("w2[", j, "]")]] +
    log(Y_mid[j] + 0.05)
}, simplify = "matrix")

posterior_mat_loglambda2_m2 <- sapply(1:ncell_m2, function(j) {
  samples_beta1_mid2 +
    samples_a21_mid2 * draws_df_mid2[[paste0("w1[", j, "]")]] +
    samples_a22_mid2 * draws_df_mid2[[paste0("w2[", j, "]")]] +
    log(Y_mid2[j] + 0.05)
}, simplify = "matrix")

posterior_mat_loglambda2_l <- sapply(1:ncell_l, function(j) {
  samples_beta1_low +
    samples_a21_low * draws_df_low[[paste0("w1[", j, "]")]] +
    samples_a22_low * draws_df_low[[paste0("w2[", j, "]")]] +
    log(Y_low[j] + 0.05)
}, simplify = "matrix")



# Now compute Bias^2, Variance, MSE
bias_var_lambda1_h <- compute_bias_var(mpp.df2.log$loglambda1, posterior_mat_loglambda1_h, cell_ids=cell_id_h)
bias_var_lambda1_m <- compute_bias_var(mpp.df2.log$loglambda1, posterior_mat_loglambda1_m, cell_ids=cell_id_m)
bias_var_lambda1_m2 <- compute_bias_var(mpp.df2.log$loglambda1, posterior_mat_loglambda1_m2, cell_ids=cell_id_m2)
bias_var_lambda1_l <- compute_bias_var(mpp.df2.log$loglambda1, posterior_mat_loglambda1_l, cell_ids=cell_id_l)

bias_var_lambda2_h <- compute_bias_var(mpp.df2.log$loglambda2, posterior_mat_loglambda2_h, cell_ids=cell_id_h)
bias_var_lambda2_m <- compute_bias_var(mpp.df2.log$loglambda2, posterior_mat_loglambda2_m, cell_ids=cell_id_m)
bias_var_lambda2_m2 <- compute_bias_var(mpp.df2.log$loglambda2, posterior_mat_loglambda2_m2, cell_ids=cell_id_m2)
bias_var_lambda2_l <- compute_bias_var(mpp.df2.log$loglambda2, posterior_mat_loglambda2_l, cell_ids=cell_id_l)


bias_var_lambda1_h$bias2
bias_var_lambda1_m$bias2
bias_var_lambda1_m2$bias2
bias_var_lambda1_l$bias2


bias_var_lambda2_h$bias2
bias_var_lambda2_m$bias2
bias_var_lambda2_m2$bias2
bias_var_lambda2_l$bias2


# Bias^2: log(lambda1)
#diff_pts_log_high <- (1/sim.pts$n)*sum( (fitted_res_log_high[cell_id_h,]$loglambda1 - mpp.df2.log$loglambda1)^2) 
#diff_pts_log_mid <- (1/sim.pts$n)*sum( (fitted_res_log_mid[cell_id_m,]$loglambda1 - mpp.df2.log$loglambda1)^2) 
#diff_pts_log_mid2 <- (1/sim.pts$n)*sum( (fitted_res_log_mid2[cell_id_m2,]$loglambda1 - mpp.df2.log$loglambda1)^2) 
#diff_pts_log_low<- (1/sim.pts$n)*sum( (fitted_res_log_low[cell_id_l,]$loglambda1 - mpp.df2.log$loglambda1)^2) 


# log(lambda2)
#diff_marks_log_high<- (1/sim.pts$n)*sum( (fitted_res_log_high[cell_id_h,]$loglambda2 - mpp.df2.log$loglambda2)^2) 
#diff_marks_log_mid<- (1/sim.pts$n)*sum( (fitted_res_log_mid[cell_id_m,]$loglambda2 - mpp.df2.log$loglambda2)^2) 
#diff_marks_log_mid2<- (1/sim.pts$n)*sum( (fitted_res_log_mid2[cell_id_m2,]$loglambda2 - mpp.df2.log$loglambda2)^2) 
#diff_marks_log_low<- (1/sim.pts$n)*sum( (fitted_res_log_low[cell_id_l,]$loglambda2 - mpp.df2.log$loglambda2)^2) 



# visualization for the lambda difference between the grids: 
# 1. create a dataframe containing all the differences between the true lambda and the fitted lambda 
viz.lambda.diff.h<- data.frame(mpp.df2$x, 
                               mpp.df2$y,
                               (fitted_res_log_high[cell_id_h,]$loglambda1 - mpp.df2.log$loglambda1),
                               (fitted_res_log_high[cell_id_h,]$loglambda2 - mpp.df2.log$loglambda2))
colnames(viz.lambda.diff.h) <- c('x','y',
                                 'loglambda1_diff',
                                 'loglambda2_diff') 

viz.lambda.diff.m<- data.frame(mpp.df2$x, mpp.df2$y,
                               (fitted_res_log_mid[cell_id_m,]$loglambda1 - mpp.df2.log$loglambda1), 
                               (fitted_res_log_mid[cell_id_m,]$loglambda2 - mpp.df2.log$loglambda2))
colnames(viz.lambda.diff.m) <- c('x','y',
                                 'loglambda1_diff',
                                 'loglambda2_diff') 

viz.lambda.diff.m2<- data.frame(mpp.df2$x, mpp.df2$y,
                                (fitted_res_log_mid2[cell_id_m2,]$loglambda1 - mpp.df2.log$loglambda1),
                                (fitted_res_log_mid2[cell_id_m2,]$loglambda2 - mpp.df2.log$loglambda2))
colnames(viz.lambda.diff.m2) <- c('x','y',
                                  'loglambda1_diff',
                                  'loglambda2_diff') 

viz.lambda.diff.l<- data.frame(mpp.df2$x, mpp.df2$y,
                               (fitted_res_log_low[cell_id_l,]$loglambda1 - mpp.df2.log$loglambda1), 
                               (fitted_res_log_low[cell_id_l,]$loglambda2 - mpp.df2.log$loglambda2))
colnames(viz.lambda.diff.l) <- c('x','y',
                                 'loglambda1_diff',
                                 'loglambda2_diff') 



# viz the differences of the two
l1_h<- ggplot(data=viz.lambda.diff.h)+
  #scale_color_gradient(low='snow2', high='red')+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda1_diff), size=2)+
  theme_classic()

l1_m<- ggplot(data=viz.lambda.diff.m)+
  #scale_color_gradient(low='snow2', high='red')+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda1_diff), size=2)+
  theme_classic()

l1_m2<- ggplot(data=viz.lambda.diff.m2)+
  #scale_color_gradient(low='snow2', high='red')+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda1_diff), size=2)+
  theme_classic()

l1_l<- ggplot(data=viz.lambda.diff.l)+
  #scale_color_gradient(low='snow2', high='red')+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda1_diff), size=2)+
  theme_classic()


combined_plot_l1<- plot_grid(l1_h, l1_m, 
                             l1_m2, l1_l, nrow=1)


ggsave(paste0('/lustre09/project/6003552/mingchi/Simulation-Two/lambda1_diff_panels_',
              as.character(iter),'.png'), 
       combined_plot_l1, width = 16, height = 4, dpi = 300)

# lambda2:
l2_h<- ggplot(data=viz.lambda.diff.h)+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda2_diff), size=2)+
  theme_classic()

l2_m<- ggplot(data=viz.lambda.diff.m)+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda2_diff), size=2)+
  theme_classic()

l2_m2<- ggplot(data=viz.lambda.diff.m2)+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda2_diff), size=2)+
  theme_classic()

l2_l<- ggplot(data=viz.lambda.diff.l)+
  scale_color_viridis_c()+
  geom_point(aes(x=x, y=y, col=loglambda2_diff), size=2)+
  theme_classic()

combined_plot_l2<- plot_grid(l2_h, l2_m, 
                             l2_m2, l2_l, nrow=1)


ggsave(paste0('/lustre09/project/6003552/mingchi/Simulation-Two/lambda2_diff_panels_',
              as.character(iter),'.png'), combined_plot_l2, width = 16, height = 4, dpi = 300)








# the visualization for the differences: add red line: y=x 
# fitted log-lambda1 vs. the true log-lambda1
png(paste0('/lustre09/project/6003552/mingchi/Simulation-Two/log_lambda1_',
           as.character(iter),'.png'),
    width=1200, height=400)




par(mfrow=c(1,4))
plot(fitted_res_log_high[cell_id_h,]$loglambda1,
     mpp.df2.log$loglambda1,
     main="50*50 grids", 
     xlab='Fitted log-lambda1',
     ylab='True log-lambda1', 
     pch=16, 
     cex=1,
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_mid[cell_id_m,]$loglambda1,
     mpp.df2.log$loglambda1,
     main="40*40 grids", 
     xlab='Fitted log-lambda1',
     ylab='True log-lambda1', 
     pch=16, 
     cex=1, 
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_mid2[cell_id_m2,]$loglambda1,
     mpp.df2.log$loglambda1,
     main="30*30 grids",
     xlab='Fitted log-lambda1',
     ylab='True log-lambda1', 
     pch=16, 
     cex=1,
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_low[cell_id_l,]$loglambda1,
     mpp.df2.log$loglambda1,
     main="20*20 grids", 
     xlab='Fitted log-lambda1',
     ylab='True log-lambda1', 
     pch=16, 
     cex=1, 
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
par(mfrow=c(1,1))
dev.off()




# fitted log-lambda2 vs. the true log-lambda2
png(paste0('/lustre09/project/6003552/mingchi/Simulation-Two/log_lambda2_',
           as.character(iter),'.png'),width=1200, height=400)



par(mfrow=c(1,4))
plot(fitted_res_log_high[cell_id_h,]$loglambda2,
     mpp.df2.log$loglambda2,
     main="50*50 grids", 
     xlab='Fitted log-lambda2',
     ylab='True log-lambda2', 
     pch=16, 
     cex=1, 
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_mid[cell_id_m,]$loglambda2,
     mpp.df2.log$loglambda2,
     main="40*40 grids",
     xlab='Fitted log-lambda2',
     ylab='True log-lambda2', 
     pch=16, 
     cex=1,
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_mid2[cell_id_m2,]$loglambda2,
     mpp.df2.log$loglambda2,
     main="30*30 grids", 
     xlab='Fitted log-lambda2',
     ylab='True log-lambda2', 
     pch=16, 
     cex=1,
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
plot(fitted_res_log_low[cell_id_l,]$loglambda2,
     mpp.df2.log$loglambda2,
     main="20*20 grids", 
     xlab='Fitted log-lambda2',
     ylab='True log-lambda2', 
     pch=16, 
     cex=1,
     cex.main=1.5,        
     cex.lab=1.5,           
     cex.axis=1.5)
abline(0,1, col="red", lwd=2, lty=2)
par(mfrow=c(1,1))
dev.off()










# save the data into the data_frame 

#save_mpp_val<-cbind(c(diff_pts_high, diff_pts_mid, diff_pts_mid2, diff_pts_low ), 
#                    c(diff_marks_high, diff_marks_mid, diff_marks_mid2, diff_marks_low ), 
#                    c(diff_pts_log_high, diff_pts_log_mid, diff_pts_log_mid2, diff_pts_log_low ), 
#                    c(diff_marks_log_high, diff_marks_log_mid, diff_marks_log_mid2, diff_marks_log_low ) )
#colnames(save_mpp_val)<- c('lambda1', 'lambda2', 'loglambda1', 'loglambda2')

save_mpp_val<-cbind(c(bias_var_lambda1_h$bias2 , 
                      bias_var_lambda1_m$bias2, 
                      bias_var_lambda1_m2$bias2, 
                      bias_var_lambda1_l$bias2 ), 
                    c(bias_var_lambda2_h$bias2 , 
                      bias_var_lambda2_m$bias2, 
                      bias_var_lambda2_m2$bias2, 
                      bias_var_lambda2_l$bias2), 
                    
                    c(bias_var_lambda1_h$variance , 
                      bias_var_lambda1_m$variance, 
                      bias_var_lambda1_m2$variance, 
                      bias_var_lambda1_l$variance ), 
                    c(bias_var_lambda2_h$variance , 
                      bias_var_lambda2_m$variance, 
                      bias_var_lambda2_m2$variance, 
                      bias_var_lambda2_l$variance),
                    
                    c(bias_var_lambda1_h$mse, 
                      bias_var_lambda1_m$mse, 
                      bias_var_lambda1_m2$mse, 
                      bias_var_lambda1_l$mse ), 
                    c(bias_var_lambda2_h$mse , 
                      bias_var_lambda2_m$mse, 
                      bias_var_lambda2_m2$mse, 
                      bias_var_lambda2_l$mse))
colnames(save_mpp_val)<- c('bias2-loglambda1', 
                           'bias2-loglambda2', 
                           'var-loglambda1', 
                           'var-loglambda2', 
                           'mse-loglambda1', 
                           'mse-loglambda2')

save_mpp_val



# write the csv files
# let's first check this in the local PC 
# then resume back to the cluster 
write.csv(save_mpp_val, paste0('/lustre09/project/6003552/mingchi/Simulation-Two/summary_grid_M2_compare_',
                               as.character(iter),'.csv'))



