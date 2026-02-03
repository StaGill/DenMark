###################################################################################
# This session summarize the simulation that need to be written in the document   #
###################################################################################
# Revised S1: 


# Simulation 1: exact GP and HSGP 

# 0. input the parameter 
args=commandArgs(trailingOnly = TRUE)

# 1. generate a GP ----

###################################################
# Poisson-Poisson Hilbert space GP approximation 
##################################################


#rm(list=ls())
#graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(bayesplot)

# 
if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}


# the functions: 
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


# the iter parameters
iter <- as.numeric(args[1])
print(iter)

set.seed(iter)

options(cmdstanr_max_rows = 20)



######################################################################
# Generating the data
######################################################################
grid_res <- 22
coords <- unname(as.matrix(expand.grid(
  x = seq(-0.99, 0.99, length.out = grid_res),
  y = seq(-0.99, 0.99, length.out = grid_res)
)))
n <- nrow(coords); n
# distance matrix 
dist_matrix <- as.matrix(dist(coords))
grid_area <- ((0.99 + 0.99) / grid_res)^2

sigma <- c(1, 1)
R <- cbind(c(1, 0.5),
           c(0.5, 1))
Astar <- t(chol(R))
A <- diag(sigma) %*% Astar

lscale <- c(0.5, 0.3)
distMat <- fields::rdist(coords)

matern32 <- function(d, sigma, lscale){
  ds <- sqrt(3) * d / lscale
  (sigma^2) * (1 + ds) * exp(-ds)
}

# generate 2 independent Gaussian Processes (w columns)
w <- sapply(1:2, function(i) {
  K <- matern32(d = distMat, sigma = 1, lscale = lscale[i])
  K <- K + diag(x = 1e-9, nrow = n, ncol = n)
  cholK <- chol(K)
  drop(crossprod(cholK, rnorm(n)))
})
# w is n x 2 (each column is one GP)
w <- as.matrix(w)

# combine through A to create latent effects z (n x 2)
# note: in your original code w was transposed in places — keep consistent:
z_mat <- (w %*% t(A))    # n x 2, first col -> z1, second -> z2
z1 <- z_mat[, 1]
z2 <- z_mat[, 2]

#####################################
# Assign parameter values 
#####################################

beta_pts   <- 1      # intercept for points (eta1)
beta_marks <- -1     # intercept for marks (eta2)
offset     <- 0.05   # same small offset used in Stan


# generate points 
mu1 <- exp(beta_pts + z1)
y1 <- rpois(n = n, lambda = grid_area*mu1)
#hist(y1, main = "Simulated points y1")

# generate marks 
mu2 <- exp(beta_marks +  log(y1 + offset) + z2)
y2 <- rpois(n = n, lambda = mu2)
#hist(y2, main = "Simulated marks y2")



# plot the point pattern and mark pattern 



df1 <- data.frame(cbind(coords, y1, y2))
colnames(df1) <- c('x', 'y', 'points', 'marks')




# 2. Fit with exact GP: Poisson-Poisson model -----

## Lastly,  prepare the stan file 
stan_file<- '/lustre09/project/6003552/mingchi/git/Exact-GP-v2.stan'
# stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/Exact-GP-v2.stan'

set_cmdstan_path(cmdstan_path())

# Step 2: Compile and Run the model

# Define fixed parameter values for simulation
stan_data <- list(
  N = n,
  Y = y1,
  M = y2,
  dist_matrix = dist_matrix,
  grid_area=grid_area
)

# define initial values 
init_fun <- function()list(
  beta0=beta_pts, 
  beta1=beta_marks, 
  #beta_y=beta_y,
  sigma=c(1,1), 
  L=matrix(c(1, 0, 0.5, 0.866), 2, 2),
  rho1=lscale[1], 
  rho2=lscale[2],
  z1=rep(0,n),
  z2=rep(0,n)
)

# Compile and fit the Stan model

mod<- cmdstan_model(stan_file, compile=T)
mod$check_syntax(pedantic=TRUE)

# default iter=1000 
cmdstan_fit_exact<- mod$sample(data=stan_data,
                               chains=4,
                               parallel_chains=4,
                               refresh=200,
                               init=init_fun, 
                               adapt_delta = 0.95,
                               max_treedepth = 10,
                               step_size = 1)

elapsed_time_exact <- cmdstan_fit_exact$time()
elapsed_time_exact
elapsed_time_exact$total/60

# summary for only the parameters 
post_summary_exact<-cmdstan_fit_exact$summary(variables = c('beta0','beta1','A[1,1]','A[2,1]','A[2,2]',
                                                            'sigma[1]', 'sigma[2]' , 'rho1','rho2', 'corr'),
                                              c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
post_summary_exact


#### compare with the true parameters ####
beta_pts
beta_marks
A
sigma
lscale
##########################################




# 3. Fit with HSGP approximation, 30*30 basis -----


#################################################################
# Preparing for Hilbert Space Approximate GP
#################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)

m1 <- 30; m2 <- 30; mstar <- m1*m2
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
input <- list(n = n,
              y1 = y1, # outcome1: points 
              y2 = y2, # outcome2: marks 
              log_grid_area=log(grid_area),
              d = 2, # dim=2
              L = L,
              mstar = mstar, 
              coords = coords,
              indices = S,
              is_centerted_PHI = 0)


init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 0.866), 2, 2),
  sigma=c(1, 1),
  beta0=beta_pts, 
  beta1=beta_marks, 
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)



str(input)

stan_file <- '/lustre09/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-v2.stan'

mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_HS_30 <- mod$sample(data = input, 
                                chains = 4,
                                parallel_chains = 4,
                                iter_warmup = 1000,
                                iter_sampling = 1000,
                                adapt_delta = 0.95,
                                max_treedepth = 10,
                                init=init_fun,
                                step_size = 1)
elapsed_time_HS_30 <- cmdstan_fit_HS_30$time()
elapsed_time_HS_30
elapsed_time_HS_30$total/60
#cmdstan_fit_HS$cmdstan_diagnose()


fit_summary_HS_30 <- cmdstan_fit_HS_30$summary(variables = c('beta0','beta1','Astar[2,1]','sigma[1]', 'sigma[2]' , 
                                                             'ell[1]','ell[2]'), 
                                               c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_30


# save the exact GP
save_df_exact<- as.data.frame(cbind(post_summary_exact, elapsed_time_exact$total))

#write.csv(save_df_exact, 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/exact_summary.csv')


write.csv(save_df_exact, 
          paste0('/lustre09/project/6003552/mingchi/Simulation-1/exact_GP_',
                 iter, '.csv'))


# save the HSGP: 30*30
save_df_HS_30 <- as.data.frame(cbind(fit_summary_HS_30, elapsed_time_HS_30$total))

write.csv(save_df_HS_30, 
          paste0('/lustre09/project/6003552/mingchi/Simulation-1/HSGP_approx_30_',
                 iter, '.csv'))

draws_df_exact <- cmdstan_fit_exact$draws(format = "df")
draws_df_HS_30 <- cmdstan_fit_HS_30$draws(format = "df")


#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP_30 <- mcmc_trace(draws_df_HS_30,  
                        pars = c("beta0", "beta1", "Astar[2,1]", "sigma[1]", "sigma[2]", "ell[1]", "ell[2]", "w1[1]"), 
                        facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines

p_exact <- mcmc_trace(draws_df_exact,  
                      pars = c("beta0", "beta1", "A[2,1]", 'corr',"sigma[1]", "sigma[2]", "rho1", "rho2", "omega1[1]"), 
                      facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines



ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_HS_30_',
                   iter, '.png'),
  plot = p_HSGP_30,
  path = "/lustre09/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)


ggsave(
  #filename = "trace_exact.png",  # or "traceplot.pdf"
  filename= paste0('trace_exact_',
                   iter, '.png'),
  plot = p_exact,
  path = "/lustre09/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)
