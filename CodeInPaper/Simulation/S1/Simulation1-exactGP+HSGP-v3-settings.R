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


#iter<- 7
#set.seed(iter)
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
#mu2 <- exp(beta_marks +  log(y1 + offset) + z2)
#y2 <- rpois(n = n, lambda = mu2)
#hist(y2, main = "Simulated marks y2")

mu2 <- exp(beta_marks + z2)
y2 <- rpois(n = n, lambda = y1*mu2)




# plot the point pattern and mark pattern 

#df1 <- data.frame(cbind(coords, y1, y2))
#colnames(df1) <- c('x', 'y', 'points', 'marks')
#par(mfrow=c(1,2))
#plot(df1$x, df1$y, cex=df1$points+0.05)
#plot(df1$x, df1$y, cex=df1$marks+0.05)
#par(mfrow=c(1,1))
#write.csv(df1, 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/Figure-S1-simulated-data.csv')



# 3. Fit with HSGP approximation, 30*30 basis -----


#################################################################
# Preparing for Hilbert Space Approximate GP
#################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)

m1 <- 25; m2 <- 25; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.02, 1.02)
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

stan_file <- '/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-vtest.stan'



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

# save one examples into samples for drawing: 
#cmdstan_fit_HS_30$save_object(
#  "C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/S1-M2-cmdstan_fit.rds"
#)



# save the HSGP: 30*30
save_df_HS_30 <- as.data.frame(cbind(fit_summary_HS_30, elapsed_time_HS_30$total))

write.csv(save_df_HS_30, 
          paste0('/project/6003552/mingchi/Simulation-1/HSGP_approx_30_',
                 iter, '.csv'))



draws_df_HS_30 <- cmdstan_fit_HS_30$draws(format = "df")
# 




#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP_30 <- mcmc_trace(draws_df_HS_30,  
                        pars = c("beta0", "beta1", "Astar[2,1]", "sigma[1]", "sigma[2]", "ell[1]", "ell[2]", "w1[1]"), 
                        facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines


ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_HS_30_',
                   iter, '.png'),
  plot = p_HSGP_30,
  path = "/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)



# 25*25: 
#################################################################
# Preparing for Hilbert Space Approximate GP
#################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)

m1 <- 10; m2 <- 10; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.02, 1.02)
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
#stan_file <- '/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file <- '/lustre09/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-v2.stan'

mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_HS_25 <- mod$sample(data = input, 
                                chains = 4,
                                parallel_chains = 4,
                                iter_warmup = 1000,
                                iter_sampling = 1000,
                                adapt_delta = 0.95,
                                max_treedepth = 10,
                                init=init_fun,
                                step_size = 1)
elapsed_time_HS_25 <- cmdstan_fit_HS_25$time()
elapsed_time_HS_25
elapsed_time_HS_25$total/60
#cmdstan_fit_HS$cmdstan_diagnose()


fit_summary_HS_25 <- cmdstan_fit_HS_25$summary(variables = c('beta0','beta1','Astar[2,1]','sigma[1]', 'sigma[2]' , 
                                                             'ell[1]','ell[2]'), 
                                               c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_25


# save the HSGP: 30*30
save_df_HS_25 <- as.data.frame(cbind(fit_summary_HS_25, elapsed_time_HS_25$total))

write.csv(save_df_HS_25, 
          paste0('/project/6003552/mingchi/Simulation-1/HSGP_approx_25_',
                 iter, '.csv'))



draws_df_HS_25 <- cmdstan_fit_HS_25$draws(format = "df")


#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP_25 <- mcmc_trace(draws_df_HS_25,  
                        pars = c("beta0", "beta1", "Astar[2,1]", "sigma[1]", "sigma[2]", "ell[1]", "ell[2]", "w1[1]"), 
                        facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines


ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_HS_25_',
                   iter, '.png'),
  plot = p_HSGP_25,
  path = "/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)



#################################################################
# Preparing for Hilbert Space Approximate GP
#################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)

m1 <- 5; m2 <- 5; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.02, 1.02)
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

#stan_file <- '/lustre09/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-v2.stan'

mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_HS_5 <- mod$sample(data = input, 
                                chains = 4,
                                parallel_chains = 4,
                                iter_warmup = 1000,
                                iter_sampling = 1000,
                                adapt_delta = 0.95,
                                max_treedepth = 10,
                                init=init_fun,
                                step_size = 1)
elapsed_time_HS_5 <- cmdstan_fit_HS_5$time()
elapsed_time_HS_5
elapsed_time_HS_5$total/60
#cmdstan_fit_HS$cmdstan_diagnose()


fit_summary_HS_5 <- cmdstan_fit_HS_5$summary(variables = c('beta0','beta1','Astar[2,1]','sigma[1]', 'sigma[2]' , 
                                                             'ell[1]','ell[2]'), 
                                               c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_5


# save the HSGP: 30*30
save_df_HS_5 <- as.data.frame(cbind(fit_summary_HS_5, elapsed_time_HS_5$total))

write.csv(save_df_HS_5, 
          paste0('/project/6003552/mingchi/Simulation-1/HSGP_approx_5_',
                 iter, '.csv'))



draws_df_HS_5 <- cmdstan_fit_HS_5$draws(format = "df")


#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP_5 <- mcmc_trace(draws_df_HS_5,  
                        pars = c("beta0", "beta1", "Astar[2,1]", "sigma[1]", "sigma[2]", "ell[1]", "ell[2]", "w1[1]"), 
                        facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines


ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_HS_5_',
                   iter, '.png'),
  plot = p_HSGP_5,
  path = "/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)




#################################################################
# Preparing for Hilbert Space Approximate GP
#################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)

m1 <- 2; m2 <- 2; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.02, 1.02)
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

#stan_file <- '/lustre09/project/6003552/mingchi/git/S1-HSGP-v2.stan'
#stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/S1-HSGP-v2.stan'

mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit_HS_3 <- mod$sample(data = input, 
                                chains = 4,
                                parallel_chains = 4,
                                iter_warmup = 1000,
                                iter_sampling = 1000,
                                adapt_delta = 0.95,
                                max_treedepth = 10,
                                init=init_fun,
                                step_size = 1)
elapsed_time_HS_3 <- cmdstan_fit_HS_3$time()
elapsed_time_HS_3
elapsed_time_HS_3$total/60
#cmdstan_fit_HS$cmdstan_diagnose()


fit_summary_HS_3 <- cmdstan_fit_HS_3$summary(variables = c('beta0','beta1','Astar[2,1]','sigma[1]', 'sigma[2]' , 
                                                             'ell[1]','ell[2]'), 
                                               c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS_3


# save the HSGP: 30*30
save_df_HS_3 <- as.data.frame(cbind(fit_summary_HS_3, elapsed_time_HS_3$total))

write.csv(save_df_HS_3, 
          paste0('/project/6003552/mingchi/Simulation-1/HSGP_approx_3_',
                 iter, '.csv'))



draws_df_HS_3 <- cmdstan_fit_HS_3$draws(format = "df")


#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP_3 <- mcmc_trace(draws_df_HS_3,  
                        pars = c("beta0", "beta1", "Astar[2,1]", "sigma[1]", "sigma[2]", "ell[1]", "ell[2]", "w1[1]"), 
                        facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines


ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_HS_3_',
                   iter, '.png'),
  plot = p_HSGP_3,
  path = "/project/6003552/mingchi/Simulation-1",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)

  
