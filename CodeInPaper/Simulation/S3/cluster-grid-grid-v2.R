###################################################################################
# This session summarize the simulation that need to be written in the document   #
###################################################################################

# Simulation 3: Grid-level data fitted with the grid-level model  

# simulate from the grid-level data and fit with the grid-level data 
# just for the most general model case 
# This is to start from M2: -----
# Simulate from M2 and fit with M2: 

args=commandArgs(trailingOnly = TRUE)

library(MASS)
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(rethinking)
library(extraDistr)
library(bayesplot)

#library(MCMCpack)
# plot the point pattern and mark pattern 
#library('plot.matrix')

if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}

iter <- as.numeric(args[1])
print(iter)




# Rejection sampling for th zero-truncated normal distribution 
#ztnormdist<- function(n, mean, sd, a){
#  samplesztn<- numeric(n)
#  i<- 1
#  while (i <= n) {
#    x <- rnorm(1, mean, sd)
#    if (x >= a ) {
#      samplesztn[i] <- x
#      i <- i + 1
#    }
#  }
#  return(samplesztn)
#}

# the functions: 
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)



# Fit with GP ==========
set.seed(iter)
grid_res <- 30
grid_size<- grid_res


L1<- 4.4
L2<- 3.6

coords <- unname(as.matrix(expand.grid(
  x = seq(-L1, L1, length.out = grid_res),
  y = seq(-L2, L2, length.out = grid_res)
)))
n <- nrow(coords); n

# distance matrix 
dist_matrix <- as.matrix(dist(coords))
#grid_area<- ((0.99+0.99)/grid_res )^2 # try to remove this 

grid_area <- ((L1 + L1) / grid_res)*((L2 + L2) / grid_res)



#================ parameters; ============== 
sigma <- c(1, 1)
R <- cbind(c(1,0.5),
           c(0.5,1))
beta_pts<- 2
beta_marks<- 1
lscale <- c(2, 1) # 0.8, 0.4
offset<- 0.05
#==============================================

Astar <- t(chol(R))
A <- diag(sigma) %*% Astar


distMat <- fields::rdist(coords)

matern32 <- function(d, sigma, lscale){
  ds <- sqrt(3)*d/lscale
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
z_mat <- (w %*% t(A))    # n x 2, first col -> z1, second -> z2
z1 <- z_mat[, 1]
z2 <- z_mat[, 2]

# 1.5 Calculate eta1 and eta2 for points and marks  
# grid_area
log_grid_area<- log(grid_area)


# generate points 
mu1 <- exp(beta_pts + z1)
y1 <- rpois(n = n, lambda = grid_area*mu1)
Y<- y1
#hist(y1, main = "Simulated points y1")

# generate marks 
mu2 <- exp(beta_marks + log(y1 + offset) + z2)
y2 <- rpois(n = n, lambda = mu2)
M<- y2
#hist(y2, main = "Simulated marks y2")




# Combine into a data frame for easy handling
simulated_data <- data.frame(coords, Y = Y, M = M)

#write.csv(simulated_data, file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/S3-sim-YM.csv')

write.csv(simulated_data, 
          paste0('/project/6003552/mingchi/Simulations-Three/S3-simulated_data_',
                 iter, '.csv'))



true.corr<- cor(Y, M/(Y+0.05))


# rough correlation 

# get a prior for the mean.range 

distMat <- fields::rdist(coords)
mrange<-max(distMat)/(2*2.75)



# 2. HSGP: use 22 basis functions -----

################################################
# HSGP approximation to one-field model -----
################################################

xRangeDat <- c(-L1,L1)
yRangeDat <- c(-L2,L2)

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

# define initial values 

init_fun <- function()list(
  Astar=matrix(c(1, 0, 
                 0.5, 0.866), 2, 2),
  sigma=c(1, 1),
  beta0=beta_pts, 
  beta1=beta_marks, 
  #beta_y=beta_y,
  ell=c(lscale[1], lscale[2]),
  betab1=rep(0, m1*m2),
  betab2=rep(0, m1*m2)
)



str(input)


#stan_file <- '/project/6003552/mingchi/git/S3-HSGP-v2.stan'
stan_file<- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation3/revised-grid-grid-simulation/S3-HSGP-v2.stan'
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



fit_summary_HS <- cmdstan_fit_HS$summary(variables = c("beta0", "beta1",  "ell[1]", 'ell[2]', 
                                                       'Astar[2,1]', 'sigma[1]', 'sigma[2]'), 
                                         c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_HS

# save the summary results ----
#draws_df_HS <- cmdstan_fit_HS$draws(format = "df")
save_df_HS<- as.data.frame(cbind(fit_summary_HS, elapsed_time$total))

write.csv(save_df_HS, 
          paste0('/project/6003552/mingchi/Simulations-Three/S3-HSGP_',
                 iter, '.csv'))


# save the Aqp4 as the file: 
#saveRDS(cmdstan_fit_HS, 
#        file = 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/S3-M2-cmdstan_fit.rds')


# check the identifiable issues with the elements in A matrix 
#varscheck<- c("sigma[1]","sigma[2]","Astar[1,1]","Astar[2,1]","Astar[2,2]")
#draws_df_HS <- cmdstan_fit_HS$draws(format = "df", variables = varscheck)
#means<-colMeans(draws_df_HS)
#sigma_mean <- c(means["sigma[1]"], means["sigma[2]"])
#Astar_mean <- matrix(c(means["Astar[1,1]"], 0,
 #                      means["Astar[2,1]"], means["Astar[2,2]"]), 2, 2, byrow=TRUE)
#A_hat <- diag(sigma_mean) %*% Astar_mean
#
# A_true is from your sim: A <- diag(sigma) %*% Astar
# check A_hat and A_true: 
#print("Estimated A_hat:")
#print(A_hat)
#print("True A:")
#print(A)
#A_mat_compare <- as.data.frame(cbind(as.vector(A_hat), as.vector(A), true.corr))
#colnames(A_mat_compare)<- c('Estimated','True', 'True-corr-Y-MY')
#rownames(A_mat_compare)<- c('A11', 'A21', 'A12','A22')
#
#write.csv(A_mat_compare, 
#          file=paste0('/project/6003552/mingchi/Simulations-Three/A_mat_compare_',iter,'.csv'))

# Save the traceplots:

## Posterior draws
draws_df_HS <- cmdstan_fit_HS$draws(format = "df")
#draws_df_HS

#library(bayesplot)
color_scheme_set("brewer-Spectral")
p_HSGP <- mcmc_trace(draws_df_HS,  
                     pars = c("beta0", "beta1",  "ell[1]", 'ell[2]', 
                              'Astar[2,1]', 'sigma[1]', 'sigma[2]', 'w1[1]'), 
                     facet_args = list(ncol = 3)) + facet_text(size = 15)+
  theme_bw() +                       # white background
  theme(panel.grid = element_blank()) # optional: remove grid lines



ggsave(
  #filename = "trace_HS.png",  # or "traceplot.pdf"
  filename= paste0('trace_S3_HS_',
                   iter, '.png'),
  plot = p_HSGP,
  path = "/project/6003552/mingchi/LMC-sim",  # change to your desired folder
  width = 12,   # in inches
  height = 8,
  dpi = 300
)




# IG: 

library(bayesplot)
draws_df_HS <- cmdstan_fit_HS$draws(format = "df")
color_scheme_set("blue")
# summary of the posterior 
samples_beta0_HS<- draws_df_HS[['beta0']]
samples_beta1_HS<- draws_df_HS[['beta1']]
samples_corr_HS<- draws_df_HS[['Astar[2,1]']]
samples_rho1_HS<- draws_df_HS[['ell[1]']]
samples_rho2_HS<- draws_df_HS[['ell[2]']]
samples_sigma1_HS<- draws_df_HS[['sigma[1]']]
samples_sigma2_HS<- draws_df_HS[['sigma[2]']]





pdf( paste0('/project/6003552/mingchi/Simulations-Three/prior_vs_posterior',iter,'.pdf'), 
     width = 10, height = 8)
par(mfrow=c(4, 2))
# beta0
#samples_beta0<- draws_df_HS[['beta0']]
plot(density(samples_beta0_HS), col='red', lwd=2, xlab='beta0', ylab='Density',main='Prior and posterior for beta0' )
lines(density(rnorm(10000, 0, sd=3)), col='blue', lwd=2)
abline(v=beta_pts, col='black', lwd=5)

# beta1
#samples_beta1_HS<- draws_df_HS[['beta1']]
plot(density(samples_beta1_HS), col='red', lwd=2,xlab='beta1', ylab='Density',main='Prior and posterior for beta1' )
lines(density(rnorm(10000, 0, 3)), col='blue', lwd=2)
abline(v=beta_marks, col='black', lwd=5)

# R matrix, a21: 
#samples_corr_HS<- draws_df_HS[['Astar[2,1]']]
#hist(samples_corr, breaks=30, xlab='corr', ylab='Density', freq = F, main='Prior and posterior for a21')
plot(density(samples_corr_HS), col='red', lwd=2,xlab='corr', ylab='Density', main='Prior and posterior for a21' )
lines(density(rlkjcorr(100000, 2,eta=1.2)[,1,2]), col='blue', lwd=2)
abline(v=Astar[2,1], col='black', lwd=5)

# sigma1: 
#samples_sigma1_HS<- draws_df_HS[['sigma[1]']]
plot(density(samples_sigma1_HS), col='red', lwd=2, xlab='sigma1', ylab='Density', main='Prior and posterior for sigma1')
#lines(density(ztnormdist(10000, 4, sd=0.8, a=0)), col='blue', lwd=2)
#lines(density(rlnorm(n, meanlog = log(4), sdlog = 0.2)), col='blue', lwd=2)
lines(density(rhnorm(100000, sigma = 1)),col='blue', lwd=2 )
abline(v=sigma[1], col='black', lwd=5)

# sigma2:
#samples_sigma2_HS<- draws_df_HS[['sigma[2]']]
plot(density(samples_sigma2_HS), col='red', lwd=2, xlab='sigma2', ylab='Density', main='Prior and posterior for sigma2')
#lines(density(ztnormdist(10000, 2, sd=0.4, a=0)), col='blue', lwd=2)
#lines(density(rlnorm(n, meanlog = log(2), sdlog = 0.2)), col='blue', lwd=2)
lines(density(rhnorm(100000, sigma = 1)),col='blue', lwd=2 )
abline(v=sigma[2], col='black', lwd=5)



# rho1:
#samples_rho1_HS<- draws_df_HS[['ell[1]']]
plot(density(samples_rho1_HS), col='red', lwd=2, xlab='rho1', ylab='Density', main='Prior and posterior for rho1')
#lines(density(rlnorm(n, meanlog = log(1.8), sdlog = 0.4)), col='blue', lwd=2)
lines(density(rinvgamma(10000, alpha = 2, beta=2)),col='blue', lwd=2 )
abline(v=lscale[1], col='black', lwd=5)


# rho2:
#samples_rho2_HS<- draws_df_HS[['ell[2]']]
plot(density(samples_rho2_HS), col='red', lwd=2, xlab='rho1', ylab='Density', main='Prior and posterior for rho2')
#lines(density(rlnorm(n, meanlog = log(2), sdlog = 0.4)), col='blue', lwd=2)
lines(density(rinvgamma(10000, alpha = 2, beta=1)),col='blue', lwd=2 )
abline(v=lscale[2], col='black', lwd=5)

par(mfrow=c(1, 1))
dev.off()




# calculate the correlation structure: 
#grid_size

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




save(corr_dist_1, 
     file=paste0('/project/6003552/mingchi/Simulations-Three/corr_dist_',iter, '.Rdata'))






#mean(corr_dist)
#par(mfrow=c(1,1))
#boxplot(corr_dist, main='Correlation between points and scaled gene expression')
#abline(h=true.corr, col='red', lwd=2)
