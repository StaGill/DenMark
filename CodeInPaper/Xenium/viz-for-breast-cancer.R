######################
# visualization of the Breast cancer dataset 
# the cell type and the ROI
######################


# install packages
library('ggplot2')
library(cowplot)
library(fields)
library(viridis)
library('rstan')
library('MASS')
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
library(loo)
library(rethinking)
library(ggplot2)
#library(spatstat)
# load the cell location and gene expression per location: 

cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/cell_spatial.csv')
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000


# 542 genes in the dataset 
colnames(cell_gene)

# extract the  proliferative invasive tumor cell type: 
cell_spatial_tumor<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]

plot(cell_spatial_tumor$x, 
     cell_spatial_tumor$y)


gene_spatial_tumor<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]



# viz the point process
p0<- ggplot(cell_spatial_tumor, aes(x=x, y=y))+
  geom_point(size=0.9)+
  labs(title='Invasive Tumor',x='',y='')+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )+
  theme_classic(base_size = 14)


df<- as.data.frame(cbind(cell_spatial_tumor$x, 
           cell_spatial_tumor$y,
           gene_spatial_tumor))

colnames(df)[1]<- 'x'
colnames(df)[2]<- 'y'


p1<- ggplot(df, aes(x = x, y = y, color = CCND1)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(low = "snow2", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none")+
  theme_classic(base_size = 14)


p2<- ggplot(df, aes(x = x, y = y, color = CDH1)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(low = "snow2", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none")+
  theme_classic(base_size = 14)


p3<- ggplot(df, aes(x = x, y = y, color = MMP2)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(low = "snow2", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none")+
  theme_classic(base_size = 14)


p4<- ggplot(df, aes(x = x, y = y, color = EGFR)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(low = "snow2", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none")+
  theme_classic(base_size = 14)



p5<- ggplot(df, aes(x = x, y = y, color = POSTN)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(low = "snow2", high = "red", trans = "sqrt") +
  labs(title="", x="", y="")+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none")+
  theme_classic(base_size = 14)

plot_grid(p0, p1, p2, p3, p4, p5, nrow=2, align='v')

# Run the model on CCND1 gene:-----

#  scale the center slice 
df<- as.data.frame(cbind( cell_spatial_tumor$x, 
            cell_spatial_tumor$y ,
            gene_spatial_tumor$CCND1 ))
colnames(df)<- c('x','y','Gene')


center_slice_cell_location<- as.data.frame(cbind(cell_spatial_tumor$x, 
                                   cell_spatial_tumor$y ))

colnames(center_slice_cell_location)<- c('x','y')



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


# viz for the grids:-----

# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

# Custom layout:
#  Plot1 | Plot2
# Make widths slightly larger than legends so they fill space equally
layout(matrix(1:2, nrow = 1), widths = c(1, 1), respect = F)

# Points
par(mar = c(4, 4, 2, 1))  # small right margin
image.plot(
  loc_mat,
  breaks = make_breaks(loc_mat, 30),
  col = viridis(30),
  main = "Gridded Cell Counts",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 2,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)

# Marks
par(mar = c(4, 4, 2, 1))  # small right margin to keep close to left plot
image.plot(
  marks_mat,
  breaks = make_breaks(marks_mat, 30),
  col = viridis(30),
  main = "Gridded Gene Expression",
  xlab = "",
  ylab = '',
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 2,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(marks_mat), ny = nrow(marks_mat), col = "black", lwd = 0.5)

# Reset
par(mfrow = c(1, 1))






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

#cmdstan_fit_HS$save_object(paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/analysis_CCDN1.rds'))




summary_gene_22<-cmdstan_fit_HS$draws(format = "df")


draws_df_HS<- cmdstan_fit_HS$draws(format = "df")


# Check the posterior for the latent grids: 
w1_summary <- cmdstan_fit_HS$summary(variables = paste0("w1[",1:grid_size^2,"]"), 
                                     c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))


w2_summary <- cmdstan_fit_HS$summary(variables = paste0("w2[",1:grid_size^2,"]"), 
                                     c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))

sigma1_summary<- cmdstan_fit_HS$summary(variables = paste0("sigma[1]"), 
                                        c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))

sigma2_summary<- cmdstan_fit_HS$summary(variables = paste0("sigma[2]"), 
                                        c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))


a21_star_summary<- cmdstan_fit_HS$summary(variables = paste0("Astar[2,1]"), 
                                          c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))


# "A" matrix: a11=sigma1 ; a21= A21star*sigma2; a22= sqrt(1-A21star^2)*sigma2



beta0_summary <- cmdstan_fit_HS$summary(variables = 'beta0', 
                                        c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))


beta1_summary <- cmdstan_fit_HS$summary(variables = 'beta1', 
                                        c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))



# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

#par(mfrow = c(1, 2))
#png("tight_two_plots.png", width = 2000, height = 1000, res = 150)

# Custom layout:
#  Plot1 | Plot2
# Make widths slightly larger than legends so they fill space equally
layout(matrix(1:2, nrow = 1), widths = c(1, 1), respect = F)

# Points
par(mfrow=c(1,1), mar = c(4, 4, 4, 2))   # small right margin
image.plot(
  loc_mat,
  breaks = make_breaks(loc_mat, 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 0.4,
  legend.mar = 3,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  cex.main = 3,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)

# Marks
par(mfrow=c(1,1), mar = c(4, 4, 4, 2))   # small right margin to keep close to left plot
image.plot(
  marks_mat,
  breaks = make_breaks(marks_mat, 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = '',
  useRaster = FALSE,
  legend.width = 0.4,
  legend.mar = 3,
  cex.main = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(marks_mat), ny = nrow(marks_mat), col = "black", lwd = 0.5)

# Reset
par(mfrow = c(1, 1))





# That is not true, change into this plots: ----



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







#################################################
# Short-range analysis: triple-positive tumors------
##################################################

# criteria of the triple-positive:
triple_pos_idx<- (cell_gene$ERBB2>=1) & (cell_gene$ESR1 >=1) & (cell_gene$PGR >=1)

# extract the triple-positive cells: 
cell_spatial_tripos<-cell_spatial[triple_pos_idx,]

gene_spatial_tripos <-cell_gene[triple_pos_idx,]


plot(cell_spatial_tripos$x, 
     cell_spatial_tripos$y)


cell_spatial_tripos$ERBB2<- gene_spatial_tripos$ERBB2
cell_spatial_tripos$ESR1<- gene_spatial_tripos$ESR1
cell_spatial_tripos$PGR<- gene_spatial_tripos$PGR


ggplot(cell_spatial_tripos, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 2) +
  labs(
    title = "Triple-positive breast cancer cells (Zoomed Region)",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(4, 4.5), ylim = c(4.7, 5.2)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )




ggplot(cell_spatial_tripos, aes(x = x, y = y, color = ERBB2)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "magma") +  # continuous color scale
  labs(
    title = "Triple-positive breast cancer cells (ERBB2 expression)",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "ERBB2"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


ggplot(cell_spatial_tripos, aes(x = x, y = y, color = ESR1)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "magma") +  # continuous color scale
  labs(
    title = "Triple-positive breast cancer cells (ESR1 expression)",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "ESR1"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


ggplot(cell_spatial_tripos, aes(x = x, y = y, color = PGR)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "magma") +  # continuous color scale
  labs(
    title = "Triple-positive breast cancer cells (PGR expression)",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "PGR"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )





#################################################
# Short-range analysis: DCIS cells, colored genes ------
##################################################







# III. Viz for the three tumor regions: -------
# In the paper: high resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis

length(cell_spatial$X)


# extract the  proliferative invasive tumor cell type: 
cell_spatial_DCIS<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
cell_spatial_invasive<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
cell_spatial_immune<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),] 

gene_spatial_DCIS<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2' ),]
gene_spatial_Invasive<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_spatial_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]



# viz the point process
p_dcis<- ggplot(cell_spatial_DCIS, aes(x=x, y=y))+
  geom_point(size=0.3)+
  labs(title='',x=' x coordinates (mm)',y='y coordinates (mm)')+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    #plot.title = element_text(size = 45),
    legend.position = "none"
  )


p_invasive<- ggplot(cell_spatial_invasive, aes(x=x, y=y))+
  geom_point(size=0.3)+
  labs(title=' ',x=' x coordinates (mm)',y='y coordinates (mm)')+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    #plot.title = element_text(size = 45),
    legend.position = "none"
  )


p_immune<- ggplot(cell_spatial_immune, aes(x=x, y=y))+
  geom_point(size=0.3)+
  labs(title=' ',x=' x coordinates (mm)',y='y coordinates (mm)')+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    #plot.title = element_text(size = 45),
    legend.position = "none"
  )



plot_grid(p_dcis,
          p_invasive,
          p_immune, nrow=1)


# the gene expression: 
cell_spatial_DCIS$MZB1<- gene_spatial_DCIS$MZB1
cell_spatial_immune$CD4<- gene_spatial_immune$CD4
cell_spatial_invasive$GATA3<- gene_spatial_Invasive$GATA3




ggplot(cell_spatial_DCIS, aes(x = x, y = y, color = MZB1)) +
  geom_point(size = 0.3) +
  scale_color_gradient(low = "snow2", high = "#8B3A6E", trans = "sqrt",
                       guide = guide_colorbar(barheight = 10, barwidth = 0.6)) +
  theme_classic() +
  labs(title = "", x = "", y = "", color = "MZB1") +
  theme(
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 45),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )



ggplot(cell_spatial_immune, aes(x = x, y = y, color = CD4)) +
  geom_point(size = 0.3) +
  scale_color_gradient(low = "snow2", high = "#6B8010", trans = "sqrt",
                       guide = guide_colorbar(barheight = 10, barwidth = 0.6)) +
  theme_classic() +
  labs(title = "", x = "", y = "", color = "CD4") +
  theme(
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 45),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )


ggplot(cell_spatial_invasive, aes(x = x, y = y, color = GATA3)) +
  geom_point(size = 0.3) +
  scale_color_gradient(low = "snow2", high = "#1C8AAD", trans = "sqrt",
                       guide = guide_colorbar(barheight = 10, barwidth = 0.6)) +
  theme_classic() +
  labs(title = "", x = "", y = "", color = "GATA3") +
  theme(
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 45),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )


# pick the sub-regions for three cell types: DCIS 1, DCIS 2 and invasive BC-----



# [6.5,7.5]*[1.8,2.8], the 1*1 window
p_zoom_dcis1<- ggplot(cell_spatial, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 1.5) +
  labs(
    title = "DCIS-1 ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(6.5, 7.5), ylim = c(1.8, 2.8)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )

# [5,6]*[0.75,1.75] window
p_zoom_dcis2<- ggplot(cell_spatial, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 1.5) +
  labs(
    title = "DCIS-2 ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(5, 6), ylim = c(0.75, 1.75)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )


# [0.5,1.5]*[2,3] window
p_zoom_invasive<- ggplot(cell_spatial, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 2) +
  labs(
    title = "Invasive tumor ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0.5, 1.5), ylim = c(2, 3)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )

# the gridded for the three regions ----
# Define grid size and coordinates 
cell_spatial_ROI3<- cell_spatial[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5& cell_spatial$y>=1.8 & cell_spatial$y<= 2.8),]
cell_spatial_ROI2<- cell_spatial[which(cell_spatial$x>=5 & cell_spatial$x<= 6& cell_spatial$y>=0.75 & cell_spatial$y<= 1.75),]
cell_spatial_ROI1<- cell_spatial[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5& cell_spatial$y>=2 & cell_spatial$y<= 3),]


center_slice_cell_location <- cell_spatial_ROI1

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
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


# viz for the grids:-----

# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

par(
  mar = c(1, 1, 1.2, 3),   # bottom, left, top, right
  mgp = c(1.5, 0.3, 0)     # axis title, labels, line
)

# Points: DCIS 1
image.plot(
  log(loc_mat + 0.05),
  breaks = make_breaks(log(loc_mat + 0.05), 30),
  col = viridis(30, option = "magma"),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 1.5,
  #legend.shrink = 0.8,
  axis.args = list(cex.axis = 2.5),  # ← increase number size (default is 1)
  xaxt = "n",
  yaxt = "n"
)

grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)


#grids.roi1<-as.vector(log(loc_mat + 0.05))

grids.roi1<- as.vector(loc_mat )


## Gridded Region 2 ----

center_slice_cell_location <- cell_spatial_ROI2

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
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


Y<- loc_100_v
M<- marks_100_v


# viz for the grids:-----

# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

par(
  mar = c(1, 1, 1.2, 3),   # bottom, left, top, right
  mgp = c(1.5, 0.3, 0)     # axis title, labels, line
)

# Points: DCIS 1
image.plot(
  log(loc_mat + 0.05),
  breaks = make_breaks(log(loc_mat + 0.05), 30),
  col = viridis(30, option = "magma"),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 1.5,
  #legend.shrink = 0.8,
  axis.args = list(cex.axis = 2.5),  # ← increase number size (default is 1)
  xaxt = "n",
  yaxt = "n"
)

grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)


#grids.roi2<-as.vector(log(loc_mat + 0.05))
grids.roi2<- as.vector(loc_mat )



## Gridded Region 3 ----


center_slice_cell_location <- cell_spatial_ROI3

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
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}



# check if points and marks are log-linearlly correlated 
loc_100_v<- as.vector(loc_100)
marks_100_v<- as.vector(marks_100)


Y<- loc_100_v
M<- marks_100_v


# viz for the grids:-----

# Transform matrices
loc_mat   <- apply(matrix(loc_100, nrow = grid_size, byrow = TRUE), 2, rev)
marks_mat <- apply(matrix(marks_100, nrow = grid_size, byrow = TRUE), 2, rev)

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

par(
  mar = c(1, 1, 1.2, 3),   # bottom, left, top, right
  mgp = c(1.5, 0.3, 0)     # axis title, labels, line
)

# Points: DCIS 1
image.plot(
  log(loc_mat + 0.05),
  breaks = make_breaks(log(loc_mat + 0.05), 30),
  col = viridis(30, option = "magma"),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 1.5,
  #legend.shrink = 0.8,
  axis.args = list(cex.axis = 2.5),  # ← increase number size (default is 1)
  xaxt = "n",
  yaxt = "n"
)

grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)



#grids.roi3<-as.vector(log(loc_mat + 0.05))

grids.roi3<- as.vector(loc_mat )




v1<- grids.roi1
v2<- grids.roi2
v3<- grids.roi3


library(ggplot2)
library(ggridges)

df <- data.frame(
  count = c(v1, v2, v3),
  ROI = factor(
    rep(c("ROI1", "ROI2", "ROI3"),
        times = c(length(v1), length(v2), length(v3))),
    levels = c("ROI1", "ROI2", "ROI3")
  )
)
ggplot(df, aes(x = count, y = ROI, fill = ROI)) +
  geom_density_ridges(
    scale = 1.2,
    alpha = 0.8,
    color = "white",
    size = 0.3
  ) +
  theme_classic() +
  labs(
    x = "Cell counts per grid",
    y = "Region of interest"
  ) +
  theme(
    legend.position = "none"
  )


# Combine data into a list
library(ggplot2)

# Combine into a single data frame
df <- data.frame(
  count = c(v1, v2, v3),
  ROI = factor(
    rep(c("ROI 1", "ROI 2", "ROI 3"),
        times = c(length(v1), length(v2), length(v3))),
    levels = c("ROI 1", "ROI 2", "ROI 3")
  )
)



library(grid)  # for unit()

ggplot(df, aes(x = count, y = ROI, fill = ROI)) +
  geom_boxplot(                              # show outliers by default
    outlier.shape = 16,                      # default solid circle
    outlier.size = 2,                        # size of outlier points
    outlier.colour = "black"                 # color of outliers
  ) +
  #geom_jitter(height = 0.2, alpha = 0.3, size = 1) + # optional: add points
  scale_fill_manual(values = c( "#006400", "#DAA520","#B8860B")) +
  labs(x = "Cell counts", y = "",
       title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.key.size = unit(2, "lines"),
    legend.position = "none"
  )







p_whole_tissue<- ggplot(cell_spatial, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 0.6) +
  scale_color_discrete(labels = function(x) gsub("_", " ", x)) +
  # First box
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = 1.8, ymax = 2.8, 
           color = "black", fill = NA, linewidth = 2) +
  # Second box
  annotate("rect", xmin = 5, xmax = 6, ymin = 0.75, ymax = 1.75, 
           color = "black", fill = NA, linewidth = 2) +
  # Third box
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 2, ymax = 3, 
           color = "black", fill = NA, linewidth = 2) +
  labs(
    title = "",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  theme_classic()  +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.key.size = unit(2, "lines"),
    legend.position = "left"
  )+
  guides(color = guide_legend(override.aes = list(size = 5)))

# plot: 
plot_grid(p_zoom_dcis1, 
          p_zoom_dcis2, 
          p_zoom_invasive, nrow=1)

# DCIS I, DCIS II, and invasive 
length(cell_spatial[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5 & cell_spatial$y>= 1.8 & cell_spatial$y <= 2.8),][,1])
length(cell_spatial[which(cell_spatial$x>=5 & cell_spatial$x<= 6 & cell_spatial$y>= 0.75 & cell_spatial$y <= 1.75),][,1])
length(cell_spatial[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y <= 3),][,1])


# color for the DCIS I, DCIS II, invasive tumor and CD8+ T cells only 
# Define your target colors
target_colors <- c(
  "DCIS_1" = "#B8860B",
  "DCIS_2" = "#DAA520",
  "Invasive_Tumor" = "#006400",
  "CD8+_T_Cells" = "#8B3A6E"
)

# Create a new variable: show only 4 cell types + "Other"
cell_spatial$celltype_legend <- ifelse(
  cell_spatial$celltype %in% names(target_colors),
  cell_spatial$celltype,
  "Others"
)

# Colors for legend (5 categories)
legend_colors <- c(
  target_colors,
  "Others" = "grey70"
)
# Split data into grey vs colored groups
df_grey  <- subset(cell_spatial, celltype_legend == "Others")
df_color <- subset(cell_spatial, celltype_legend != "Others")

ggplot() +
  # 1. Grey points FIRST
  geom_point(data = df_grey, aes(x = x, y = y, color = celltype_legend),
             size = 0.6) +
  # 2. Colored points SECOND (on top)
  geom_point(data = df_color, aes(x = x, y = y, color = celltype_legend),
             size = 0.6) +
  
  # Manual colors
  scale_color_manual(
    values = legend_colors,
    labels = function(x) gsub("_", " ", x)
  ) +
  
  # Boxes
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = 1.8, ymax = 2.8,
           color = "black", fill = NA, linewidth = 2) +
  annotate("rect", xmin = 5, xmax = 6, ymin = 0.75, ymax = 1.75,
           color = "black", fill = NA, linewidth = 2) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 2, ymax = 3,
           color = "black", fill = NA, linewidth = 2) +
  
  labs(
    title = "",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 30),
    #legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))


# [6.5,7.5]*[1.8,2.8], the 1*1 window
ggplot() +
  geom_point(data = df_grey, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # 2. Colored points SECOND (on top)
  geom_point(data = df_color, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # Manual colors
  scale_color_manual(
    values = legend_colors,
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "DCIS-1 ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(6.5, 7.5), ylim = c(1.8, 2.8)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )

# [5,6]*[0.75,1.75] window

# DCIS-II
ggplot() +
  geom_point(data = df_grey, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # 2. Colored points SECOND (on top)
  geom_point(data = df_color, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # Manual colors
  scale_color_manual(
    values = legend_colors,
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "DCIS-2 ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(5, 6), ylim = c(0.75, 1.75)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )


# [0.5,1.5]*[2,3] window
ggplot() +
  geom_point(data = df_grey, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # 2. Colored points SECOND (on top)
  geom_point(data = df_color, aes(x = x, y = y, color = celltype_legend),
             size = 1.5) +
  # Manual colors
  scale_color_manual(
    values = legend_colors,
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Invasive tumor ROI",
    x = "X coordinate (mm)",
    y = "Y coordinate (mm)",
    color = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0.5, 1.5), ylim = c(2, 3)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )





# Create data
library(ggplot2)

# Create data
cell_counts <- data.frame(
  category = c("DCIS 1", "DCIS 2", "Invasive tumor", "Other"),
  count = c(5031, 3812, 6645, 167780 - 5031 - 3812 - 6645)
)

# Calculate percentages and positions
cell_counts$percentage <- round(cell_counts$count / sum(cell_counts$count) * 100, 1)
cell_counts$fraction <- cell_counts$count / sum(cell_counts$count)
cell_counts$ymax <- cumsum(cell_counts$fraction)
cell_counts$ymin <- c(0, head(cell_counts$ymax, n = -1))
cell_counts$labelPosition <- (cell_counts$ymax + cell_counts$ymin) / 2

# Create labels
cell_counts$label <- paste0(cell_counts$percentage, "%")

# Define colors
colors <- c("DCIS 1" = "#B8860B",        # Dark goldenrod
            "DCIS 2" = "#DAA520",        # Goldenrod
            "Invasive tumor" = "#006400",  # Dark green
            "Other" = "#D3D3D3")         # Light gray

# Create donut chart
ggplot(cell_counts, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = category)) +
  geom_rect(color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  xlim(c(0, 8)) +
  scale_fill_manual(values = colors) +
  labs(title = "",
       fill = "ROI") +
  theme_void() +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.box.margin = margin(l = -80, t = 0, b = 0, r = 0),
    legend.margin = margin(t = -10)  
  )

# test if invasive region is homogeneity 
#invasive_idx_region<-cell_spatial[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y<= 3),]


#window<- owin(xrange=c(0.5, 1.5),
#              yrange=c(2, 3))

#pp<- ppp(invasive_idx_region$x, 
#         invasive_idx_region$y, 
#         window=window)
#quadrat.test(pp) # test for CSR, p-value=1e-6<0.05


# test the immune cell distribution at those regions:




#######################################################
# long-range analysis viz (for specific region) --------
#######################################################


# extract the  proliferative invasive tumor cell type: 
cell_spatial_DCIS<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
cell_spatial_invasive<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
cell_spatial_immune<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),] 

gene_spatial_DCIS<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2' ),]
gene_spatial_Invasive<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_spatial_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]








