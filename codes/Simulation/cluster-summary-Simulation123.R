###########################################################################
# Summarize the results from Simulation 1, Simulation 2, and Simulation 3
###########################################################################

# Step 0: some package used in this file
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(fields)
library(reshape2)
library(patchwork)
# quantile functions:
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)



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


# I. Summarize the results from Simulation 1-----
# simulation is to compare exact GP and HSGP, on general model (M4) 
# Based on 15*15 grids 
# compare: 
# (1) computational time 
# (2) fitting parameters 


#rm(list=ls())




# Combine the boxplots, final codes for simulation 1: -----
beta_pts   <- 1      # intercept for points (eta1)
beta_marks <- -1     # intercept for marks (eta2)
offset     <- 0.05   # same small offset used in Stan
sigma <- c(1, 1)
R <- cbind(c(1, 0.5),
           c(0.5, 1))
lscale <- c(0.5, 0.3)
Astar <- t(chol(R))
A <- diag(sigma) %*% Astar


# exact GP
post.summary.exact.GP<- data.frame()
post.summary.exact.GP<- NULL
upper.summary.exact.GP<- data.frame()
upper.summary.exact.GP<- NULL
lower.summary.exact.GP<- data.frame()
lower.summary.exact.GP<- NULL
time.consume.exact.GP<- data.frame()
time.consume.exact.GP<- NULL

# high grids
post.summary.HSGP.h<- data.frame()
post.summary.HSGP.h<- NULL
upper.summary.HSGP.h<- data.frame()
upper.summary.HSGP.h<- NULL
lower.summary.HSGP.h<- data.frame()
lower.summary.HSGP.h<- NULL
time.consume.HSGP.h<- data.frame()
time.consume.HSGP.h<- NULL

# middle grids
post.summary.HSGP.m<- data.frame()
post.summary.HSGP.m<- NULL
upper.summary.HSGP.m<- data.frame()
upper.summary.HSGP.m<- NULL
lower.summary.HSGP.m<- data.frame()
lower.summary.HSGP.m<- NULL
time.consume.HSGP.m<- data.frame()
time.consume.HSGP.m<- NULL

# middle grids
post.summary.HSGP.m2<- data.frame()
post.summary.HSGP.m2<- NULL
upper.summary.HSGP.m2<- data.frame()
upper.summary.HSGP.m2<- NULL
lower.summary.HSGP.m2<- data.frame()
lower.summary.HSGP.m2<- NULL
time.consume.HSGP.m2<- data.frame()
time.consume.HSGP.m2<- NULL

# small grids
post.summary.HSGP.l<- data.frame()
post.summary.HSGP.l<- NULL
upper.summary.HSGP.l<- data.frame()
upper.summary.HSGP.l<- NULL
lower.summary.HSGP.l<- data.frame()
lower.summary.HSGP.l<- NULL
time.consume.HSGP.l<- data.frame()
time.consume.HSGP.l<- NULL


# high # of basis functions 
sim_list<- c(3,4,6:12,14,16:96,99:100)

for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/results/exact_GP_',i,'.csv')
  post.summary.exact.GP<- cbind(post.summary.exact.GP, read.csv(file=path, header=T)[,3])
  lower.summary.exact.GP <- cbind(lower.summary.exact.GP, read.csv(file=path, header=T)[,5])
  upper.summary.exact.GP <- cbind(upper.summary.exact.GP, read.csv(file=path, header=T)[,9])
  time.consume.exact.GP <- cbind(time.consume.exact.GP, read.csv(file=path, header=T)[,13])
}


for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/revise-V3/HSGP_approx_30_',i,'.csv')
  post.summary.HSGP.h <- cbind(post.summary.HSGP.h, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.h <- cbind(lower.summary.HSGP.h, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.h <- cbind(upper.summary.HSGP.h, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.h <- cbind(time.consume.HSGP.h, read.csv(file=path, header=T)[,13])
}

# middle # of basis functions
for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/revise-V3/HSGP_approx_25_',i,'.csv')
  post.summary.HSGP.m <- cbind(post.summary.HSGP.m, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.m <- cbind(lower.summary.HSGP.m, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.m <- cbind(upper.summary.HSGP.m, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.m <- cbind(time.consume.HSGP.m, read.csv(file=path, header=T)[,13])
}


# low # of basis functions
for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/revise-V3/HSGP_approx_5_',i,'.csv')
  post.summary.HSGP.m2 <- cbind(post.summary.HSGP.m2, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.m2 <- cbind(lower.summary.HSGP.m2, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.m2 <- cbind(upper.summary.HSGP.m2, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.m2 <- cbind(time.consume.HSGP.m2, read.csv(file=path, header=T)[,13])
}


# boundary factors: 
# high boundary factors 
for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/revise-V3/HSGP_approx_3_',i,'.csv')
  post.summary.HSGP.l <- cbind(post.summary.HSGP.l, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.l <- cbind(lower.summary.HSGP.l, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.l <- cbind(upper.summary.HSGP.l, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.l <- cbind(time.consume.HSGP.l, read.csv(file=path, header=T)[,13])
}




# ---- Combine posterior means for a21 (row 5) ----
a21_basis_values <- list(
  gl    = as.numeric(post.summary.HSGP.l[3, ]),
  gm2    = as.numeric(post.summary.HSGP.m2[3, ]),
  gm    = as.numeric(post.summary.HSGP.m[3, ]),
  gh    = as.numeric(post.summary.HSGP.h[3, ]),
  gexact = as.numeric(post.summary.exact.GP[10, ])
)

# ---- Combine time consumption (take row 1, since all rows are identical) ----
time_values_b <- list(
  gtl    = as.numeric(time.consume.HSGP.l[1, ])/60,
  gtm2    = as.numeric(time.consume.HSGP.m2[1, ])/60,
  gtm    = as.numeric(time.consume.HSGP.m[1, ])/60,
  gth    = as.numeric(time.consume.HSGP.h[1, ])/60,
  gtexact = as.numeric(time.consume.exact.GP[1, ])/60
)

# ---- Labels ----
#method_labels <- c( "4 basis", "25 basis", "100 basis", '625 basis', 'exact GP')


# convert into ggplot2:
# ---- Method labels in your desired order ----
method_labels <- c("4", "25", "100", "625", "exact GP")

# ---- Personalized colors ----
my_colors <- c(
  "4"   = "#C8EAF5",  # green
  "25"  = "#C8EAF5",  # orange
  "100" = "#C8EAF5",  # blue
  "625" = "#C8EAF5",  # pink
  "exact GP"  = "#0096C9"   # gray
)

my_colors_t <- c(
  "4"   = "#E2A7CC",  # green
  "25"  = "#E2A7CC",  # orange
  "100" = "#E2A7CC",  # blue
  "625" = "#E2A7CC",  # pink
  "exact GP"  = "#9B5678"   # gray
)



# ---- Combine posterior means into data frame ----
a21_df <- data.frame(
  method = factor(rep(method_labels, each = length(a21_basis_values$gl)),
                  levels = method_labels),
  a21 = c(a21_basis_values$gl,
          a21_basis_values$gm2,
          a21_basis_values$gm,
          a21_basis_values$gh,
          a21_basis_values$gexact)
)

# ---- Combine computation times into data frame ----
time_df <- data.frame(
  method = factor(rep(method_labels, each = length(time_values_b$gtl)),
                  levels = method_labels),
  time = c(time_values_b$gtl,
           time_values_b$gtm2,
           time_values_b$gtm,
           time_values_b$gth,
           time_values_b$gtexact)
)

# ---- Plot posterior means ----
S1_basis <- ggplot(a21_df, aes(x = method, y = a21, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5,size=1.2) +
  scale_fill_manual(values = my_colors) +
  labs(x = 'Basis', y = "",
            title = '')+
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )

# ---- Plot computation times ----
S1_basis_t <- ggplot(time_df, aes(x = method, y = time, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size=1.2) +
  scale_y_log10() +  # log scale for skewed times
  scale_fill_manual(values = my_colors_t) +
  labs(x = 'Basis', y = "",
       title = '')+
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )

# ---- Combine side by side ----
plot_grid(S1_basis, S1_basis_t, ncol = 2)




# for the boundary factors: 

# high c
post.summary.HSGP.hc<- data.frame()
post.summary.HSGP.hc<- NULL
upper.summary.HSGP.hc<- data.frame()
upper.summary.HSGP.hc<- NULL
lower.summary.HSGP.hc<- data.frame()
lower.summary.HSGP.hc<- NULL
time.consume.HSGP.hc<- data.frame()
time.consume.HSGP.hc<- NULL

# middle c
post.summary.HSGP.mc<- data.frame()
post.summary.HSGP.mc<- NULL
upper.summary.HSGP.mc<- data.frame()
upper.summary.HSGP.mc<- NULL
lower.summary.HSGP.mc<- data.frame()
lower.summary.HSGP.mc<- NULL
time.consume.HSGP.mc<- data.frame()
time.consume.HSGP.mc<- NULL

# small c
post.summary.HSGP.lc<- data.frame()
post.summary.HSGP.lc<- NULL
upper.summary.HSGP.lc<- data.frame()
upper.summary.HSGP.lc<- NULL
lower.summary.HSGP.lc<- data.frame()
lower.summary.HSGP.lc<- NULL
time.consume.HSGP.lc<- data.frame()
time.consume.HSGP.lc<- NULL



# boundary factor
sim_list<- sim_list

for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/backup-backup-V3/HSGP_approx_30_',i,'.csv')
  post.summary.HSGP.lc <- cbind(post.summary.HSGP.lc, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.lc <- cbind(lower.summary.HSGP.lc, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.lc <- cbind(upper.summary.HSGP.lc, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.lc <- cbind(time.consume.HSGP.lc, read.csv(file=path, header=T)[,13])
}

# middle # of basis functions
for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/backup-backup-V3/HSGP_approx_25_',i,'.csv')
  post.summary.HSGP.mc <- cbind(post.summary.HSGP.mc, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.mc <- cbind(lower.summary.HSGP.mc, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.mc <- cbind(upper.summary.HSGP.mc, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.mc <- cbind(time.consume.HSGP.mc, read.csv(file=path, header=T)[,13])
}

# middle # of basis functions
for (i in sim_list){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation1/V3/backup-backup-V3/HSGP_approx_5_',i,'.csv')
  post.summary.HSGP.hc <- cbind(post.summary.HSGP.hc, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP.hc <- cbind(lower.summary.HSGP.hc, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP.hc <- cbind(upper.summary.HSGP.hc, read.csv(file=path, header=T)[,9])
  time.consume.HSGP.hc <- cbind(time.consume.HSGP.hc, read.csv(file=path, header=T)[,13])
}




# convert into ggplot2:
# ---- Method labels in your desired order ----
boundary_labels <- c("1.02", "1.05", "1.5", 'Exact GP')

time_values_b <- list(
  gtlc    = as.numeric(time.consume.HSGP.lc[1, ])/60,
  gtmc    = as.numeric(time.consume.HSGP.mc[1, ])/60,
  gthc    = as.numeric(time.consume.HSGP.hc[1, ])/60,
  gtexact = as.numeric(time.consume.exact.GP[1, ])/60
  
)


# ---- Personalized colors ----
my_colors <- c(
  "1.02"   = "#D5E6A8",  # green
  "1.05"  = "#D5E6A8",  # orange
  "1.5" = "#D5E6A8",  # blue
  'Exact GP'='#8BA04E'
)

my_colors_t <- c(
  "1.02"   = "#FFF193",  # green
  "1.05"  = "#FFF193",  # orange
  "1.5" = "#FFF193",  # blue
  'Exact GP'='#E8B92E'
)

# ---- Combine posterior means for a21 (row 5) ----
a21_b_values <- list(
  gtlc    = as.numeric(post.summary.HSGP.lc[3, ]),
  gtmc    = as.numeric(post.summary.HSGP.mc[3, ]),
  gthc    = as.numeric(post.summary.HSGP.hc[3, ]),
  gexact = as.numeric(post.summary.exact.GP[10, ])
)


# ---- Combine posterior means into data frame ----
a21_b_df <- data.frame(
  method = factor(rep(boundary_labels, each = length(a21_b_values$gtlc)),
                  levels = boundary_labels),
  a21 = c(a21_b_values$gtlc,
          a21_b_values$gtmc,
          a21_b_values$gthc,
          a21_b_values$gexact)
)


# ---- Combine computation times into data frame ----
time_b_df <- data.frame(
  method = factor(rep(boundary_labels, each = length(a21_b_values$gexact)),
                  levels = boundary_labels),
  time = c(time_values_b$gtlc,
           time_values_b$gtmc,
           time_values_b$gthc,
           time_values_b$gtexact)
)



# ---- Plot posterior means ----
S1_bound <- ggplot(a21_b_df, aes(x = method, y = a21, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size=1.2) +
  scale_fill_manual(values = my_colors) +
  #geom_hline(yintercept = R[2,1], color = "#000", linetype = "dashed", size = 1, alpha=.5) +
  #labs(x = 'Boundary factor', y = "",
  #     title = expression("Posterior mean of "~rho)) +
  labs(x = 'Boundary factor', y = "",
            title = '')+
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )


S1_bound_t<- ggplot(time_b_df, aes(x = method, y = time, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size=1.2) +
  #scale_y_log10() +  # log scale for skewed times
  scale_fill_manual(values = my_colors_t) +
  scale_y_log10() +
  #labs(x = 'Boundary factor', y = "",
  #     title = "Computation time (Minutes)") +
  labs(x = 'Boundary factor', y = "",
       title = "") +
  theme_classic(base_size = 14) +
  #ylim(1, 3)+
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )


# for this figure (Figure D, part II), we use the arrow for the fourth panel: 
y_arrow_start <- 2.5       # arrow tail
y_arrow_end   <- 5      # arrow head (top)
label_y       <- 4.5      # label above arrow

# the figure for the boxplot (time)
ggplot() +
  # ---- Boxplots for methods 1–3 ----
geom_boxplot(
  data = subset(time_b_df, method != 'Exact GP'),
  aes(x = method, y = time, fill = method),
  outlier.shape = NA, alpha = 0.5, size = 1.2
) +
  
  # ---- Arrow for method 4 instead of boxplot ----
geom_segment(
  aes(
    x = 'Exact GP', xend = 'Exact GP',
    y = y_arrow_start, yend = y_arrow_end
  ),
  arrow = arrow(length = unit(0.32, "inches")),
  linewidth = 2.5, 
  color='#E8B92E'
) +
  geom_text(
    aes(x = 'Exact GP', y = label_y, label = "100+"),
    nudge_x = 0.32, 
    size = 15, fontface = "bold"
  ) +
  
  scale_fill_manual(values = my_colors_t) +
  scale_y_log10() +
  labs(x = "Boundary factor", y = "", title = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 40),
    legend.position = "none"
  )


plot_grid(S1_bound, S1_bound_t)







# II. Summarize the results from Simulation 2-----

L1<- 4.4
L2<- 3.6

# the viz for the point-level and grid-level data:
S2_viz<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/Figure-S2-simulated-data.csv')

colnames(S2_viz)[4]<- 'Gene'
#my_cols <- c('#FFE5E6',"#fdcfd1", "#fc9ea3",   "#fa6e76", "#f93d48", '#ED1B2F')
#my_cols <- c( "#FFC1C4",   "#FF9AA2", "#FA6E76",'#F93D48',
#              '#ED1B2F','#C1272D', '#A10E1A','#800000' )

#my_cols <- c(
#  "#FF9AA2",  # Very Low
# "#FA4A50",  # Medium-Low
#  "#F93D48",  # Medium
#  "#ED1B2F",  # Medium-High
# # "#C1272D",  # High
#  "#8B0000"   # Peak
#)


my_cols <- c(
  "#EAC0BD",  # Very Low
  "#EC6756",  # Medium
  "#D42D24",  # Medium-High
  "#BC1119",  # High
  "#9C0824"   # Peak
)


my_cols_pts <- c(
  "#5F5F5F",  # Medium
  "#454545",  # Medium-High
  "#2C2C2C",  # High
  "#1B1B1B"   # Peak
)


# the simulated marked point process: 
S2_mpp<- ggplot(S2_viz, aes(x = x, y = y, color = Gene )) +
  geom_jitter(size = 3, width = 0.1, height = 0.1) +
  scale_color_gradientn(
    colors = my_cols,
    breaks = seq(0, 12, by = 4),
    labels = seq(0, 12, by = 4),
    guide = guide_colorbar(
      barwidth = 2,
      barheight = 20,
      title = "Gene\nExpression"
    )
  ) +
  labs(title = "", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 45)
  )


# the gridded marks:
grid_size_large <- 50 

grid_summary_large <- grid_pts_centroid(grid_size=grid_size_large, 
                                        dataset=S2_viz[,-1])

grid_data_large <- buildgridpp(grid_size = grid_size_large, 
                               dataset = S2_viz[,-1])

Y_high <-  grid_data_large[[2]]
M_high <-  grid_data_large[[3]]
n<-  grid_summary_large[[3]]
grid_area_high <- grid_summary_large[[2]]

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}

nr <- nrow(M_high)
nc <- ncol(M_high)

# X coordinates: map columns to -L1 ... L1
x_coords <- seq(-L1, L1, length.out = nc)

# Y coordinates: map rows to -L2 ... L2
y_coords <- seq(-L2, L2, length.out = nr)



# use ggplot2 for the gridded points and gridded gene expression: 
# Create a reusable ggplot2 function for spatial grid plots

# Updated function with actual coordinate units on axes
plot_spatial_grid_gg <- function(x_coords, y_coords, data_matrix, 
                                 title = "Spatial Plot",
                                 n_breaks = 10,
                                 color_option = "magma",
                                 legend_title = "",
                                 show_grid = TRUE,
                                 grid_color = "black",
                                 grid_size = 0.3,
                                 n_axis_breaks = 5) {
  
  # Convert matrix to long format for ggplot2
  df <- expand.grid(x = x_coords, y = y_coords)
  df$z <- as.vector(data_matrix)
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile(color = if(show_grid) grid_color else NA, 
              linewidth = if(show_grid) grid_size else 0) +
    scale_fill_viridis_c(option = color_option, 
                         name = legend_title,
                         breaks = pretty(df$z, n = n_breaks)) +
    # X-axis with actual coordinate values
    scale_x_continuous(
      breaks = pretty(x_coords, n = n_axis_breaks),
      expand = c(0, 0)
    ) +
    # Y-axis with actual coordinate values
    scale_y_continuous(
      breaks = pretty(y_coords, n = n_axis_breaks),
      expand = c(0, 0)
    ) +
    labs(title = title, 
         x = "X Coordinate", 
         y = "Y Coordinate") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.8, "cm"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 13, margin = margin(t = 10)),
      axis.title.y = element_text(size = 13, margin = margin(r = 10)),
      axis.ticks = element_line(),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}



# ===== Plot 1: Marked Point Process =====
S2_mpp <- ggplot(S2_viz, aes(x = x, y = y, color = Gene)) +
  geom_jitter(size = 1.5, width = 0.1, height = 0.1) +
  scale_color_gradientn(
    colors = my_cols,
    breaks = seq(0, 12, by = 4),
    labels = seq(0, 12, by = 4),
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 15,
      title = "Gene\nExpression"
    )
  ) +
  labs(title = "Simulated Cell Locations\n& Expression", 
       x = "X Coordinate", 
       y = "Y Coordinate") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, margin = margin(r = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  coord_fixed()

# ===== Plot 2: Gridded Cell Counts =====
p_gridded_cells <- plot_spatial_grid_gg(
  x_coords = x_coords,
  y_coords = y_coords,
  data_matrix = Y_high,
  title = "Gridded Cell Counts",
  n_breaks = 10,
  color_option = "magma",
  legend_title = "Count",
  n_axis_breaks = 5
)

# ===== Plot 3: Gridded Gene Expression =====
p_gridded_avg_express <- plot_spatial_grid_gg(
  x_coords = x_coords,
  y_coords = y_coords,
  data_matrix = M_high / (Y_high + 0.05),
  title = "Avg Gridded Expression",
  n_breaks = 10,
  color_option = "magma",
  legend_title = "Expression",
  n_axis_breaks = 5
)

# ===== Option 1: Using cowplot (Recommended) =====
p_combined <- plot_grid(
  S2_mpp,
  p_gridded_cells,
  p_gridded_avg_express, 
  align = 'hv',
  axis = 'tb',
  ncol = 3,
  rel_widths = c(1, 1, 1)  # Equal widths
)

print(p_combined)





# Simulation 2 revise the criteria ---
file_path <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation2/general-model/results-V2/'


# the summarized version of the dataset
datanum<- 100

loglambda1_high<- rep(NA, datanum)
loglambda1_mid<- rep(NA, datanum)
loglambda1_mid2<- rep(NA, datanum)
loglambda1_low<- rep(NA, datanum)


loglambda2_high<- rep(NA, datanum)
loglambda2_mid<- rep(NA, datanum)
loglambda2_mid2<- rep(NA, datanum)
loglambda2_low<- rep(NA, datanum)


# read csv 
for (i in 1:100){
  file_df<- read.csv(paste0(file_path, 'summary_grid_M2_compare_',i,'.csv'))
  
  # loglambda1
  loglambda1_high[i] <- file_df[,2][1]
  loglambda1_mid[i] <- file_df[,2][2]
  loglambda1_mid2[i] <- file_df[,2][3]
  loglambda1_low[i] <- file_df[,2][4]
  
  # loglambda2
  loglambda2_high[i] <- file_df[,3][1]
  loglambda2_mid[i] <- file_df[,3][2]
  loglambda2_mid2[i] <- file_df[,3][3]
  loglambda2_low[i] <- file_df[,3][4]
}



# Sample data setup
labels <- c("400", "900", "1600", "2500")
df <- data.frame(
  grid = factor(rep(labels, each = datanum), levels = labels),
  loglambda1 = c(loglambda1_low, loglambda1_mid2, loglambda1_mid, loglambda1_high),
  loglambda2 = c(loglambda2_low, loglambda2_mid2, loglambda2_mid, loglambda2_high)
)

df_long <- df %>%
  pivot_longer(cols = starts_with("loglambda"),
               names_to = "parameter",
               values_to = "bias")

my_colors <- c(
  "400"  = "#999999",
  "900"  = "#999999",
  "1600"  = "#999999",
  "2500" = "#999999"
)

# Plot lambda1 with proper Greek letters
S2_p1 <- ggplot(df_long %>% filter(parameter == "loglambda1"),
             aes(x = grid, y = bias, fill = grid)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, size=1.2) +
  scale_fill_manual(values = my_colors) +
  geom_hline(yintercept = 0, color = "#000", linetype = "dashed", size = 2, alpha=.5) +
  ylim(-1.5, 0)+
  labs(x = 'Grids', y = "",
       title = expression("Bias for "*log(lambda[1]))) +   # <- correct
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 45),
        axis.text.y = element_text(size = 45),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 50))



# Plot lambda2
S2_p2 <- ggplot(df_long %>% filter(parameter == "loglambda2"),
             aes(x = grid, y = bias, fill = grid)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, size=1.2) +
  geom_hline(yintercept = 0, color = "#000", linetype = "dashed", size = 2, alpha=.5) +
  scale_fill_manual(values = my_colors) +
  labs(x = 'Grids', y = '',
       title = expression("Bias for "*log(lambda[2]) )) +   # <- correct
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 45),
        axis.text.y = element_text(size = 45),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 50))

# Combine side by side
plot_grid(S2_p1, S2_p2, labels = c("", ""), ncol = 2)









# III. Summarize the results from Simulation 3-----
# simulate from grids and fit with grids 
# checking the convergence 

# from the dataset 
file_path <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation3/revised-grid-grid-simulation/results/'

# the parameter recovery: ----
#================ parameters; ================ 
sigma <- c(1, 1)
R <- cbind(c(1,0.5),
           c(0.5,1))
beta_pts<- 2
beta_marks<- 1
lscale <- c(2, 1) # 0.8, 0.4
offset<- 0.05
#==============================================

# the half-length of both side 
L1<- 4.4
L2<- 3.6

# save the posterior summary 
post.summary.HSGP<- data.frame()
post.summary.HSGP<- NULL

upper.summary.HSGP<- data.frame()
upper.summary.HSGP<- NULL

lower.summary.HSGP<- data.frame()
lower.summary.HSGP<- NULL

time.consume.HSGP<- data.frame()
time.consume.HSGP<- NULL


Rhat.HSGP<- data.frame()
Rhat.HSGP<- NULL


for (i in 1:100){
  path<- paste0('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation3/revised-grid-grid-simulation/results/S3-HSGP_',i,'.csv')
  post.summary.HSGP<- cbind(post.summary.HSGP, read.csv(file=path, header=T)[,3])
  lower.summary.HSGP <- cbind(lower.summary.HSGP, read.csv(file=path, header=T)[,5])
  upper.summary.HSGP <- cbind(upper.summary.HSGP, read.csv(file=path, header=T)[,9])
  Rhat.HSGP <- cbind(Rhat.HSGP, read.csv(file=path, header=T)[,10]  )
}


mean.beta0.HS  <- mean(post.summary.HSGP[1,])
mean.beta1.HS<- mean(post.summary.HSGP[2,])
mean.rho1.HS <- mean(post.summary.HSGP[3,])
mean.rho2.HS <- mean(post.summary.HSGP[4,])
mean.a21.HS<- mean(post.summary.HSGP[5,])
mean.sigma1.HS<- mean(post.summary.HSGP[6,])
mean.sigma2.HS<- mean(post.summary.HSGP[7,])


# True values for reference lines
true_vals <- list(
  beta0  = beta_pts,
  beta1  = beta_marks,
  rho1   = lscale[1],
  rho2   = lscale[2],
  a21    = R[2,1],
  sigma1 = sigma[1],
  sigma2 = sigma[2]
)

# Parameter names and row mapping
params <- c("beta0", "beta1", "rho1", "rho2", "a21", "sigma1", "sigma2")
rows    <- c(1,       2,       3,      4,      5,     6,        7)

# Set up plotting area: 2 rows × 4 cols
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
# Parameter names and row mapping
params_labels <- expression(
  beta[0],
  beta[1],
  phi[1],
  phi[2],
  rho,
  sigma[1],
  sigma[2]
)
rows    <- c(1, 2, 3,4,5,6,7)

# Set up plotting area: 2 rows × 4 cols
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))

# Loop through parameters
for (i in seq_along(params)) {
  p <- params_labels[i]
  r <- rows[i]
  ylim_val <- NULL
  
  hist(Rhat.HSGP[r, ],
          main = p,
          ylab = "Distribution",
       xlab='Rhat',
          )
}

# Fill last (8th) panel with empty plot
plot.new()



# switch the Rhat statistics into ggplot2 version: 
params <- c("beta0", "beta1", "rho1", "rho2", "a21", "sigma1", "sigma2")
rows   <- c(1, 2, 3, 4, 5, 6, 7)

param_data <- lapply(rows, function(r) Rhat.HSGP[r, ])
param_df <- do.call(rbind, param_data)

# Add parameter column
param_df <- data.frame(
  value = unlist(param_data),
  param = rep(params, each = ncol(Rhat.HSGP))
)

# Custom colors
cols <- c(
  #beta0   = "#44C8F5",
  #beta1   = "#44C8F5",
  #rho1    = "#B2D235",
  #rho2    = "#B2D235",
  #a21    = "#F7941D",
  #sigma1  = "#C768A9",
  #sigma2  = "#C768A9"
  beta0 = '#D4D4D4', 
  beta1 = '#D4D4D4',
  rho1='#D4D4D4',
  rho2='#D4D4D4',
  a21 ='#D4D4D4',
  sigma1='#D4D4D4',
  sigma2='#D4D4D4'
)

# ggplot
ggplot(param_df, aes(x = param, y = value, fill = param)) +
  geom_boxplot(outlier.shape = NA, alpha=.5) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(
    labels = c(
      expression(rho),
      expression(beta[0]),
      expression(beta[1]),
      expression(phi[1]),
      expression(phi[2]),
      expression(sigma[1]),
      expression(sigma[2])
    )
  ) +
   scale_y_continuous(
     limits = c(1, 1.01),
     breaks = c(1.00, 1.01)
   ) +
  labs(y = "Rhat posterior", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.position = "none",
    plot.margin   = margin(t = 20, r = 20, b = 60, l = 20)
)




 
 
 
 



# save the figure into a pdf:

# convert the boxplot into ggplot2 -----
# Prepare data in long format for ggplot
params <- c("beta0", "beta1", "rho1", "rho2", "a21", "sigma1", "sigma2")
rows   <- c(1, 2, 3, 4, 5, 6, 7)

param_data <- lapply(rows, function(r) post.summary.HSGP[r, ] )
param_df <- do.call(rbind, param_data)


# True values (convert list to data frame)
true_df <- data.frame(
  param = params,
  value = sapply(params, function(p) true_vals[[p]])
)

# Add parameter column
param_df <- data.frame(
  value = unlist(param_data),
  param = rep(params, each = ncol(post.summary.HSGP))
)


# adjusted parameters:
param_adjust_df <- data.frame(
  value = unlist(param_data)-rep(true_df$value, each=ncol(post.summary.HSGP)),
  param = rep(params, each = ncol(post.summary.HSGP))
)



# Custom colors
cols <- c(
  beta0 = '#D4D4D4', 
  beta1 = '#D4D4D4',
  rho1='#D4D4D4',
  rho2='#D4D4D4',
  a21 ='#D4D4D4',
  sigma1='#D4D4D4',
  sigma2='#D4D4D4'
  )

# ggplot
p_S3_recovery <- ggplot(param_df, aes(x = param, y = value, fill = param)) +
  geom_boxplot(outlier.shape = NA, alpha=.5) +
  geom_point(data = true_df, aes(x = param, y = value),
             shape = 8, size = 4, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(
    labels = c(
      expression(rho),
      expression(beta[0]),
      expression(beta[1]),
      expression(phi[1]),
      expression(phi[2]),
      expression(sigma[1]),
      expression(sigma[2])
    )
  ) +
  labs(y = "Posterior mean distribution", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.position = "none",
    plot.margin   = margin(t = 20, r = 20, b = 60, l = 20)
  )

print(p_S3_recovery)
# the adjusted boxplot for the parameter recovery: 

# ggplot
p_S3_adjust_recovery <- ggplot(param_adjust_df, aes(x = param, y = value, fill = param)) +
  geom_boxplot(outlier.shape = NA, alpha=.5) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(
    labels = c(
      expression(rho),
      expression(beta[0]),
      expression(beta[1]),
      expression(phi[1]),
      expression(phi[2]),
      expression(sigma[1]),
      expression(sigma[2])
    )
  ) +
  labs(y = "posterior-true difference", x=NULL) +
  theme_classic() +
  theme(
    #axis.text.x = element_blank(),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.position = "none",
    plot.margin = margin(t = 40, r = 40, b = 40, l = 40)
)


print(p_S3_adjust_recovery)


 

# the correlation decay plot in S1 (take one as an example): ----
my_cols_mark <- c(
  "#FDE0DD",  # Very Low (light rose)
  "#FBB4B9",  # Low (soft pink-red)
  "#F768A1",  # Medium (bright pinkish-red)
  "#C51B8A",  # High (deep magenta-red)
  "#7A0177"   # Peak (dark crimson-purple)
)
# the viz for Y and M 
viz.S1<-read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/Figure-S1-simulated-data.csv')
viz.S1.df<-as.data.frame(viz.S1)


# iter<- 7
# true parameters: -----
# 
#beta_pts   <- 1      # intercept for points (eta1)
#beta_marks <- -1     # intercept for marks (eta2)
#offset     <- 0.05   # same small offset used in Stan
#grid_res <- 22
# [-1, 1]*[-1,1]
sigma <- c(1, 1)
R <- cbind(c(1, 0.5),
           c(0.5, 1))
lscale <- c(0.5, 0.3)
Astar <- t(chol(R))
A <- diag(sigma) %*% Astar


my_cols <- c(
  "#EAC0BD",  # Very Low
  "#EC6756",  # Medium
  "#D42D24",  # Medium-High
  "#BC1119",  # High
  "#9C0824"   # Peak
)




my_cols_pts <- c(
 # "#B4D4DA",  # Medium
#  "#3F8BBA",  # Medium-High
 # "#26456E"   # Peak
  '#D0D0D0',
  '#7A7A7A',
  '#1E1E1E'
  
)


# viz by the continuous scale: 
ggplot(viz.S1.df, aes(x = x, y = y, color = points)) +
  geom_point(size=8, shape=15) +
  labs(y = "", x = NULL) +
  theme_classic() +
  scale_color_gradientn(colors = my_cols_pts,   # use your custom palette as a gradient
                        #breaks = seq(0, 2, by = 1),
                        #labels = seq(0, 2, by = 1),
                        breaks = seq(0, 2, by = 1),
                        labels = seq(0, 2, by = 1),
                        limits = c(0, 2), 
                        guide = guide_colorbar(
                          barwidth = 2,   # bar width
                          barheight = 20,  # bar height
                          title = "Cell\nCounts",   # <-- this sets the legend title
                          guide = guide_colorbar(
                            barwidth = 2,
                            barheight = 20,
                            title.vjust = 8
                          ))) +
  theme(
    #axis.text.x = element_text(size = 45),
    #axis.text.y = element_text(size = 45),
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 45)
)+ theme(
  legend.title = element_text(size = 40, margin = margin(b = 20)),  # add bottom margin
  legend.text  = element_text(size = 40)
)




# The avg gene expression: 
# the avg gene expression:
viz.S1.df$avg_mark<- viz.S1.df$marks/(viz.S1.df$points+0.05)


ggplot(viz.S1.df, aes(x = x, y = y, color = avg_mark)) +
  geom_point(size=8, shape=15) +
  labs(y = "", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 45)
  )+
  scale_color_gradientn(colors = my_cols,   # use your custom palette as a gradient
                        breaks = seq(0, 8, by = 2),
                        labels = seq(0, 8, by = 2),
                        limits = c(0, 8),
                        guide = guide_colorbar(
                          barwidth = 2,   # bar width
                          barheight = 20,  # bar height
                          title = "Avg\nExpression",   # <-- this sets the legend title
                          guide = guide_colorbar(
                            barwidth = 2,
                            barheight = 20,
                            title.vjust = 8
                          ))) +theme(
                            legend.title = element_text(size = 40, margin = margin(b = 20)),  # add bottom margin
                            legend.text  = element_text(size = 40)
                          )




# rearrange for the viz----
# Calculate average gene expression
viz.S1.df$avg_mark <- viz.S1.df$marks / (viz.S1.df$points + 0.05)

# ===== Plot 1: Cell Counts (WITH axes and ticks) =====
S_a1 <- ggplot(viz.S1.df, aes(x = x, y = y, color = points)) +
  geom_point(size = 6, shape = 15) +
  labs(
    title = "Cell Counts",
    x = "X Coordinate (mm)", 
    y = "Y Coordinate (mm)"
  ) +
  theme_classic() +
  scale_color_gradientn(
    colors = my_cols_pts,
    breaks = seq(0, 2, by = 1),
    labels = seq(0, 2, by = 1),
    limits = c(0, 2), 
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 15,
      title = ""
    )
  ) +
  theme(
    axis.text.x = element_text(size = 11),      # Show x-axis text
    axis.text.y = element_text(size = 11),      # Show y-axis text
    axis.title.x = element_text(size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, margin = margin(r = 10)),
    axis.ticks = element_line(),                # Show ticks
    legend.title = element_text(size = 14, margin = margin(b = 10)),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  coord_fixed()

# ===== Plot 2: Average Expression (NO axes labels/ticks, but keep axes) =====
S_a2 <- ggplot(viz.S1.df, aes(x = x, y = y, color = avg_mark)) +
  geom_point(size = 6, shape = 15) +
  labs(
    title = "Average Expression",
    x = "",    # Empty x label
    y = ""     # Empty y label
  ) +
  theme_classic() +
  scale_color_gradientn(
    colors = my_cols,
    breaks = seq(0, 8, by = 2),
    labels = seq(0, 8, by = 2),
    limits = c(0, 8),
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 15,
      title = ""
    )
  ) +
  theme(
    axis.text.x = element_blank(),              # Hide x-axis text
    axis.text.y = element_blank(),              # Hide y-axis text
    axis.title.x = element_blank(),             # Hide x-axis title
    axis.title.y = element_blank(),             # Hide y-axis title
    axis.ticks = element_blank(),               # Hide ticks
    axis.line = element_line(),                 # Keep axis lines
    legend.title = element_text(size = 14, margin = margin(b = 10)),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  coord_fixed()

# ===== Combine plots in 1x2 layout =====
p_combined <- plot_grid(
  p1, p2,
  align = 'hv',
  axis = 'tb',
  ncol = 2,
  rel_widths = c(1, 1)
)

print(p_combined)




# Plot for the correlation decay: 


# read the RDS file 
#cmdstan_fit_HS_decay <-readRDS('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/S3-M2-cmdstan_fit.rds')
cmdstan_fit_HS_decay <- readRDS(
  "C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/S1-M2-cmdstan_fit.rds"
)


# Now both work:
fit_summary_decay <- cmdstan_fit_HS_decay$summary(
  variables = c("beta0", "beta1", "ell[1]", "ell[2]",
                "Astar[2,1]", "sigma[1]", "sigma[2]"),
  ~c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail")
)



summary_gene_22 <- cmdstan_fit_HS_decay$draws(format = "df")


# caluclate the range parameter: 

fit_summary_decay <- cmdstan_fit_HS_decay$summary(variables = c("beta0", "beta1",   "ell[1]", 'ell[2]', 
                                                       'Astar[2,1]', 'sigma[1]','sigma[2]'), 
                                         c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary_decay


summary_gene_22 <- cmdstan_fit_HS_decay$draws(format = "df")


# the kernel: 
matern32 <- function(d, sigma, rho) {
  sigma^2*(1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}

L1<- 1
L2<- 1

max_dist<- 2*sqrt( L1^2 + L2^2 )
dist_choice<- seq(0, max_dist, by=0.1)




# mean and quantiles of phi1 and phi2

# corr decay plot -----
dist_choice<- seq(0, max_dist, by=0.01)

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


# the true ranges for two latent fields: 
post_rho1_true<- rep(NA, length(dist_choice))
post_rho2_true<- rep(NA, length(dist_choice))
for(k in 1:length(dist_choice)){
  post_rho1_true[k]<- matern32(d=dist_choice[k], sigma=1, rho= lscale[1])
  post_rho2_true[k]<- (A[2,1]^2*matern32(d=dist_choice[k], sigma=1, rho= lscale[1])+
                          A[2,2]^2*matern32(d=dist_choice[k], sigma=1, rho= lscale[2]))/(A[2,1]^2+A[2,2]^2)
}
# find the length where the length of scale==0:
true_range_1<-dist_choice[which(abs(post_rho1_true-0.05)==min(abs(post_rho1_true-0.05)))]
true_range_2<-dist_choice[which(abs(post_rho2_true-0.05)==min(abs(post_rho2_true-0.05)))]





df_rho1<- cbind( dist_choice ,post_mean_rho1, post_upper_rho1, post_lower_rho1, post_rho1_true)
colnames(df_rho1)<- c('dist_choice',  'post_mean_rho1', 'post_upper_rho1', 'post_lower_rho1', 'post_true_rho1')
est_range_1<-dist_choice[which(abs(df_rho1[,2]-0.05)==min(abs(df_rho1[,2]-0.05)))]


eta1_corr<-ggplot(data=df_rho1, aes(x=dist_choice))+
  geom_ribbon(aes(ymin=post_lower_rho1, ymax= post_upper_rho1), fill = "#cbcfd2")+
  geom_line(aes(y=post_mean_rho1), color='#000', linewidth=2)+
  geom_line(aes(y=post_true_rho1), color='#ED1B2F', linewidth=2)+
  geom_hline(yintercept=0.05, linetype='dashed', col='#000', linewidth=2, alpha=.5)+
  xlim(0, max_dist)+
  #geom_vline(xintercept=est_range_1, linetype='dashed', col='#44C8F5', linewidth=1.05, alpha=.5)+
  #geom_vline(xintercept=true_range_1, linetype='dashed', col='#F7941D', linewidth=1.05, alpha=.5)+
  labs(x='Distance (Millimeter)', y=expression( 'Posterior distribution of '~ r[1]) )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -0, hjust = 0))+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.position = "none"
  )


df_rho2<- cbind( dist_choice ,post_mean_rho2, post_upper_rho2, post_lower_rho2, post_rho2_true)
colnames(df_rho2)<- c('dist_choice',  'post_mean_rho2', 'post_upper_rho2', 'post_lower_rho2', 'post_true_rho2')
est_range_2<-dist_choice[which(abs(df_rho2[,2]-0.05)==min(abs(df_rho2[,2]-0.05)))]


eta2_corr<- ggplot(data=df_rho2, aes(x=dist_choice))+
  geom_ribbon(aes(ymin=post_lower_rho2, ymax= post_upper_rho2), fill = "#cbcfd2")+
  geom_line(aes(y=post_mean_rho2), color='#000', linewidth=2)+
  geom_line(aes(y=post_true_rho2), color='#ED1B2F', linewidth=2)+
  #geom_vline(xintercept=est_range_2, linetype='dashed', col='#44C8F5', linewidth=1.05, alpha=.5)+
  #geom_vline(xintercept=true_range_2, linetype='dashed', col='#F7941D', linewidth=1.05, alpha=.5)+
  xlim(0, max_dist)+
  geom_hline(yintercept=0.05, linetype='dashed', col='black', linewidth=2, alpha=.5)+
  labs(x='Distance (Millimeter)', y=expression( 'Posterior distribution of '~ r[2]))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -0, hjust = 0))+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.position = "none"
  )


S3_right_side <- plot_grid(eta1_corr, eta2_corr)



# the other check is the true correlation decay plots: 



S3_subfig<- plot_grid(p_S3_recovery,S3_right_side )





# Summarize all figures together: -------

# check the figure list: 


# A
fig1<- S_a1
fig2<- S_a2

#B 
fig3<- eta1_corr
fig4<- eta2_corr

# C 
fig5<- S1_basis
fig6<- S1_basis_t

# D
fig7<- S1_bound
fig8<- S1_bound_t


# E
fig9<- S2_mpp
fig10<- p_gridded_cells
fig11<- p_gridded_avg_express

# F
fig12<- S2_p1
fig13<- S2_p2
# Panel A: figs 1 & 2
panel_A <- plot_grid(fig1, fig2,
                     ncol = 2,
                     align = "hv",
                     rel_widths = c(1, 1))

# Panel B: figs 3 & 4
panel_B <- plot_grid(fig3, fig4,
                     ncol = 2,
                     align = "hv")

# Panel C: figs 5 & 6
panel_C <- plot_grid(fig5, fig6,
                     ncol = 2,
                     align = "hv")

# Panel D: figs 7 & 8
panel_D <- plot_grid(fig7, fig8,
                     ncol = 2,
                     align = "hv")

# Panel E: figs 9, 10, 11
panel_E <- plot_grid(fig9, fig10, fig11,
                     ncol = 3,
                     align = "hv")

# Panel F: figs 12, 13
panel_F <- plot_grid(fig12, fig13,
                     ncol = 2,
                     align = "hv")

## ---- Add labels A–F -------------------------------------------------------
panel_A <- ggdraw() + draw_label("A", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_A)

panel_B <- ggdraw() + draw_label("B", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_B)

panel_C <- ggdraw() + draw_label("C", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_C)

panel_D <- ggdraw() + draw_label("D", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_D)

panel_E <- ggdraw() + draw_label("E", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_E)

panel_F <- ggdraw() + draw_label("F", x = 0, y = 1, hjust = -0.5, vjust = 1.5, size = 14) +
  draw_plot(panel_F)

## ---- Final 3 × 2 layout ---------------------------------------------------
final_plot <- plot_grid(
  plot_grid(panel_A, panel_C, ncol = 2, rel_widths = c(1, 1)),
  plot_grid(panel_B, panel_D, ncol = 2, rel_widths = c(1, 1)),
  plot_grid(panel_E, panel_F, ncol = 2, rel_widths = c(1, 1)),
  ncol = 1,
  rel_heights = c(1, 1, 1)
)


ggsave("C:/Users/xumc7/OneDrive/Desktop/BIOS702 Protocol/Protocol/final_figure.png", final_plot,
       width = 12, height = 16, dpi = 600)





# the simulated-based correlation distribution （Corr(Y, M/Y)） ----
file_path <- 'C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Simulation3/grid-grid-result/'


# the summarized version of the dataset
true_corr<- rep(NA, 100)
post_mean_corr<- rep(NA, 100)


# read csv 
for (i in c(1:100)){
  file_df<- read.csv(paste0(file_path, 'A_mat_compare_',i,'.csv'))
  # lambda1:
  true_corr[i] <- file_df[,4][1]
}



avg_true_corr<-mean(true_corr)

for (i in c(1:100)){
   load(paste0(file_path, 'corr_dist_',i,'.Rdata'))
  # lambda1:
  post_mean_corr[i] <- mean(corr_dist_1)
}


post_mean_corr

par(mfrow = c(1, 1))

# Density plot
plot(density( na.omit(post_mean_corr)), 
     main = "Density of Posterior Mean Correlation",
     xlab = "Posterior Mean Correlation",
     lwd = 2)

# Add vertical line for the true value
abline(v = avg_true_corr, col = "red", lwd = 2, lty = 2)







# Appendix codes in Simulation 2-----

# for the gridded gene expression:
par(mfrow = c(1,1), mar = c(4,4,4,5))
image.plot(
  x = x_coords,
  y = y_coords,
  M_high ,
  breaks = make_breaks(M_high, 10),
  col = viridis(10+2, option = "magma")[-c(1,2)],
  main = "Gridded Gene Expression",
  xlab = "",
  ylab = '',
  useRaster = F, 
  xaxt = "n",
  yaxt = "n",
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
  legend.width = 1,
  cex.main=3,
  cex.axis = 2.5 
)
grid(nx = ncol(M_high), ny = nrow(M_high), col = "black", lwd = 0.5)




# for the gene expression in Figure A: 

ggplot(viz.S1.df, aes(x = x, y = y, color = marks)) +
  geom_point(size=8, shape=15) +
  labs(y = "", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 45)
  )+
  scale_color_gradientn(colors = my_cols,   # use your custom palette as a gradient
                        breaks = seq(0, 8, by = 2),
                        labels = seq(0, 8, by = 2),
                        limits = c(0, 8),
                        guide = guide_colorbar(
                          barwidth = 2,   # bar width
                          barheight = 20,  # bar height
                          title = "Gene\nExpression",   # <-- this sets the legend title
                          guide = guide_colorbar(
                            barwidth = 2,
                            barheight = 20,
                            title.vjust = 8
                          ))) +theme(
                            legend.title = element_text(size = 40, margin = margin(b = 20)),  # add bottom margin
                            legend.text  = element_text(size = 40)
                          )

