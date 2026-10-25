## ---- CONFIGURE ME ------------------------------------------------------------
## Launched by cluster-long-whole-DCIS.sh and cluster-long-whole-Tumor.sh, which
## export the three variables below; the defaults apply when they are not set.
DATA_DIR   <- Sys.getenv("DATA_DIR", "data")       # holds Dataset/xenium-breast-cancer/rep1/
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", "output")   # reads qc-whole/ (qc-select-genes-whole.R); writes BC-whole/
REPO_ROOT  <- Sys.getenv("REPO_ROOT", ".")         # root of the DenMark checkout
STAN_DIR   <- file.path(REPO_ROOT, "CodeInPaper", "Xenium-BreastCancer", "whole-section")   # holds M2-latent-field-mod-genes.stan, not in this repository (see "Missing files" in CodeInPaper/Xenium-BreastCancer/README.md)
## -----------------------------------------------------------------------------

####################################################
# DenMark general framework - WHOLE TISSUE by group
# ---------------------------------------------------
# Runs DenMark on all cells of one cell-type group,
# over the whole tissue (no ROI box).
# args: [1] array task = index into this group's QC list
#       [2] group id (1=Immune, 2=DCIS, 3=Invasive, 4=Tumor)
#
# Division B (fig B): groups 1, 2, 3  (Immune / DCIS / Invasive)
# Division C (fig C): groups 1, 4     (Immune / Tumor)
####################################################

args = commandArgs(trailingOnly = TRUE)

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

if (length(args) < 2){
  stop("Two arguments required: array task and group id.", call.=FALSE)
}

# quantiles:
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)

options(cmdstanr_max_rows = 20)

task     <- as.numeric(args[1])   # index into the group's QC-passed gene list
group_id <- as.numeric(args[2])   # 1=Immune, 2=DCIS, 3=Invasive, 4=Tumor

# --- the 4 unique cell-type groups -------------------
groups <- list(
  Immune   = c("CD8+_T_Cells"),
  DCIS     = c("DCIS_1", "DCIS_2"),
  Invasive = c("Invasive_Tumor"),
  Tumor    = c("DCIS_1", "DCIS_2", "Invasive_Tumor")
)
if (group_id < 1 || group_id > length(groups))
  stop(sprintf("group id must be 1..%d", length(groups)), call.=FALSE)

grp_name <- names(groups)[group_id]
members  <- groups[[group_id]]
cat("Running group:", grp_name, "(", paste(members, collapse = ", "), ")\n")

# --- map array task -> gene id via QC-passed list ----
qc <- readRDS(paste0(file.path(OUTPUT_DIR, 'qc-whole/'),
                     grp_name, '_qc_genes.rds'))
if (task > qc$n_pass) {
  message(sprintf("task %d exceeds %d QC-passed genes for %s - exiting.",
                  task, qc$n_pass, grp_name))
  quit(save = "no", status = 0)
}
iter <- qc$gene_ids[task]         # true gene id; expression column = 1 + iter
print(task); print(iter); print(grp_name)


# [Step 0]: load the data -----------------------------------
cell_gene    <- read.csv(file.path(DATA_DIR, 'Dataset/xenium-breast-cancer/rep1/cell_by_gene_expression.csv'))
cell_spatial <- read.csv(file.path(DATA_DIR, 'Dataset/xenium-breast-cancer/rep1/cell_spatial.csv'))
cell_spatial$x <- cell_spatial$x / 1000
cell_spatial$y <- cell_spatial$y / 1000


# [Step 0C]: select the group's cells (whole tissue) --------
sel_idx <- which(cell_spatial$celltype %in% members)

center_slice_cell_location <- cell_spatial[sel_idx, ]
gene_matrix                <- cell_gene[sel_idx, ]

if (nrow(center_slice_cell_location) < 20) {
  message(sprintf("Only %d cells for %s - skipping.",
                  nrow(center_slice_cell_location), grp_name))
  quit(save = "no", status = 0)
}

center_slice_cell_location$Gene <- gene_matrix[, (1 + iter)]


# 0D. hyperparameters ---------------------------------------
grid_size    <- 50
basis_choice <- 25
boundary_choice <- 1.5


# [Step 1]: window normalization (group's spatial extent) ---
L1 <- max(abs(center_slice_cell_location$x -
               (min(center_slice_cell_location$x) + max(center_slice_cell_location$x)) / 2))
L2 <- max(abs(center_slice_cell_location$y -
               (min(center_slice_cell_location$y) + max(center_slice_cell_location$y)) / 2))

coords <- unname(as.matrix(expand.grid(x = seq(-L1, L1, length.out = grid_size),
                                       y = seq(-L2, L2, length.out = grid_size))))
distMat <- fields::rdist(coords)
mrange  <- max(distMat) / (2 * 2.75)
grid_area <- (2 * L1 * 2 * L2) / (grid_size^2)
log_grid_area <- log(grid_area)


# [Step 2]: gridded N and M ---------------------------------
loc <- center_slice_cell_location
eg_gene_express <- center_slice_cell_location$Gene

init_x <- min(loc$x); init_y <- min(loc$y)
step_x <- (max(loc$x) - min(loc$x)) / grid_size
step_y <- (max(loc$y) - min(loc$y)) / grid_size

loc_100   <- matrix(NA, nrow = grid_size, ncol = grid_size)
marks_100 <- matrix(NA, nrow = grid_size, ncol = grid_size)
for (i in 1:grid_size){
  for (j in 1:grid_size){
    x_bound_1 <- init_x + (i - 1) * step_x
    x_bound_2 <- init_x + i * step_x
    y_bound_1 <- init_y + (j - 1) * step_y
    y_bound_2 <- init_y + j * step_y

    pts_index <- intersect(which(loc$x >= x_bound_1 & loc$x <= x_bound_2),
                           which(loc$y >= y_bound_1 & loc$y <= y_bound_2))
    marks_100[i, j] <- sum(eg_gene_express[pts_index])
    loc_100[i, j]   <- length(loc[pts_index, ][, 1])
  }
}

Y <- as.vector(loc_100)
M <- as.vector(marks_100)


# [Step 3]: HSGP preparation --------------------------------
xRangeDat <- c(-L1, L1)
yRangeDat <- c(-L2, L2)

m1 <- basis_choice; m2 <- basis_choice; mstar <- m1 * m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(boundary_choice, boundary_choice)
L <- c * Lstar
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[, 2:1]))

input <- list(n = length(Y),
              y1 = Y,
              y2 = M,
              d = 2,
              mstar = mstar,
              coords = coords,
              log_grid_area = log_grid_area,
              mrange = mrange,
              indices = S,
              L = L,
              is_centerted_PHI = 0)

init_fun <- function() list(
  Astar = matrix(c(1, 0, 0.5, 1), 2, 2),
  sigma = c(1, 1),
  beta0 = 0,
  beta1 = -1,
  ell = c(1, 1),
  betab1 = rep(0, m1 * m2),
  betab2 = rep(0, m1 * m2)
)

stan_file <- file.path(STAN_DIR, 'M2-latent-field-mod-genes.stan')
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
cmdstan_fit_HS <- mod$sample(data = input,
                             chains = 4,
                             parallel_chains = 4,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             adapt_delta = 0.95,
                             init = init_fun,
                             max_treedepth = 10,
                             step_size = 1)

elapsed_time <- cmdstan_fit_HS$time()
elapsed_time$total / 60

fit_summary_HS <- cmdstan_fit_HS$summary(
  variables = c("beta0", "beta1", "ell[1]", 'ell[2]',
                'Astar[2,1]', 'sigma[1]', 'sigma[2]'),
  c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))

save_df <- as.data.frame(fit_summary_HS)


# [Step 4]: write csv into a per-group folder ---------------
out_dir <- paste0(file.path(OUTPUT_DIR, 'BC-whole/BC-whole-'), grp_name, '/')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(save_df,
          paste0(out_dir, colnames(gene_matrix)[1 + iter], '_M2_', iter, '.csv'))
