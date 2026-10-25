## ---- CONFIGURE ME ------------------------------------------------------------
## Sourced by 02_experiment_FPR.R and 04_experiment_TPR.R, which set REPO_ROOT
## first. This file reads no data and writes no output, so it needs only REPO_ROOT.
if (!exists("REPO_ROOT")) REPO_ROOT <- "."   # root of the DenMark checkout (holds inst/stan/)
## -----------------------------------------------------------------------------

## ============================================================
## 00_functions.R
## Core utilities for the density <-> gene correlation study.
##
## Key idea: in your model the "correlation between cell density and
## gene expression" is rho = Astar[2,1] (the off-diagonal of the
## Cholesky factor of the 2x2 correlation matrix). For a 2x2 corr matrix
## R = Astar %*% t(Astar), we have R[2,1] = Astar[2,1], so Astar[2,1]
## IS the correlation. We simulate data so that the latent-field
## correlation equals a chosen rho, fit your (fixed) Stan model, and
## flag a gene as a "density-correlated gene" (DCG) when 0 lies OUTSIDE
## the 95% posterior CI of Astar[2,1].
## ============================================================

library(cmdstanr)    # install.packages("cmdstanr", repos=c("https://stan-dev.r-universe.dev"))
library(posterior)


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

STAN_FILE <- file.path(REPO_ROOT, "inst", "stan", "DenMark.stan")   # your fixed Stan file (in working dir)

## Compile once, reuse everywhere.
compile_model <- function(stan_file = STAN_FILE) cmdstan_model(stan_file)

## ---- 1. spatial grid (centred at 0 -- required by the HSGP basis) ----
make_grid <- function(n_per_side = 12, half_width = 0.5) {
  s <- seq(-half_width, half_width, length.out = n_per_side)
  g <- as.matrix(expand.grid(x = s, y = s))
  colnames(g) <- c("x", "y")
  g                                      # n = n_per_side^2 rows
}

## ---- 2. Matern 3/2 unit-variance kernel ----
matern32 <- function(D, ell) {
  s3 <- sqrt(3) * D / ell
  (1 + s3) * exp(-s3)
}

## ---- 3. simulate ONE marked point pattern (cells + one gene) -------
## This mirrors the model's generative process exactly:
##   z1 = sigma1 * w1
##   z2 = sigma2 * (rho * w1 + sqrt(1 - rho^2) * w2)   <- co-regionalisation
##   y1 ~ Poisson(exp(log_grid_area + beta0 + z1))     <- cell density
##   y2 ~ Poisson(exp(beta1 + log(y1 + 0.05) + z2))    <- gene expression
## so the latent correlation between z1 and z2 equals `rho` = Astar[2,1].
simulate_mpp <- function(coords, rho,
                         beta0 = 2, beta1 = 0,
                         sigma = c(0.6, 0.6),
                         ell   = c(0.25, 0.25),
                         log_grid_area = 0,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n  <- nrow(coords)
  D  <- as.matrix(dist(coords))
  U1 <- chol(matern32(D, ell[1]) + diag(1e-8, n))   # upper-triangular
  U2 <- chol(matern32(D, ell[2]) + diag(1e-8, n))
  w1 <- as.vector(crossprod(U1, rnorm(n)))          # ~ N(0, K(ell1)), var 1
  w2 <- as.vector(crossprod(U2, rnorm(n)))          # ~ N(0, K(ell2)), var 1
  z1 <- sigma[1] * w1
  z2 <- sigma[2] * (rho * w1 + sqrt(1 - rho^2) * w2)
  y1 <- rpois(n, exp(log_grid_area + beta0 + z1))
  y2 <- rpois(n, exp(beta1 + log(y1 + 0.05) + z2))
  list(coords = coords, y1 = y1, y2 = y2, z1 = z1, z2 = z2,
       rho = rho, log_grid_area = log_grid_area,
       params = list(beta0 = beta0, beta1 = beta1, sigma = sigma, ell = ell))
}

## ---- 4. simulate ONE tissue (shared cell field) + many genes -------
## Used for the FDR experiment: one cell-density field, many genes, each
## with its own target rho (positive / negative / zero).
simulate_tissue <- function(coords, gene_rho,
                            beta0 = 2, beta1 = 0,
                            sigma = c(0.6, 0.6),
                            ell   = c(0.25, 0.25),
                            log_grid_area = 0,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords); G <- length(gene_rho)
  D  <- as.matrix(dist(coords))
  U1 <- chol(matern32(D, ell[1]) + diag(1e-8, n))
  U2 <- chol(matern32(D, ell[2]) + diag(1e-8, n))
  w1 <- as.vector(crossprod(U1, rnorm(n)))          # the ONE shared density field
  z1 <- sigma[1] * w1
  y1 <- rpois(n, exp(log_grid_area + beta0 + z1))
  Y2 <- matrix(NA_integer_, n, G)
  for (g in seq_len(G)) {
    w2 <- as.vector(crossprod(U2, rnorm(n)))
    z2 <- sigma[2] * (gene_rho[g] * w1 + sqrt(1 - gene_rho[g]^2) * w2)
    Y2[, g] <- rpois(n, exp(beta1 + log(y1 + 0.05) + z2))
  }
  list(coords = coords, y1 = y1, Y2 = Y2, gene_rho = gene_rho,
       log_grid_area = log_grid_area)
}

## ---- 5. HSGP multi-indices (1:M in each of d dims) -----------------
build_indices <- function(M, d = 2) {
  idx <- as.matrix(expand.grid(rep(list(1:M), d)))
  dimnames(idx) <- NULL
  storage.mode(idx) <- "integer"
  idx                                    # mstar x d, with mstar = M^d
}

## ---- 6. assemble the data list the Stan model expects -------------
make_stan_data <- function(coords, y1, y2, log_grid_area = 0,
                           M = 5, c_factor = 1.5,
                           mrange = 0.3, is_centered_PHI = 0) {
  d <- ncol(coords)
  cc <- coords; dimnames(cc) <- NULL
  S  <- apply(abs(cc), 2, max)
  L  <- c_factor * S
  indices <- build_indices(M, d)
  storage.mode(y1) <- "integer"
  storage.mode(y2) <- "integer"
  list(
    d                = d,
    L                = as.array(as.numeric(L)),
    mstar            = nrow(indices),
    n                = nrow(cc),
    coords           = cc,
    log_grid_area    = log_grid_area,
    y1               = y1,
    y2               = y2,
    mrange           = mrange,
    indices          = indices,
    is_centerted_PHI = as.integer(is_centered_PHI)   # name matches your Stan exactly
  )
}

## ---- 7. fit + extract the rho = Astar[2,1] posterior --------------
## detected == TRUE  <=>  0 is OUTSIDE the 95% credible interval.
fit_rho <- function(stan_data, mod,
                    chains = 2, parallel_chains = 1,
                    iter_warmup = 400, iter_sampling = 400,
                    adapt_delta = 0.9, seed = 1, prob = 0.95) {
  fit <- mod$sample(
    data = stan_data, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  rho <- as.numeric(fit$draws("Astar[2,1]", format = "draws_matrix"))
  a   <- (1 - prob) / 2
  q   <- quantile(rho, c(a, 1 - a))
  dg  <- fit$summary("Astar[2,1]", "rhat", "ess_bulk")
  list(
    rho_mean  = mean(rho),
    rho_lower = unname(q[1]),
    rho_upper = unname(q[2]),
    detected  = !(q[1] <= 0 & q[2] >= 0),
    sign      = sign(mean(rho)),
    rhat      = dg$rhat,
    ess_bulk  = dg$ess_bulk,
    n_diverg  = sum(fit$diagnostic_summary("divergences")$num_divergent)
  )
}

## ---- 8. one full FPR replicate: simulate (own field) -> fit -> test
run_replicate <- function(coords, rho, sim_seed, fit_seed, mod,
                          sim_args = list(), stan_args = list(), fit_args = list()) {
  sim <- do.call(simulate_mpp,
                 c(list(coords = coords, rho = rho, seed = sim_seed), sim_args))
  sd  <- do.call(make_stan_data,
                 c(list(coords = sim$coords, y1 = sim$y1, y2 = sim$y2,
                        log_grid_area = sim$log_grid_area), stan_args))
  res <- do.call(fit_rho, c(list(stan_data = sd, mod = mod, seed = fit_seed), fit_args))
  res$true_rho <- rho
  res$sim_seed <- sim_seed
  res
}
