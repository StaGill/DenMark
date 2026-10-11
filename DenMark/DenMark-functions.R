###############################################################################
#                                                                             #
#   DenMark: Density-dependent Marked point process framework (R version)     #
#                                                                             #
#   User-facing functions to grid single-cell resolution spatial              #
#   transcriptomics (ST) data and fit the DenMark models with Stan            #
#   (cmdstanr backend).                                                       #
#                                                                             #
#   Models (Stan files shipped alongside this script):                        #
#     M2  'DenMark-M2.stan'  full DenMark, two dependent latent fields        #
#           log lambda1(g) = log|g| + beta0 + a11 w1(g)                       #
#           log lambda2(g) = beta1 + log(Y(g)+0.05) + a21 w1(g) + a22 w2(g)   #
#           rho = Astar[2,1] = density-expression correlation                 #
#     M1  'DenMark-M1.stan'  simplified DenMark, independent fields (a21 = 0) #
#                                                                             #
#   Mingchi Xu, Alexandra M. Schmidt, Qihuang Zhang                           #
###############################################################################
#
#  CONTENTS
#  --------------------------------------------------------------------------
#  Data preparation
#    denmark_filter_genes()   QC filter (variance + detection rate)
#    buildgridpp()            grid cells: cell counts Y and expression M
#    grid_pts_centroid()      grid centroids, grid area, number of grids
#    denmark_prepare()        build the Stan data list (HSGP inputs)
#  Model fitting
#    denmark_model()          compile (and cache) M1 / M2
#    denmark_init()           valid initial values for M1 / M2
#    DenMark()                fit DenMark to one gene  -> object of class "denmark"
#    DenMarkmain()            backward-compatible wrapper of the old interface
#    denmark_scan()           fit many genes, call PCG / NCG (DCG discovery)
#  Post-processing
#    denmark_summary()        posterior summary of the key parameters
#    denmark_draws()          harmonised posterior draws (rho, a11, a21, a22 ...)
#    denmark_fields()         latent-field decomposition per grid
#    denmark_decay()          spatial correlation decay curves + practical ranges
#    denmark_compare()        WAIC model comparison (e.g. M2 vs M1)
#    denmark_check_hsgp()     rule-of-thumb check of basis number / boundary
#    denmark_save(), denmark_load()  keep fits across R sessions
#  Visualisation (ggplot2; patchwork optional for multi-panel figures)
#    plot_denmark_grid(), plot_denmark_fields(), plot_denmark_decay(),
#    plot_denmark_scan()
#  Utilities
#    denmark_simulate()       simulate a marked point pattern with known rho
#    q2.5(), q25(), q50(), q75(), q97.5(), matern32()
#
#  REQUIREMENTS
#  --------------------------------------------------------------------------
#    cmdstanr (+ a working CmdStan, see https://mc-stan.org/cmdstanr/),
#    posterior, loo, ggplot2.   Optional: patchwork (multi-panel figures).
#
#  QUICK START
#  --------------------------------------------------------------------------
#    source("DenMark/DenMark-functions.R")
#    fit <- DenMark(dataset = df, gene = "CDH1", model = "M2",
#                   grid_size = 10, basis = 25, boundfactor = 1.5)
#    fit                           # print method: rho, WAIC, diagnostics
#    denmark_summary(fit)
#    plot_denmark_fields(fit)
###############################################################################


## ------------------------------------------------------------------------- ##
## 0. Packages and internal state                                            ##
## ------------------------------------------------------------------------- ##

local({
  need <- c("cmdstanr", "posterior", "loo", "ggplot2")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    warning("DenMark: please install the missing package(s): ",
            paste(miss, collapse = ", "),
            "\n  cmdstanr: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
            call. = FALSE)
})

# Internal environment: remembers where this file lives (to find the Stan
# files) and caches compiled models so that a gene scan compiles only once.
.denmark_env <- new.env(parent = emptyenv())
.denmark_env$models <- list()
.denmark_env$dir <- tryCatch({
  f <- sys.frame(1)$ofile                      # set when the file is source()d
  if (is.null(f)) getwd() else dirname(normalizePath(f))
}, error = function(e) getwd())


## ------------------------------------------------------------------------- ##
## 0b. Small utilities                                                       ##
## ------------------------------------------------------------------------- ##

# posterior quantile helpers (kept for compatibility with earlier scripts)
q2.5  <- function(x) stats::quantile(x, probs = 0.025, names = FALSE)
q25   <- function(x) stats::quantile(x, probs = 0.25,  names = FALSE)
q50   <- function(x) stats::quantile(x, probs = 0.50,  names = FALSE)
q75   <- function(x) stats::quantile(x, probs = 0.75,  names = FALSE)
q97.5 <- function(x) stats::quantile(x, probs = 0.975, names = FALSE)

#' Matern (nu = 3/2) covariance, same parameterisation as the Stan models:
#'   C(d) = sigma^2 (1 + sqrt(3) d / ell) exp(-sqrt(3) d / ell)
#' The argument name `rho` (length-scale) is kept for backward compatibility.
matern32 <- function(d, sigma = 1, rho) {
  u <- sqrt(3) * d / rho
  sigma^2 * (1 + u) * exp(-u)
}

# Distance at which a Matern-3/2 correlation drops to `threshold`
# (0.05 -> practical range = 2.7389 * ell).
.matern32_range_factor <- function(threshold = 0.05) {
  u <- stats::uniroot(function(u) (1 + u) * exp(-u) - threshold,
                      c(1e-8, 100), tol = 1e-10)$root
  u / sqrt(3)
}

# Standardise the input data frame to columns x, y, gene.
.denmark_xyz <- function(dataset, gene = NULL, x = "x", y = "y") {
  dataset <- as.data.frame(dataset)
  if (!all(c(x, y) %in% names(dataset))) {
    if (ncol(dataset) < 3)
      stop("`dataset` needs columns '", x, "', '", y, "' and one expression column.")
    warning("Columns '", x, "'/'", y, "' not found; using columns 1 and 2 as x and y.",
            call. = FALSE)
    names(dataset)[1:2] <- c(x, y)
  }
  if (is.null(gene)) {
    others <- setdiff(names(dataset), c(x, y))
    if (length(others) < 1) stop("No expression column found in `dataset`.")
    gene <- others[1]
  }
  if (!gene %in% names(dataset)) stop("Gene '", gene, "' is not a column of `dataset`.")
  out <- data.frame(x = as.numeric(dataset[[x]]),
                    y = as.numeric(dataset[[y]]),
                    gene = as.numeric(dataset[[gene]]))
  ok <- stats::complete.cases(out)
  if (any(!ok)) {
    warning(sum(!ok), " cell(s) with missing coordinates/expression removed.", call. = FALSE)
    out <- out[ok, , drop = FALSE]
  }
  if (any(out$gene < 0)) stop("Expression values must be non-negative counts.")
  attr(out, "gene") <- gene
  out
}

# grid_size may be a single number (same for x and y) or c(nx, ny)
.grid_dims <- function(grid_size) {
  grid_size <- as.integer(round(grid_size))
  if (length(grid_size) == 1) grid_size <- rep(grid_size, 2)
  if (length(grid_size) != 2 || any(grid_size < 2))
    stop("`grid_size` must be one integer >= 2 or a vector c(nx, ny).")
  grid_size
}

.bbox <- function(x, y, bbox = NULL) {
  if (is.null(bbox)) return(c(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)))
  bbox <- as.numeric(bbox)
  if (length(bbox) != 4 || bbox[2] <= bbox[1] || bbox[4] <= bbox[3])
    stop("`bbox` must be c(xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax.")
  stats::setNames(bbox, c("xmin", "xmax", "ymin", "ymax"))
}


## ------------------------------------------------------------------------- ##
## 1. Data preparation                                                       ##
## ------------------------------------------------------------------------- ##

#' Quality-control filter for genes
#'
#' Keeps genes whose variance across cells is above the `var_quantile`
#' quantile AND that are detected (count > 0) in at least `min_nonzero` of
#' the cells. This is the filter used in the DenMark paper
#' (drop the lowest 5% variance genes, require >= 5% non-zero cells).
#'
#' @param expr cells x genes matrix / data.frame of counts (numeric columns only)
#' @return character vector of the genes kept
denmark_filter_genes <- function(expr, var_quantile = 0.05, min_nonzero = 0.05) {
  mat <- as.matrix(expr)
  storage.mode(mat) <- "double"
  if (is.null(colnames(mat))) stop("`expr` must have gene names as column names.")
  gene_var <- apply(mat, 2, stats::var)
  keep_var <- gene_var >= stats::quantile(gene_var, var_quantile, na.rm = TRUE)
  keep_det <- colMeans(mat > 0) >= min_nonzero
  colnames(mat)[keep_var & keep_det]
}


#' Discretise cell locations and expression on a regular grid
#'
#' Each cell is assigned to exactly one grid cell (cells on an interior edge
#' go to the upper/right cell; cells on the outer boundary are kept).
#'
#' @param grid_size number of grid cells per axis, or c(nx, ny)
#' @param dataset data.frame with columns x, y and the expression column(s)
#' @param gene expression column to use (default: first non-coordinate column)
#' @param bbox optional c(xmin, xmax, ymin, ymax); default = data range.
#'   Cells outside the box are dropped.
#' @return a list. The first three elements keep the order of the original
#'   DenMark function, so `res[[1]]`, `res[[2]]`, `res[[3]]` still work:
#'   \item{loc_mark_mat}{N x 2 matrix, columns 'point' (= Y) and 'mark' (= M)}
#'   \item{Ymat}{nx x ny matrix of cell counts (rows = x, columns = y)}
#'   \item{Mmat}{nx x ny matrix of summed expression}
#'   \item{grid}{data.frame: grid_index, i, j, x, y (centroids), Y, M}
#'   \item{bbox, step, grid_area, grid_size, gene}
#'   The vectors are ordered with x varying fastest, i.e. grid_index =
#'   i + (j - 1) * nx, identical to as.vector(Ymat).
buildgridpp <- function(grid_size, dataset, gene = NULL, bbox = NULL) {
  dat <- .denmark_xyz(dataset, gene)
  gene <- attr(dat, "gene")
  dims <- .grid_dims(grid_size); nx <- dims[1]; ny <- dims[2]
  bb <- .bbox(dat$x, dat$y, bbox)

  inside <- dat$x >= bb["xmin"] & dat$x <= bb["xmax"] &
            dat$y >= bb["ymin"] & dat$y <= bb["ymax"]
  if (any(!inside)) dat <- dat[inside, , drop = FALSE]
  if (nrow(dat) == 0) stop("No cells inside the bounding box.")

  dx <- unname(bb["xmax"] - bb["xmin"]) / nx
  dy <- unname(bb["ymax"] - bb["ymin"]) / ny
  if (dx <= 0 || dy <= 0) stop("All cells share the same x or y coordinate.")
  ix <- pmin(pmax(floor((dat$x - bb["xmin"]) / dx) + 1L, 1L), nx)
  iy <- pmin(pmax(floor((dat$y - bb["ymin"]) / dy) + 1L, 1L), ny)
  idx <- ix + (iy - 1L) * nx

  Y <- tabulate(idx, nbins = nx * ny)
  M <- numeric(nx * ny)
  s <- rowsum(dat$gene, idx)
  M[as.integer(rownames(s))] <- s[, 1]
  if (any(abs(M - round(M)) > 1e-8))
    warning("Expression values are not integers; summed counts are rounded for the Poisson model.",
            call. = FALSE)
  M <- round(M)

  gi <- rep(seq_len(nx), times = ny)
  gj <- rep(seq_len(ny), each = nx)
  grid <- data.frame(grid_index = seq_len(nx * ny), i = gi, j = gj,
                     x = unname(bb["xmin"]) + (gi - 0.5) * dx,
                     y = unname(bb["ymin"]) + (gj - 0.5) * dy,
                     Y = Y, M = M)

  loc_mark_mat <- cbind(point = Y, mark = M)
  list(loc_mark_mat = loc_mark_mat,
       Ymat = matrix(Y, nx, ny),
       Mmat = matrix(M, nx, ny),
       grid = grid,
       bbox = bb, step = c(dx = dx, dy = dy), grid_area = dx * dy,
       grid_size = c(nx = nx, ny = ny), gene = gene)
}


#' Grid centroids, grid area and number of grids
#'
#' @return list(grid_points, grid_area, N, step, bbox); the first three
#'   elements keep the order of the original function. grid_points are the
#'   true centroids of the grid cells (x varying fastest).
grid_pts_centroid <- function(grid_size, dataset, bbox = NULL) {
  dataset <- as.data.frame(dataset)
  if (!all(c("x", "y") %in% names(dataset))) names(dataset)[1:2] <- c("x", "y")
  dims <- .grid_dims(grid_size); nx <- dims[1]; ny <- dims[2]
  bb <- .bbox(dataset$x, dataset$y, bbox)
  dx <- unname(bb["xmax"] - bb["xmin"]) / nx
  dy <- unname(bb["ymax"] - bb["ymin"]) / ny
  grid_points <- expand.grid(x = unname(bb["xmin"]) + (seq_len(nx) - 0.5) * dx,
                             y = unname(bb["ymin"]) + (seq_len(ny) - 0.5) * dy)
  list(grid_points = grid_points, grid_area = dx * dy, N = nrow(grid_points),
       step = c(dx = dx, dy = dy), bbox = bb)
}


#' HSGP basis-function index matrix (all pairs 1..basis x 1..basis)
denmark_hsgp_indices <- function(basis) {
  basis <- as.integer(basis)
  if (length(basis) == 1) basis <- rep(basis, 2)
  S <- as.matrix(expand.grid(S1 = seq_len(basis[1]), S2 = seq_len(basis[2])))
  dimnames(S) <- NULL
  storage.mode(S) <- "integer"
  S
}


#' Build the Stan data list for DenMark
#'
#' The region of interest (bounding box) is moved to [-L1, L1] x [-L2, L2]
#' (centred at the origin), the grid centroids are the HSGP input
#' coordinates and the HSGP domain is boundfactor * (L1, L2).
#'
#' @param dataset data.frame with x, y and expression column(s)
#' @param gene expression column (default: first non-coordinate column)
#' @param grid_size grid cells per axis (or c(nx, ny))
#' @param basis number of HSGP basis functions per axis (or c(m1, m2))
#' @param boundfactor HSGP boundary factor c (>= 1; 1.5 is a good default)
#' @param bbox optional ROI c(xmin, xmax, ymin, ymax); default = data range
#' @param L1,L2 optional half-ranges. Only used to warn if they disagree with
#'   the data (the old interface required them; they are now computed).
#' @return list(stan_data, grid, info)
denmark_prepare <- function(dataset, gene = NULL, grid_size = 10, basis = 25,
                            boundfactor = 1.5, bbox = NULL, L1 = NULL, L2 = NULL) {
  if (boundfactor < 1) stop("`boundfactor` must be >= 1.")
  g <- buildgridpp(grid_size, dataset, gene = gene, bbox = bbox)
  bb <- g$bbox
  centre <- c(x = mean(bb[c("xmin", "xmax")]), y = mean(bb[c("ymin", "ymax")]))
  half <- c(L1 = unname(diff(bb[c("xmin", "xmax")])) / 2,
            L2 = unname(diff(bb[c("ymin", "ymax")])) / 2)
  for (k in c("L1", "L2")) {
    user <- get(k)
    if (!is.null(user) && abs(user - half[[k]]) > 0.01 * half[[k]])
      warning(sprintf("%s = %.4g supplied but the data give %.4g; the data value is used.",
                      k, user, half[[k]]), call. = FALSE)
  }

  coords <- cbind(g$grid$x - centre["x"], g$grid$y - centre["y"])
  dimnames(coords) <- NULL
  S <- denmark_hsgp_indices(basis)
  dmax <- sqrt(diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2)

  stan_data <- list(
    n = nrow(coords),
    y1 = as.integer(g$grid$Y),               # outcome 1: cell counts
    y2 = as.integer(g$grid$M),               # outcome 2: summed expression
    d = 2L,
    mstar = nrow(S),                         # total number of basis functions
    coords = coords,                         # centred grid centroids
    log_grid_area = log(g$grid_area),        # log |g| (area offset)
    mrange = dmax / (2 * 2.75),              # scale of the inv-gamma prior on ell
    indices = S,
    L = as.numeric(boundfactor * half),      # HSGP domain half-widths
    is_centerted_PHI = 0L)

  grid <- g$grid
  grid$cx <- coords[, 1]; grid$cy <- coords[, 2]
  list(stan_data = stan_data,
       grid = grid,
       info = list(gene = g$gene, grid_size = g$grid_size, step = g$step,
                   grid_area = g$grid_area, bbox = bb, centre = centre,
                   L1 = half[["L1"]], L2 = half[["L2"]], basis = basis,
                   boundfactor = boundfactor, max_dist = dmax,
                   n_cells = sum(g$grid$Y), Ymat = g$Ymat, Mmat = g$Mmat))
}


## ------------------------------------------------------------------------- ##
## 2. Model compilation, initial values, fitting                             ##
## ------------------------------------------------------------------------- ##

.denmark_stan_path <- function(model, stan_file = NULL) {
  if (!is.null(stan_file)) {
    if (!file.exists(stan_file)) stop("Stan file not found: ", stan_file)
    return(normalizePath(stan_file))
  }
  fname <- switch(model, M2 = "DenMark-M2.stan", M1 = "DenMark-M1.stan")
  dirs <- unique(c(getOption("denmark.stan_dir"), .denmark_env$dir,
                   getwd(), file.path(getwd(), "DenMark")))
  for (d in dirs) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("Cannot find '", fname, "'. Pass `stan_file =` or set ",
       "options(denmark.stan_dir = '<folder with the .stan files>').")
}

# Recognise which model a Stan file implements from its parameter names.
.denmark_detect_model <- function(stan_file) {
  code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  if (grepl("cholesky_factor_corr", code)) "M2" else "M1"
}

#' Compile (once) and return a DenMark cmdstanr model
#'
#' @param model "M2" (full DenMark) or "M1" (simplified, independent fields)
#' @param stan_file optional path to a Stan file (overrides `model` lookup)
denmark_model <- function(model = c("M2", "M1"), stan_file = NULL,
                          force_recompile = FALSE, quiet = TRUE) {
  model <- match.arg(model)
  path <- .denmark_stan_path(model, stan_file)
  if (!force_recompile && !is.null(.denmark_env$models[[path]]))
    return(.denmark_env$models[[path]])
  mod <- cmdstanr::cmdstan_model(path, compile = TRUE, quiet = quiet,
                                 force_recompile = force_recompile)
  .denmark_env$models[[path]] <- mod
  mod
}


#' Initial values for DenMark
#'
#' Returns a function (as expected by cmdstanr) giving valid starting values.
#' For M2 the Cholesky factor of the correlation matrix must be lower
#' triangular with unit-norm rows: Astar = [[1, 0], [rho0, sqrt(1 - rho0^2)]].
#' The length-scales start at the prior mean (mrange) and the intercepts at
#' their empirical values, which shortens warm-up.
denmark_init <- function(model = c("M2", "M1"), stan_data, rho0 = 0.3) {
  model <- match.arg(model)
  if (abs(rho0) >= 1) stop("`rho0` must lie in (-1, 1).")
  mstar <- stan_data$mstar
  y1 <- stan_data$y1; y2 <- stan_data$y2
  beta0 <- log(max(mean(y1), 0.1)) - stan_data$log_grid_area
  beta1 <- log(max(sum(y2), 1) / sum(y1 + 0.05))
  ell <- rep(stan_data$mrange, 2)
  if (model == "M2") {
    function() list(Astar = matrix(c(1, rho0, 0, sqrt(1 - rho0^2)), 2, 2),
                    sigma = c(1, 0.5), beta0 = beta0, beta1 = beta1, ell = ell,
                    betab1 = rep(0, mstar), betab2 = rep(0, mstar))
  } else {
    function() list(a11 = 1, a22 = 0.5, beta0 = beta0, beta1 = beta1, ell = ell,
                    betab1 = rep(0, mstar), betab2 = rep(0, mstar))
  }
}


#' Fit DenMark to one gene
#'
#' @param dataset data.frame with columns x, y and expression column(s)
#'   (one row per cell; expression = raw counts)
#' @param gene name of the expression column (default: first non-x/y column)
#' @param model "M2" (full DenMark; estimates rho) or "M1" (simplified)
#' @param grid_size grid cells per axis, or c(nx, ny)
#' @param basis HSGP basis functions per axis, or c(m1, m2)
#' @param boundfactor HSGP boundary factor (>= 1)
#' @param bbox optional ROI c(xmin, xmax, ymin, ymax)
#' @param chains,parallel_chains,iter_warmup,iter_sampling,adapt_delta,
#'   max_treedepth,seed,refresh  passed to cmdstanr `$sample()`
#' @param init NULL (sensible defaults via denmark_init()) or anything
#'   accepted by cmdstanr
#' @param rho0 starting value of rho for M2
#' @param stan_file optional path to a custom Stan file
#' @param output_dir optional folder for the CmdStan csv files (created)
#' @param waic compute WAIC from the pointwise log-likelihood
#' @param keep_fit keep the CmdStanMCMC object (set FALSE in large scans)
#' @param quiet suppress sampler progress output
#' @param ... further arguments passed to `$sample()`
#' @return object of class "denmark" (a list), see print.denmark()
DenMark <- function(dataset, gene = NULL, model = c("M2", "M1"),
                    grid_size = 10, basis = 25, boundfactor = 1.5, bbox = NULL,
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, max_treedepth = 10,
                    seed = 2024, refresh = 400, init = NULL, rho0 = 0.3,
                    stan_file = NULL, output_dir = NULL, waic = TRUE,
                    keep_fit = TRUE, quiet = FALSE, L1 = NULL, L2 = NULL, ...) {
  model <- match.arg(model)
  if (!is.null(stan_file)) model <- .denmark_detect_model(stan_file)

  prep <- denmark_prepare(dataset, gene = gene, grid_size = grid_size,
                          basis = basis, boundfactor = boundfactor, bbox = bbox,
                          L1 = L1, L2 = L2)
  sd <- prep$stan_data
  if (sum(sd$y2) == 0) stop("Gene '", prep$info$gene, "' has zero counts in the ROI.")

  mod <- denmark_model(model, stan_file = stan_file)
  if (is.null(init)) init <- denmark_init(model, sd, rho0 = rho0)
  if (!is.null(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  run <- function() mod$sample(
    data = sd, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, refresh = if (quiet) 0 else refresh, init = init,
    output_dir = output_dir, show_messages = !quiet, ...)
  fit <- if (quiet) suppressMessages(run()) else run()

  out <- list(fit = fit, model = model, gene = prep$info$gene,
              stan_data = sd, grid = prep$grid, info = prep$info,
              settings = list(chains = chains, iter_warmup = iter_warmup,
                              iter_sampling = iter_sampling,
                              adapt_delta = adapt_delta,
                              max_treedepth = max_treedepth, seed = seed))
  class(out) <- "denmark"

  out$summary <- denmark_summary(out)
  out$diagnostics <- .denmark_diagnostics(out)
  out$time <- fit$time()$total / 60
  if (waic) {
    ll <- posterior::as_draws_matrix(fit$draws("log_lik"))
    out$waic <- suppressWarnings(loo::waic(matrix(ll, nrow = nrow(ll))))
  }
  out$hsgp <- denmark_check_hsgp(out, warn = !quiet)

  # fields kept for backward compatibility with DenMarkmain()
  out$rhat <- stats::setNames(out$summary$rhat, out$summary$variable)
  out$Ymat <- prep$info$Ymat
  out$Mmat <- prep$info$Mmat

  if (!keep_fit) {                      # keep only what post-processing needs
    out$draws_cache <- .denmark_draw_matrix(out)
    out$fit <- NULL
  }
  out
}


#' Backward-compatible wrapper of the original DenMark interface
#'
#' `denmark_version` may be a path to 'DenMark-M2.stan' / 'DenMark-M1.stan'
#' (as in the earlier tutorial) or the strings "M2" / "M1". L1 and L2 are no
#' longer needed (they are computed from the data) and are only checked.
DenMarkmain <- function(grid_size, dataset, basis, boundfactor,
                        L1 = NULL, L2 = NULL, denmark_version = "M2", ...) {
  if (denmark_version %in% c("M1", "M2")) {
    DenMark(dataset = dataset, model = denmark_version, grid_size = grid_size,
            basis = basis, boundfactor = boundfactor, L1 = L1, L2 = L2, ...)
  } else {
    DenMark(dataset = dataset, stan_file = denmark_version, grid_size = grid_size,
            basis = basis, boundfactor = boundfactor, L1 = L1, L2 = L2, ...)
  }
}


.denmark_key_vars <- function(model) {
  if (model == "M2")
    c("beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]", "sigma[1]", "sigma[2]")
  else
    c("beta0", "beta1", "ell[1]", "ell[2]", "a11", "a22")
}

.denmark_diagnostics <- function(obj) {
  fit <- obj$fit
  dg <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  s <- obj$summary
  list(num_divergent = if (is.null(dg)) NA_integer_ else sum(dg$num_divergent),
       num_max_treedepth = if (is.null(dg)) NA_integer_ else sum(dg$num_max_treedepth),
       min_ebfmi = if (is.null(dg)) NA_real_ else min(dg$ebfmi),
       max_rhat = max(s$rhat, na.rm = TRUE),
       min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
       min_ess_tail = min(s$ess_tail, na.rm = TRUE))
}


## ------------------------------------------------------------------------- ##
## 3. Post-processing                                                        ##
## ------------------------------------------------------------------------- ##

#' Posterior summary of the key DenMark parameters
#'
#' M2: beta0, beta1, ell[1], ell[2], Astar[2,1] (= rho), sigma[1], sigma[2]
#' M1: beta0, beta1, ell[1], ell[2], a11, a22
#' @param variables optional override of the variables to summarise
#' @return data.frame with mean, sd, quantiles, rhat, ess_bulk, ess_tail
denmark_summary <- function(obj, variables = NULL,
                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  stopifnot(inherits(obj, "denmark"))
  if (is.null(variables)) variables <- .denmark_key_vars(obj$model)
  if (!is.null(obj$fit)) {
    dr <- obj$fit$draws(variables = variables)
  } else {
    if (!all(variables %in% colnames(obj$draws_cache)))
      stop("The fit was not kept (keep_fit = FALSE); only key variables are available.")
    dr <- obj$draws_cache[, variables, drop = FALSE]
  }
  qfun <- function(x) stats::setNames(stats::quantile(x, probs, names = FALSE),
                                      paste0("q", probs * 100))
  s <- posterior::summarise_draws(dr, mean = mean, sd = stats::sd, qfun,
                                  rhat = posterior::rhat,
                                  ess_bulk = posterior::ess_bulk,
                                  ess_tail = posterior::ess_tail)
  s <- as.data.frame(s)
  s$label <- c(beta0 = "intercept, cell density", beta1 = "intercept, expression",
               `ell[1]` = "length-scale, density field",
               `ell[2]` = "length-scale, expression field",
               `Astar[2,1]` = "rho: density-expression correlation",
               `sigma[1]` = "SD, density field", `sigma[2]` = "SD, expression field",
               a11 = "SD, density field", a22 = "SD, expression field")[s$variable]
  s
}

# draws_matrix of the key scalar parameters (used when keep_fit = FALSE)
.denmark_draw_matrix <- function(obj) {
  posterior::as_draws_matrix(obj$fit$draws(variables = .denmark_key_vars(obj$model)))
}

#' Harmonised posterior draws
#'
#' One row per posterior draw with columns beta0, beta1, ell1, ell2, rho,
#' sigma1, sigma2 and the coregionalisation coefficients
#'   a11 = sigma1, a21 = sigma2 * rho, a22 = sigma2 * sqrt(1 - rho^2).
#' For M1, rho = a21 = 0 and sigma1 = a11, sigma2 = a22.
denmark_draws <- function(obj) {
  stopifnot(inherits(obj, "denmark"))
  dm <- if (!is.null(obj$fit)) .denmark_draw_matrix(obj) else obj$draws_cache
  dm <- unclass(dm)
  g <- function(v) as.numeric(dm[, v])
  if (obj$model == "M2") {
    rho <- g("Astar[2,1]"); s1 <- g("sigma[1]"); s2 <- g("sigma[2]")
    data.frame(beta0 = g("beta0"), beta1 = g("beta1"),
               ell1 = g("ell[1]"), ell2 = g("ell[2]"),
               rho = rho, sigma1 = s1, sigma2 = s2,
               a11 = s1, a21 = s2 * rho, a22 = s2 * sqrt(pmax(0, 1 - rho^2)))
  } else {
    a11 <- g("a11"); a22 <- g("a22")
    data.frame(beta0 = g("beta0"), beta1 = g("beta1"),
               ell1 = g("ell[1]"), ell2 = g("ell[2]"),
               rho = 0, sigma1 = a11, sigma2 = a22,
               a11 = a11, a21 = 0, a22 = a22)
  }
}


#' Latent-field decomposition for every grid cell
#'
#' Posterior means (and optional quantiles) of
#'   eta1    = beta0 + a11 w1             log cell intensity per unit area
#'   a11w1   = a11 w1                     density latent field
#'   eta2    = beta1 + a21 w1 + a22 w2    log expression per cell
#'   a21w1   = a21 w1                     density-shared part of expression
#'   a22w2   = a22 w2                     expression-specific part
#' The full linear predictors are log_lambda1 = log|g| + eta1 and
#' log_lambda2 = log(Y + 0.05) + eta2.
#' @return data.frame: one row per grid cell (grid_index, i, j, x, y, Y, M, ...)
denmark_fields <- function(obj, probs = NULL) {
  stopifnot(inherits(obj, "denmark"))
  if (is.null(obj$fit)) stop("denmark_fields() needs the fit (use keep_fit = TRUE).")
  n <- obj$stan_data$n
  W1 <- unclass(posterior::as_draws_matrix(obj$fit$draws("w1")))
  W2 <- unclass(posterior::as_draws_matrix(obj$fit$draws("w2")))
  W1 <- W1[, paste0("w1[", seq_len(n), "]"), drop = FALSE]
  W2 <- W2[, paste0("w2[", seq_len(n), "]"), drop = FALSE]
  d <- denmark_draws(obj)
  comp <- list(a11w1 = d$a11 * W1,               # row k scaled by draw k
               a21w1 = d$a21 * W1,
               a22w2 = d$a22 * W2)
  comp$eta1 <- d$beta0 + comp$a11w1
  comp$eta2 <- d$beta1 + comp$a21w1 + comp$a22w2
  out <- obj$grid[, c("grid_index", "i", "j", "x", "y", "Y", "M")]
  for (nm in c("eta1", "a11w1", "eta2", "a21w1", "a22w2")) {
    out[[nm]] <- colMeans(comp[[nm]])
    if (!is.null(probs)) for (p in probs)
      out[[sprintf("%s_q%s", nm, p * 100)]] <- apply(comp[[nm]], 2, stats::quantile, p)
  }
  out$log_lambda1 <- obj$stan_data$log_grid_area + out$eta1
  out$log_lambda2 <- log(out$Y + 0.05) + out$eta2
  out
}


#' Spatial correlation decay of the two latent processes
#'
#' Density process:    r1(d) = M(d; ell1)
#' Expression process: r2(d) = [a21^2 M(d; ell1) + a22^2 M(d; ell2)] / (a21^2 + a22^2)
#' with M the Matern-3/2 correlation. The practical range is the distance at
#' which the correlation falls to `threshold` (default 0.05), computed per
#' posterior draw.
#' @param distances distances at which to evaluate the curves (data units);
#'   default 0 .. max distance between grid centroids
#' @return list(curves = data.frame, ranges = data.frame)
denmark_decay <- function(obj, distances = NULL, probs = c(0.025, 0.975),
                          threshold = 0.05) {
  stopifnot(inherits(obj, "denmark"))
  d <- denmark_draws(obj)
  dmax <- obj$info$max_dist
  if (is.null(distances)) distances <- seq(0, dmax, length.out = 200)
  m1 <- outer(d$ell1, distances, function(l, r) matern32(r, 1, l))   # draws x dist
  m2 <- outer(d$ell2, distances, function(l, r) matern32(r, 1, l))
  w21 <- d$a21^2; w22 <- d$a22^2
  r2 <- (w21 * m1 + w22 * m2) / (w21 + w22)
  summ <- function(mat, field) data.frame(
    field = field, distance = distances, mean = colMeans(mat),
    lower = apply(mat, 2, stats::quantile, probs[1]),
    upper = apply(mat, 2, stats::quantile, probs[2]))
  curves <- rbind(summ(m1, "density"), summ(r2, "expression"))

  k <- .matern32_range_factor(threshold)
  range1 <- k * d$ell1
  range2 <- vapply(seq_len(nrow(d)), function(s) {
    if (w21[s] == 0) return(k * d$ell2[s])
    f <- function(r) (w21[s] * matern32(r, 1, d$ell1[s]) +
                      w22[s] * matern32(r, 1, d$ell2[s])) / (w21[s] + w22[s]) - threshold
    stats::uniroot(f, c(0, k * max(d$ell1[s], d$ell2[s]) * 1.01), tol = 1e-8)$root
  }, numeric(1))
  rsum <- function(v, field) data.frame(
    field = field, mean = mean(v), median = stats::median(v),
    lower = stats::quantile(v, probs[1], names = FALSE),
    upper = stats::quantile(v, probs[2], names = FALSE))
  list(curves = curves,
       ranges = rbind(rsum(range1, "density"), rsum(range2, "expression")),
       range_draws = data.frame(density = range1, expression = range2),
       threshold = threshold)
}


#' WAIC comparison of DenMark fits
#'
#' @param ... named "denmark" objects, e.g. denmark_compare(M2 = f2, M1 = f1)
#' @return data.frame ordered by WAIC (lower is better), with the difference
#'   in elpd relative to the best model and its standard error
denmark_compare <- function(...) {
  fits <- list(...)
  if (is.null(names(fits)) || any(names(fits) == ""))
    names(fits) <- vapply(fits, function(f) paste(f$model, f$gene), character(1))
  w <- lapply(fits, function(f) {
    if (is.null(f$waic)) stop("Fit without WAIC; refit with waic = TRUE.")
    f$waic
  })
  tab <- data.frame(model = names(w),
                    waic = vapply(w, function(x) x$estimates["waic", "Estimate"], 0),
                    waic_se = vapply(w, function(x) x$estimates["waic", "SE"], 0),
                    elpd_waic = vapply(w, function(x) x$estimates["elpd_waic", "Estimate"], 0),
                    p_waic = vapply(w, function(x) x$estimates["p_waic", "Estimate"], 0),
                    row.names = NULL)
  if (length(w) > 1) {
    cmp <- loo::loo_compare(w)
    tab$elpd_diff <- cmp[tab$model, "elpd_diff"]
    tab$se_diff <- cmp[tab$model, "se_diff"]
  }
  tab[order(tab$waic), , drop = FALSE]
}


#' Rule-of-thumb check of the HSGP approximation
#'
#' Uses the 1-D Matern-3/2 recommendations of Riutort-Mayol et al. (2023,
#' Statistics and Computing): with S the half-range of the domain, c the
#' boundary factor and m the number of basis functions, the smallest
#' length-scale that is well approximated is ell_min = 3.42 * c * S / m, and
#' the boundary factor should satisfy c >= max(1.2, 4.5 * ell / S).
#' Applied per axis; a warning is issued when the posterior median
#' length-scale is below ell_min (increase `basis`) or c is too small.
#' The rule targets an accurate approximation of the latent function itself
#' and is conservative for estimating rho; treat it as guidance.
denmark_check_hsgp <- function(obj, warn = TRUE) {
  d <- denmark_draws(obj)
  m <- obj$info$basis; if (length(m) == 1) m <- rep(m, 2)
  S <- c(obj$info$L1, obj$info$L2); cfac <- obj$info$boundfactor
  ell <- c(density = stats::median(d$ell1), expression = stats::median(d$ell2))
  tab <- expand.grid(field = names(ell), axis = c("x", "y"), stringsAsFactors = FALSE)
  tab$ell_median <- ell[tab$field]
  ax <- match(tab$axis, c("x", "y"))
  tab$ell_min <- 3.42 * cfac * S[ax] / m[ax]
  tab$c_needed <- pmax(1.2, 4.5 * tab$ell_median / S[ax])
  tab$basis_ok <- tab$ell_median >= tab$ell_min
  tab$boundary_ok <- cfac >= tab$c_needed
  if (warn && (!all(tab$basis_ok) || !all(tab$boundary_ok))) {
    msg <- character(0)
    if (!all(tab$basis_ok))
      msg <- c(msg, sprintf("posterior length-scale below %.3g: consider increasing `basis`",
                            max(tab$ell_min[!tab$basis_ok])))
    if (!all(tab$boundary_ok))
      msg <- c(msg, sprintf("consider boundfactor >= %.2f",
                            max(tab$c_needed[!tab$boundary_ok])))
    message("DenMark HSGP check (", obj$gene, ", rule of thumb): ",
            paste(msg, collapse = "; "))
  }
  tab
}


#' print method
print.denmark <- function(x, digits = 3, ...) {
  cat(sprintf("DenMark fit  |  model %s  |  gene %s\n", x$model, x$gene))
  cat(sprintf("  grid %d x %d (%d grids, %d cells)  |  basis %s  |  boundfactor %.2f\n",
              x$info$grid_size[1], x$info$grid_size[2], x$stan_data$n,
              x$info$n_cells, paste(x$info$basis, collapse = "x"), x$info$boundfactor))
  s <- x$summary
  if (x$model == "M2") {
    r <- s[s$variable == "Astar[2,1]", ]
    cat(sprintf("  rho (density-expression correlation): %.3f  [95%% CrI %.3f, %.3f]\n",
                r$mean, r$q2.5, r$q97.5))
  }
  if (!is.null(x$waic))
    cat(sprintf("  WAIC: %.2f (SE %.2f)\n", x$waic$estimates["waic", "Estimate"],
                x$waic$estimates["waic", "SE"]))
  dg <- x$diagnostics
  cat(sprintf("  max R-hat %.3f | min ESS bulk %.0f | divergences %s | time %.2f min\n",
              dg$max_rhat, dg$min_ess_bulk, dg$num_divergent, x$time))
  cat("\n")
  print(s[, c("variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat", "ess_bulk")],
        digits = digits, row.names = FALSE)
  invisible(x)
}


#' Save / load a DenMark fit
#'
#' cmdstanr keeps posterior draws in CSV files (by default in a temporary
#' folder that disappears when R closes). denmark_save() reads everything
#' into memory first, so the saved object can be reloaded and post-processed
#' in a later session.
denmark_save <- function(obj, file) {
  stopifnot(inherits(obj, "denmark"))
  if (!is.null(obj$fit)) {
    obj$fit$draws()
    try(obj$fit$sampler_diagnostics(), silent = TRUE)
    try(obj$fit$init(), silent = TRUE)
    try(obj$fit$profiles(), silent = TRUE)
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, file)
  invisible(file)
}

denmark_load <- function(file) {
  obj <- readRDS(file)
  stopifnot(inherits(obj, "denmark"))
  obj
}


## ------------------------------------------------------------------------- ##
## 4. DCG discovery: scan many genes                                         ##
## ------------------------------------------------------------------------- ##

#' Fit DenMark gene by gene and call density-correlated genes (DCG)
#'
#' A gene is a positively density-correlated gene (PCG) when the
#' `cred_level` credible interval of rho lies above 0, a negatively
#' density-correlated gene (NCG) when it lies below 0, and "not DCG"
#' otherwise.
#'
#' @param data data.frame: columns x, y and one column per gene (cells in rows)
#' @param genes genes to scan (default: all non-coordinate columns)
#' @param fit_simple also fit M1 and report the WAIC difference (doubles time)
#' @param save_dir optional folder; one .rds per gene is written and genes
#'   already present are loaded instead of refitted (resumable scans)
#' @param cred_level credible level used to call DCGs
#' @param ... further arguments passed to DenMark() (grid_size, basis, ...)
#' @return data.frame, one row per gene
denmark_scan <- function(data, genes = NULL, fit_simple = FALSE, save_dir = NULL,
                         cred_level = 0.95, quiet = TRUE, verbose = TRUE, ...) {
  data <- as.data.frame(data)
  if (is.null(genes)) genes <- setdiff(names(data), c("x", "y"))
  genes <- intersect(genes, names(data))
  if (!length(genes)) stop("None of the requested genes are columns of `data`.")
  if (!is.null(save_dir)) dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  a <- (1 - cred_level) / 2

  one_gene <- function(g) {
    f_rds <- if (!is.null(save_dir)) file.path(save_dir, paste0("denmark_", make.names(g), ".rds"))
    if (!is.null(f_rds) && file.exists(f_rds)) return(readRDS(f_rds))
    t0 <- Sys.time()
    row <- tryCatch({
      f2 <- DenMark(data[, c("x", "y", g)], gene = g, model = "M2",
                    keep_fit = FALSE, quiet = quiet, ...)
      rho <- denmark_draws(f2)$rho
      r <- data.frame(gene = g,
                      rho_mean = mean(rho), rho_median = stats::median(rho),
                      rho_lower = stats::quantile(rho, a, names = FALSE),
                      rho_upper = stats::quantile(rho, 1 - a, names = FALSE),
                      prob_positive = mean(rho > 0),
                      total_counts = sum(f2$stan_data$y2),
                      max_rhat = f2$diagnostics$max_rhat,
                      min_ess_bulk = f2$diagnostics$min_ess_bulk,
                      divergences = f2$diagnostics$num_divergent,
                      waic_M2 = f2$waic$estimates["waic", "Estimate"],
                      time_min = f2$time, error = NA_character_)
      if (fit_simple) {
        f1 <- DenMark(data[, c("x", "y", g)], gene = g, model = "M1",
                      keep_fit = FALSE, quiet = quiet, ...)
        r$waic_M1 <- f1$waic$estimates["waic", "Estimate"]
        r$delta_waic <- r$waic_M1 - r$waic_M2     # > 0 favours the full model
        r$time_min <- r$time_min + f1$time
      }
      r
    }, error = function(e) data.frame(gene = g, error = conditionMessage(e)))
    if (!is.null(f_rds) && is.na(row$error[1])) saveRDS(row, f_rds)
    if (verbose) {
      if (is.na(row$error[1])) {
        message(sprintf("[%s] rho = %+.3f [%+.3f, %+.3f]  R-hat %.3f  (%.1f min)",
                        g, row$rho_mean, row$rho_lower, row$rho_upper,
                        row$max_rhat, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      } else message(sprintf("[%s] failed: %s", g, row$error))
    }
    row
  }

  rows <- lapply(genes, one_gene)
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA; r[cols]
  }))
  res$class <- ifelse(is.na(res$rho_lower), NA_character_,
               ifelse(res$rho_lower > 0, "PCG",
               ifelse(res$rho_upper < 0, "NCG", "not DCG")))
  res$converged <- !is.na(res$max_rhat) & res$max_rhat < 1.01 &
                   !is.na(res$divergences) & res$divergences == 0
  res <- res[order(-res$rho_mean), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "cred_level") <- cred_level
  res
}


## ------------------------------------------------------------------------- ##
## 5. Visualisation                                                          ##
## ------------------------------------------------------------------------- ##

.theme_denmark <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 1),
                   legend.key.height = ggplot2::unit(0.9, "lines"))
}

.tile_plot <- function(df, value, title, limits = NULL, legend = NULL,
                       xlab = "x", ylab = "y", diverging = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                        fill = .data[[value]])) +
    ggplot2::geom_tile(colour = "grey30", linewidth = 0.15) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = title, x = xlab, y = ylab, fill = legend) +
    .theme_denmark()
  if (diverging) {
    lim <- if (is.null(limits)) max(abs(df[[value]])) * c(-1, 1) else limits
    p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                      limits = lim)
  } else {
    p + ggplot2::scale_fill_viridis_c(option = "magma", limits = limits)
  }
}

.combine <- function(plots, ncol) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = ncol)
  } else {
    message("Install 'patchwork' to combine the panels; returning a list of plots.")
    plots
  }
}

#' Gridded cell counts (Y) and summed expression (M)
plot_denmark_grid <- function(obj, xlab = "x", ylab = "y") {
  g <- if (inherits(obj, "denmark")) obj$grid else obj$grid
  gene <- if (inherits(obj, "denmark")) obj$gene else obj$gene
  .combine(list(
    .tile_plot(g, "Y", "Gridded cell counts", legend = "cells", xlab = xlab, ylab = ylab),
    .tile_plot(g, "M", paste0("Gridded expression: ", gene), legend = "counts",
               xlab = xlab, ylab = ylab)), ncol = 2)
}

#' Latent-field decomposition maps (posterior means)
#'
#' Row 1 (cell density): beta0 + a11 w1 and a11 w1 (shared colour scale)
#' Row 2 (expression): beta1 + a21 w1 + a22 w2, a21 w1 and a22 w2 (shared scale)
plot_denmark_fields <- function(obj, fields = NULL, xlab = "x", ylab = "y") {
  if (is.null(fields)) fields <- denmark_fields(obj)
  lim1 <- range(c(fields$eta1, fields$a11w1))
  lim2 <- range(c(fields$eta2, fields$a21w1, fields$a22w2))
  ttl <- list(eta1  = quote(beta[0] + a[11] * w[1]),
              a11w1 = quote(a[11] * w[1]),
              eta2  = quote(beta[1] + a[21] * w[1] + a[22] * w[2]),
              a21w1 = quote(a[21] * w[1]),
              a22w2 = quote(a[22] * w[2]))
  leg <- c(eta1 = "density", a11w1 = "density", eta2 = "expression",
           a21w1 = "expression", a22w2 = "expression")
  p <- lapply(names(ttl), function(v)
    .tile_plot(fields, v, as.expression(ttl[[v]]),
               limits = if (v %in% c("eta1", "a11w1")) lim1 else lim2,
               legend = leg[[v]], xlab = xlab, ylab = ylab))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    row1 <- (p[[1]] | p[[2]] | patchwork::plot_spacer()) +
      patchwork::plot_layout(guides = "collect")
    row2 <- (p[[3]] | p[[4]] | p[[5]]) + patchwork::plot_layout(guides = "collect")
    row1 / row2
  } else p
}

#' Spatial correlation decay curves with posterior bands
plot_denmark_decay <- function(obj, decay = NULL, xlab = "Distance") {
  if (is.null(decay)) decay <- denmark_decay(obj)
  cv <- decay$curves
  cv$field <- factor(cv$field, c("density", "expression"),
                     c("r1: cell density field", "r2: expression field"))
  ggplot2::ggplot(cv, ggplot2::aes(x = .data$distance)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), linewidth = 1) +
    ggplot2::geom_hline(yintercept = decay$threshold, linetype = 2, colour = "grey40") +
    ggplot2::facet_wrap(~ field) +
    ggplot2::labs(x = xlab, y = "Posterior spatial correlation") +
    .theme_denmark()
}

#' Forest plot of rho across genes from denmark_scan()
plot_denmark_scan <- function(scan, top = NULL) {
  s <- scan[!is.na(scan$rho_mean), , drop = FALSE]
  if (!is.null(top) && nrow(s) > top)
    s <- s[order(-abs(s$rho_mean)), ][seq_len(top), ]
  s$gene <- factor(s$gene, levels = s$gene[order(s$rho_mean)])
  s$class <- factor(s$class, c("PCG", "NCG", "not DCG"))
  ggplot2::ggplot(s, ggplot2::aes(x = .data$rho_mean, y = .data$gene,
                                  colour = .data$class)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$rho_lower, xmax = .data$rho_upper),
                            height = 0.25) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(PCG = "#B2182B", NCG = "#2166AC",
                                            `not DCG` = "grey50"), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::labs(x = expression("Density-expression correlation " * rho *
                                 " (posterior mean, credible interval)"),
                  y = NULL, colour = NULL) +
    .theme_denmark()
}


## ------------------------------------------------------------------------- ##
## 6. Simulation (for testing and the tutorial)                              ##
## ------------------------------------------------------------------------- ##

# Smooth zero-mean, unit-variance Gaussian field with (approximately)
# Matern-3/2 covariance, via random Fourier features: the Matern-3/2
# spectral density in 2-D is a bivariate t with 3 degrees of freedom.
.rff_field <- function(ell, n_features = 500) {
  z <- matrix(stats::rnorm(2 * n_features), ncol = 2)
  u <- stats::rchisq(n_features, df = 3)
  omega <- (sqrt(3) / ell) * z / sqrt(u / 3)
  phase <- stats::runif(n_features, 0, 2 * pi)
  amp <- stats::rnorm(n_features)
  function(x, y) {
    arg <- outer(x, omega[, 1]) + outer(y, omega[, 2]) +
      matrix(phase, length(x), n_features, byrow = TRUE)
    as.numeric(cos(arg) %*% amp) * sqrt(2 / n_features)
  }
}

#' Simulate a marked point pattern from the DenMark data-generating process
#'
#' Cells follow an inhomogeneous Poisson process with
#'   log lambda1(s) = beta0 + a11 w1(s)
#' and each cell's count for gene k is Poisson with
#'   log mu_k(s) = beta1 + sigma2 * (rho_k w1(s) + sqrt(1 - rho_k^2) w2k(s)).
#'
#' @param rho named vector of true correlations, one gene per element
#' @param width,height size of the window [0, width] x [0, height]
#' @param n_cells expected number of cells
#' @return data.frame with x, y and one count column per gene; the true
#'   values are stored in attr(, "truth")
denmark_simulate <- function(rho = c(geneA = 0.6, geneB = 0, geneC = -0.6),
                             width = 6, height = 4, n_cells = 5000,
                             ell1 = 1, ell2 = 1, a11 = 1, sigma2 = 0.8,
                             beta1 = 0.5, seed = 1) {
  if (is.null(names(rho))) names(rho) <- paste0("gene", seq_along(rho))
  set.seed(seed)
  w1 <- .rff_field(ell1)
  # thinning: simulate a homogeneous process and keep with prob lambda/lambda_max
  xs <- stats::runif(60000, 0, width); ys <- stats::runif(60000, 0, height)
  f1 <- w1(xs, ys)
  lam <- exp(a11 * f1)
  keep <- stats::runif(length(lam)) < lam / max(lam)
  xs <- xs[keep]; ys <- ys[keep]; f1 <- f1[keep]
  if (length(xs) > n_cells) {
    s <- sample.int(length(xs), n_cells); xs <- xs[s]; ys <- ys[s]; f1 <- f1[s]
  }
  out <- data.frame(x = xs, y = ys)
  z1 <- (f1 - mean(f1)) / stats::sd(f1)
  for (g in names(rho)) {
    w2 <- .rff_field(ell2)(xs, ys)
    # orthogonalise w2 against w1 over the observed cells, so that the
    # realised correlation of the two fields equals rho (two independent
    # smooth fields can be correlated by chance on a finite window)
    w2 <- stats::residuals(stats::lm(w2 ~ z1))
    w2 <- w2 / stats::sd(w2)
    eta <- beta1 + sigma2 * (rho[[g]] * z1 + sqrt(1 - rho[[g]]^2) * w2)
    out[[g]] <- stats::rpois(length(xs), exp(eta))
  }
  attr(out, "truth") <- list(rho = rho, ell1 = ell1, ell2 = ell2, a11 = a11,
                             sigma2 = sigma2, beta1 = beta1,
                             window = c(width = width, height = height))
  out
}

###############################################################################
# End of DenMark-functions.R
###############################################################################
