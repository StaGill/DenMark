## ------------------------------------------------------------------------- ##
## 2. Model compilation, initial values, fitting                             ##
## ------------------------------------------------------------------------- ##

.denmark_stan_path <- function(stan_file = NULL, vararea = FALSE) {
  if (!is.null(stan_file)) {
    if (!file.exists(stan_file)) stop("Stan file not found: ", stan_file)
    return(normalizePath(stan_file))
  }
  fname <- paste0("DenMark", if (vararea) "-vararea", ".stan")
  dirs <- unique(c(getOption("denmark.stan_dir"),
                   system.file("stan", package = "DenMark")))
  for (d in dirs) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("Cannot find '", fname, "'. Pass `stan_file =` or set ",
       "options(denmark.stan_dir = '<folder with the .stan files>').")
}

#' Compile (once) and return the DenMark cmdstanr model
#'
#' `DenMark.stan`, or `DenMark-vararea.stan` with `vararea = TRUE`, is looked
#' up in `getOption("denmark.stan_dir")` and then in the `stan` folder of the
#' installed package. The compiled model is cached for the R session. The
#' executables of the packaged models are written to
#' `tools::R_user_dir("DenMark", which = "cache")` (one subfolder per package
#' version, because the library may be read-only), so that later sessions
#' reuse them; any other Stan file is compiled next to itself.
#'
#' @param stan_file Optional path to a Stan file, compiled instead of the
#'   packaged model.
#' @param vararea Use the variable-area model `DenMark-vararea.stan`, the same
#'   model with one area offset `log_grid_area` per grid cell (used with a
#'   tissue boundary).
#' @param force_recompile Recompile even if a compiled model is available.
#' @param quiet Suppress the compiler output.
#' @return A cmdstanr `CmdStanModel` object.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
#' mod <- denmark_model()
#' mod$exe_file()
denmark_model <- function(stan_file = NULL, vararea = FALSE, force_recompile = FALSE,
                          quiet = TRUE) {
  path <- .denmark_stan_path(stan_file, vararea)
  if (!force_recompile && !is.null(.denmark_env$models[[path]]))
    return(.denmark_env$models[[path]])
  cache <- NULL
  pkg_stan <- system.file("stan", package = "DenMark")
  if (nzchar(pkg_stan) && normalizePath(dirname(path), winslash = "/") ==
      normalizePath(pkg_stan, winslash = "/")) {
    cache <- file.path(tools::R_user_dir("DenMark", which = "cache"),
                       as.character(utils::packageVersion("DenMark")))
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  }
  mod <- cmdstanr::cmdstan_model(path, compile = TRUE, quiet = quiet,
                                 force_recompile = force_recompile, dir = cache)
  .denmark_env$models[[path]] <- mod
  mod
}


#' Initial values for DenMark
#'
#' Returns a function (as expected by cmdstanr) giving valid starting values.
#' The Cholesky factor of the correlation matrix must be lower triangular
#' with unit-norm rows: `Astar = [[1, 0], [rho0, sqrt(1 - rho0^2)]]`.
#' The length-scales start at the prior mean (`mrange`) and the intercepts at
#' their empirical values (for `beta0`, the number of cells per unit of
#' tissue area), which shortens warm-up.
#'
#' @param stan_data Stan data list, as in `denmark_prepare()$stan_data`.
#' @param rho0 Starting value of rho, in (-1, 1).
#' @return A function without arguments that returns a list of initial
#'   values.
#' @export
#' @examples
#' toy <- denmark_toy_data()
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' prep <- denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10)
#' init <- denmark_init(prep$stan_data)
#' str(init())
denmark_init <- function(stan_data, rho0 = 0.3) {
  if (abs(rho0) >= 1) stop("`rho0` must lie in (-1, 1).")
  mstar <- stan_data$mstar
  y1 <- stan_data$y1; y2 <- stan_data$y2
  area <- exp(rep_len(stan_data$log_grid_area, length(y1)))
  beta0 <- log(max(sum(y1), 1) / sum(area))
  beta1 <- log(max(sum(y2), 1) / sum(y1 + 0.05))
  ell <- rep(stan_data$mrange, 2)
  function() list(Astar = matrix(c(1, rho0, 0, sqrt(1 - rho0^2)), 2, 2),
                  sigma = c(1, 0.5), beta0 = beta0, beta1 = beta1, ell = ell,
                  betab1 = rep(0, mstar), betab2 = rep(0, mstar))
}


#' Fit DenMark to one gene
#'
#' Grids the data, builds the HSGP inputs, compiles the Stan model (once per
#' session), sets valid initial values and samples with cmdstanr.
#'
#' The DenMark model (`DenMark.stan`) has two dependent latent fields,
#' `log lambda1(g) = log|g| + beta0 + a11 w1(g)` and
#' `log lambda2(g) = beta1 + log(Y(g) + 0.05) + a21 w1(g) + a22 w2(g)`;
#' `rho = Astar[2,1]` is the density-expression correlation.
#'
#' A tissue section is rarely a rectangle. With `boundary` (an alpha shape of
#' all cells of the section, from [denmark_alpha_shape()]), the grid cells
#' outside the tissue are dropped and every remaining grid cell enters with
#' its area inside the tissue, through the variable-area model
#' `DenMark-vararea.stan` (the same model with one `log_grid_area` per grid
#' cell). Without `boundary`, the whole bounding box is gridded with a common
#' grid area.
#'
#' The print method reports rho, the sampler diagnostics and the posterior
#' summary.
#'
#' @param dataset Data frame with columns `x`, `y` and the expression
#'   column(s), one row per cell; expression values are raw counts.
#' @param gene Name of the expression column (default: the first column other
#'   than `x` and `y`).
#' @param grid_size Grid cells per axis, or `c(nx, ny)`.
#' @param basis HSGP basis functions per axis, or `c(m1, m2)`.
#' @param boundfactor HSGP boundary factor (>= 1).
#' @param bbox Optional region of interest `c(xmin, xmax, ymin, ymax)`.
#' @param chains,parallel_chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth,seed,refresh
#'   Passed to the cmdstanr `$sample()` method.
#' @param init `NULL` (valid defaults from [denmark_init()]) or anything
#'   accepted by cmdstanr.
#' @param rho0 Starting value of rho.
#' @param stan_file Optional path to a custom Stan file with the parameters
#'   of `DenMark.stan`; with a `boundary` it must declare
#'   `vector[n] log_grid_area`.
#' @param output_dir Optional folder for the CmdStan CSV files (created if
#'   needed).
#' @param keep_fit Keep the `CmdStanMCMC` object (set `FALSE` in large
#'   scans).
#' @param quiet Suppress the sampler progress output.
#' @param L1,L2 Optional half-ranges of the region of interest; only checked
#'   against the data (see [denmark_prepare()]).
#' @param boundary Optional tissue polygon from [denmark_alpha_shape()]; grid
#'   cells outside the tissue are dropped and boundary grid cells use their
#'   clipped area (variable-area model).
#' @param min_area_frac See [denmark_prepare()].
#' @param ... For `DenMark()`, further arguments passed to `$sample()`;
#'   ignored by `print()`.
#' @return An object of class `denmark`: a list with the cmdstanr fit
#'   (`fit`), `gene`, `stan_data`, `grid`, `info`, the sampler `settings`,
#'   the posterior `summary` (see [denmark_summary()]), the sampler
#'   `diagnostics`, the run `time` in minutes, the HSGP check `hsgp` (see
#'   [denmark_check_hsgp()]) and, for backward compatibility with
#'   [DenMarkmain()], `rhat`, `Ymat` and `Mmat`. With `keep_fit = FALSE` the
#'   fit is dropped and the draws of the key parameters are kept in
#'   `draws_cache`.
#' @seealso [denmark_scan()] to fit many genes.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10, chains = 2,
#'                iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' fit
DenMark <- function(dataset, gene = NULL,
                    grid_size = 10, basis = 25, boundfactor = 1.5, bbox = NULL,
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, max_treedepth = 10,
                    seed = 2024, refresh = 400, init = NULL, rho0 = 0.3,
                    stan_file = NULL, output_dir = NULL,
                    keep_fit = TRUE, quiet = FALSE, L1 = NULL, L2 = NULL,
                    boundary = NULL, min_area_frac = 0, ...) {
  .denmark_removed_args(...)
  prep <- denmark_prepare(dataset, gene = gene, grid_size = grid_size,
                          basis = basis, boundfactor = boundfactor, bbox = bbox,
                          L1 = L1, L2 = L2, boundary = boundary,
                          min_area_frac = min_area_frac)
  sd <- prep$stan_data
  if (sum(sd$y2) == 0) stop("Gene '", prep$info$gene, "' has zero counts in the ROI.")
  if (!is.null(stan_file) && prep$info$vararea &&
      !any(grepl("vector\\[n\\][[:space:]]*log_grid_area",
                 readLines(stan_file, warn = FALSE))))
    stop("With a tissue boundary the Stan model needs `vector[n] log_grid_area` ",
         "(use the '-vararea' models).")

  mod <- denmark_model(stan_file = stan_file, vararea = prep$info$vararea)
  if (is.null(init)) init <- denmark_init(sd, rho0 = rho0)
  if (!is.null(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  run <- function() mod$sample(
    data = sd, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, refresh = if (quiet) 0 else refresh, init = init,
    output_dir = output_dir, show_messages = !quiet, ...)
  fit <- if (quiet) suppressMessages(run()) else run()

  out <- list(fit = fit, gene = prep$info$gene,
              stan_data = sd, grid = prep$grid, info = prep$info,
              settings = list(chains = chains, iter_warmup = iter_warmup,
                              iter_sampling = iter_sampling,
                              adapt_delta = adapt_delta,
                              max_treedepth = max_treedepth, seed = seed))
  class(out) <- "denmark"

  out$summary <- denmark_summary(out)
  out$diagnostics <- .denmark_diagnostics(out)
  out$time <- fit$time()$total / 60
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
#' Keeps scripts written for the original interface running. `L1` and `L2`
#' are no longer needed (they are computed from the data) and are only
#' checked.
#'
#' @param grid_size Grid cells per axis, or `c(nx, ny)`.
#' @param dataset Data frame with columns `x`, `y` and one expression column.
#' @param basis HSGP basis functions per axis, or `c(m1, m2)`.
#' @param boundfactor HSGP boundary factor (>= 1).
#' @param L1,L2 Optional half-ranges, only checked against the data.
#' @param denmark_version `"M2"` (the default) fits the DenMark model; any
#'   other value is taken as the path to a Stan file (as in the earlier
#'   tutorial) and passed to [DenMark()] as `stan_file`.
#' @param ... Further arguments passed to [DenMark()].
#' @return An object of class `denmark`, see [DenMark()].
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- toy[toy$celltype == "Tumor", c("x", "y", "geneA")]
#' fit <- DenMarkmain(grid_size = 10, dataset = tumor,
#'                    basis = 10, boundfactor = 1.5, boundary = tissue,
#'                    chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                    seed = 2024, quiet = TRUE)
#' fit$rhat
DenMarkmain <- function(grid_size, dataset, basis, boundfactor,
                        L1 = NULL, L2 = NULL, denmark_version = "M2", ...) {
  if (identical(denmark_version, "M1"))
    stop("The simplified model (M1) was removed in DenMark 0.9.0; ",
         "DenMark() fits the main model.")
  if (identical(denmark_version, "M2")) {
    DenMark(dataset = dataset, grid_size = grid_size,
            basis = basis, boundfactor = boundfactor, L1 = L1, L2 = L2, ...)
  } else {
    DenMark(dataset = dataset, stan_file = denmark_version, grid_size = grid_size,
            basis = basis, boundfactor = boundfactor, L1 = L1, L2 = L2, ...)
  }
}


.denmark_key_vars <- function() {
  c("beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]", "sigma[1]", "sigma[2]")
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
