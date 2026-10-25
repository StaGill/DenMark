## ------------------------------------------------------------------------- ##
## 3. Post-processing                                                        ##
## ------------------------------------------------------------------------- ##

#' Posterior summary of the key DenMark parameters
#'
#' The key parameters are `beta0`, `beta1`, `ell[1]`, `ell[2]`, `Astar[2,1]`
#' (= rho), `sigma[1]` and `sigma[2]`.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @param variables Optional override of the variables to summarize.
#' @param probs Probabilities of the posterior quantiles to report.
#' @return Data frame with the columns `variable`, `mean`, `sd`, one column
#'   per quantile (`q2.5`, `q25`, ...), `rhat`, `ess_bulk`, `ess_tail` and a
#'   descriptive `label`.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' denmark_summary(fit)
denmark_summary <- function(obj, variables = NULL,
                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  stopifnot(inherits(obj, "denmark"))
  if (is.null(variables)) variables <- .denmark_key_vars()
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
               `sigma[1]` = "SD, density field", `sigma[2]` = "SD, expression field")[s$variable]
  s
}

# draws_matrix of the key scalar parameters (used when keep_fit = FALSE)
.denmark_draw_matrix <- function(obj) {
  posterior::as_draws_matrix(obj$fit$draws(variables = .denmark_key_vars()))
}

#' Harmonized posterior draws
#'
#' One row per posterior draw with the columns `beta0`, `beta1`, `ell1`,
#' `ell2`, `rho`, `sigma1`, `sigma2` and the coregionalization coefficients
#' `a11 = sigma1`, `a21 = sigma2 * rho` and `a22 = sigma2 * sqrt(1 - rho^2)`.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @return Data frame with one row per posterior draw.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' draws <- denmark_draws(fit)
#' quantile(draws$rho, c(0.025, 0.5, 0.975))
denmark_draws <- function(obj) {
  stopifnot(inherits(obj, "denmark"))
  dm <- if (!is.null(obj$fit)) .denmark_draw_matrix(obj) else obj$draws_cache
  dm <- unclass(dm)
  g <- function(v) as.numeric(dm[, v])
  rho <- g("Astar[2,1]"); s1 <- g("sigma[1]"); s2 <- g("sigma[2]")
  data.frame(beta0 = g("beta0"), beta1 = g("beta1"),
             ell1 = g("ell[1]"), ell2 = g("ell[2]"),
             rho = rho, sigma1 = s1, sigma2 = s2,
             a11 = s1, a21 = s2 * rho, a22 = s2 * sqrt(pmax(0, 1 - rho^2)))
}


#' Latent-field decomposition for every grid cell
#'
#' Posterior means (and optional quantiles) of
#' * `eta1 = beta0 + a11 w1`, the log cell intensity per unit area;
#' * `a11w1 = a11 w1`, the density latent field;
#' * `eta2 = beta1 + a21 w1 + a22 w2`, the log expression per cell;
#' * `a21w1 = a21 w1`, the density-shared part of expression;
#' * `a22w2 = a22 w2`, the expression-specific part.
#'
#' The full linear predictors are `log_lambda1 = log|g| + eta1` and
#' `log_lambda2 = log(Y + 0.05) + eta2`.
#'
#' @param obj A `denmark` object fitted with `keep_fit = TRUE`.
#' @param probs Optional probabilities; for each, the posterior quantiles of
#'   the five fields are added in columns such as `eta1_q2.5`.
#' @return Data frame with one row per grid cell: `grid_index`, `i`, `j`,
#'   `x`, `y`, `Y`, `M`, the five fields, `log_lambda1` and `log_lambda2`.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' fields <- denmark_fields(fit)
#' head(fields)
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
#' Density process: `r1(d) = M(d; ell1)`. Expression process:
#' `r2(d) = [a21^2 M(d; ell1) + a22^2 M(d; ell2)] / (a21^2 + a22^2)`, with `M`
#' the Matern-3/2 correlation. The practical range is the distance at which
#' the correlation falls to `threshold` (default 0.05), computed per
#' posterior draw.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @param distances Distances at which to evaluate the curves (data units);
#'   default: 200 points from 0 to the maximum distance between grid
#'   centroids.
#' @param probs Lower and upper probabilities of the posterior band.
#' @param threshold Correlation level that defines the practical range.
#' @return A list with `curves` (data frame: field, distance, posterior mean,
#'   lower and upper band), `ranges` (data frame: posterior mean, median and
#'   interval of the practical range of each field), `range_draws` (the
#'   practical ranges per draw) and `threshold`.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' decay <- denmark_decay(fit)
#' decay$ranges
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


#' Rule-of-thumb check of the HSGP approximation
#'
#' Uses the one-dimensional Matern-3/2 recommendations of Riutort-Mayol et
#' al. (2023, Statistics and Computing): with `S` the half-range of the
#' domain, `c` the boundary factor and `m` the number of basis functions, the
#' smallest length-scale that is well approximated is
#' `ell_min = 3.42 * c * S / m`, and the boundary factor should satisfy
#' `c >= max(1.2, 4.5 * ell / S)`. Applied per axis; a message is issued when
#' the posterior median length-scale is below `ell_min` (increase `basis`)
#' or `c` is too small. The rule targets an accurate approximation of the
#' latent function itself and is conservative for estimating rho; treat it
#' as guidance.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @param warn Issue a message when the rule is not met.
#' @return Data frame with one row per field and axis: `field`, `axis`,
#'   `ell_median`, `ell_min`, `c_needed`, `basis_ok` and `boundary_ok`.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' denmark_check_hsgp(fit, warn = FALSE)
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


#' @rdname DenMark
#' @param x A `denmark` object.
#' @param digits Number of significant digits in the printed summary table.
#' @export
print.denmark <- function(x, digits = 3, ...) {
  cat(sprintf("DenMark fit  |  gene %s\n", x$gene))
  cat(sprintf("  grid %d x %d (%d grids, %d cells)  |  basis %s  |  boundfactor %.2f\n",
              x$info$grid_size[1], x$info$grid_size[2], x$stan_data$n,
              x$info$n_cells, paste(x$info$basis, collapse = "x"), x$info$boundfactor))
  s <- x$summary
  r <- s[s$variable == "Astar[2,1]", ]
  cat(sprintf("  rho (density-expression correlation): %.3f  [95%% CrI %.3f, %.3f]\n",
              r$mean, r$q2.5, r$q97.5))
  dg <- x$diagnostics
  cat(sprintf("  max R-hat %.3f | min ESS bulk %.0f | divergences %s | time %.2f min\n",
              dg$max_rhat, dg$min_ess_bulk, dg$num_divergent, x$time))
  cat("\n")
  print(s[, c("variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat", "ess_bulk")],
        digits = digits, row.names = FALSE)
  invisible(x)
}


#' Save and load a DenMark fit
#'
#' cmdstanr keeps posterior draws in CSV files (by default in a temporary
#' folder that disappears when R closes). `denmark_save()` reads everything
#' into memory first, so the saved object can be reloaded and post-processed
#' in a later session.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @param file Path of the `.rds` file.
#' @return `denmark_save()` returns `file` invisibly; `denmark_load()` returns
#'   the `denmark` object.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' f <- tempfile(fileext = ".rds")
#' denmark_save(fit, f)
#' fit2 <- denmark_load(f)
#' denmark_summary(fit2)
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

#' @rdname denmark_save
#' @export
denmark_load <- function(file) {
  obj <- readRDS(file)
  stopifnot(inherits(obj, "denmark"))
  obj
}
