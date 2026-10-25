## ------------------------------------------------------------------------- ##
## 0b. Small utilities                                                       ##
## ------------------------------------------------------------------------- ##

#' Posterior quantile helpers
#'
#' Shorthands for the 2.5th, 25th, 50th, 75th and 97.5th percentiles of a
#' vector of posterior draws, kept for compatibility with earlier scripts.
#'
#' @param x Numeric vector, for example the posterior draws of one parameter.
#' @return A single unnamed number: the sample quantile of `x` from
#'   `stats::quantile()` with its default method.
#' @export
#' @examples
#' x <- 0:100
#' c(q2.5(x), q25(x), q50(x), q75(x), q97.5(x))
q2.5  <- function(x) stats::quantile(x, probs = 0.025, names = FALSE)
#' @rdname q2.5
#' @export
q25   <- function(x) stats::quantile(x, probs = 0.25,  names = FALSE)
#' @rdname q2.5
#' @export
q50   <- function(x) stats::quantile(x, probs = 0.50,  names = FALSE)
#' @rdname q2.5
#' @export
q75   <- function(x) stats::quantile(x, probs = 0.75,  names = FALSE)
#' @rdname q2.5
#' @export
q97.5 <- function(x) stats::quantile(x, probs = 0.975, names = FALSE)

#' Matern (nu = 3/2) covariance
#'
#' Matern covariance with smoothness 3/2, in the same parameterization as the
#' Stan models: `C(d) = sigma^2 (1 + sqrt(3) d / ell) exp(-sqrt(3) d / ell)`.
#' The argument name `rho` (the length-scale `ell`) is kept for backward
#' compatibility.
#'
#' @param d Distances.
#' @param sigma Marginal standard deviation.
#' @param rho Length-scale `ell`; not the density-expression correlation.
#' @return Numeric vector of covariances, one per element of `d`.
#' @export
#' @examples
#' matern32(c(0, 0.5, 1, 2), rho = 1)
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
  if (!is.numeric(dataset[[gene]]))
    stop("Column '", gene, "' is not numeric: pass `gene =`, and drop non-numeric ",
         "columns such as a cell-type label.")
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

# Arguments removed in 0.9.0, when the package kept only the DenMark model:
# fail before any gridding or compilation instead of inside $sample().
.denmark_removed_args <- function(...) {
  bad <- intersect(...names(), c("model", "waic", "fit_simple"))
  if (length(bad))
    stop("`", bad[1], "` was removed in DenMark 0.9.0: the package fits only the ",
         "DenMark model and no longer computes WAIC.", call. = FALSE)
}
