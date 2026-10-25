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

# Irregular, lobed tissue shape inside [0, width] x [0, height]
.tissue_radius <- function(theta)
  0.95 * (1 + 0.22 * sin(3 * theta + 0.5) + 0.12 * cos(2 * theta) +
            0.08 * sin(5 * theta + 1)) / 1.42

.in_tissue <- function(x, y, width, height) {
  u <- (x - width / 2) / (width / 2); v <- (y - height / 2) / (height / 2)
  sqrt(u^2 + v^2) < .tissue_radius(atan2(v, u))
}

#' Simulate a marked point pattern from the DenMark data-generating process
#'
#' Cells of the analyzed type follow an inhomogeneous Poisson process with
#' `log lambda1(s) = beta0 + a11 w1(s)`, and each cell's count for gene `k` is
#' Poisson with
#' `log mu_k(s) = beta1 + sigma2 * (rho_k w1(s) + sqrt(1 - rho_k^2) w2k(s))`.
#'
#' The fields are smooth, zero-mean, unit-variance Gaussian fields with
#' approximately Matern-3/2 covariance, drawn by random Fourier features. The
#' fields are standardized, and each `w2k` is orthogonalized against `w1`, on
#' a regular lattice over the tissue (that is, over space rather than over the
#' cells), so that the spatial correlation of the two fields equals `rho_k`,
#' which is what rho measures in the model.
#'
#' With `tissue = "irregular"`, the cells lie in a lobed tissue section inside
#' the window, so that the analysis needs a tissue boundary (see
#' [denmark_alpha_shape()]). The cells of other types are drawn from the
#' candidate locations that the thinning rejected, so they are densest where
#' the analyzed cells are sparse; their gene counts are Poisson with mean 0.3.
#'
#' @param rho Named vector of true correlations, one gene per element.
#' @param width,height Size of the window `[0, width] x [0, height]`.
#' @param n_cells Number of cells of the analyzed type: the thinned point
#'   pattern is subsampled to `n_cells` cells when it has more.
#' @param ell1,ell2 Length-scales of the density and expression fields.
#' @param a11 Standard deviation of the log cell intensity.
#' @param sigma2 Standard deviation of the log expression.
#' @param beta1 Intercept of the log expression.
#' @param seed Random seed, set with `set.seed()`; this resets the global
#'   random number stream.
#' @param tissue `"rectangle"` (cells fill the window) or `"irregular"` (cells
#'   lie in a lobed tissue section).
#' @param n_other Number of cells of other types; together with the analyzed
#'   cells they define the tissue boundary.
#' @param celltype Label of the analyzed cell type (used when `n_other > 0`;
#'   the other cells are labeled `"Other"`).
#' @return Data frame with the columns `x`, `y`, `celltype` (when
#'   `n_other > 0`) and one count column per gene; the true values are stored
#'   in `attr(, "truth")`, including the true tissue outline for
#'   `tissue = "irregular"`.
#' @export
#' @examples
#' sim <- denmark_simulate(rho = c(up = 0.5, down = -0.5), n_cells = 500)
#' head(sim)
#' attr(sim, "truth")$rho
#'
#' sec <- denmark_simulate(n_cells = 1000, tissue = "irregular", n_other = 800,
#'                         celltype = "Tumor")
#' table(sec$celltype)
denmark_simulate <- function(rho = c(geneA = 0.6, geneB = 0, geneC = -0.6),
                             width = 6, height = 4, n_cells = 5000,
                             ell1 = 1, ell2 = 1, a11 = 1, sigma2 = 0.8,
                             beta1 = 0.5, seed = 1,
                             tissue = c("rectangle", "irregular"),
                             n_other = 0, celltype = "Target") {
  tissue <- match.arg(tissue)
  if (is.null(names(rho))) names(rho) <- paste0("gene", seq_along(rho))
  set.seed(seed)
  w1 <- .rff_field(ell1)
  # candidate locations: uniform over the window, restricted to the tissue
  xs <- stats::runif(80000, 0, width); ys <- stats::runif(80000, 0, height)
  if (tissue == "irregular") {
    inside <- .in_tissue(xs, ys, width, height)
    xs <- xs[inside]; ys <- ys[inside]
  }
  # thinning: keep a candidate with probability lambda / max(lambda)
  f1 <- w1(xs, ys)
  lam <- exp(a11 * f1)
  keep <- stats::runif(length(lam)) < lam / max(lam)
  xs_o <- xs[!keep]; ys_o <- ys[!keep]                 # pool for other cells
  xs <- xs[keep]; ys <- ys[keep]; f1 <- f1[keep]
  if (length(xs) > n_cells) {
    s <- sample.int(length(xs), n_cells); xs <- xs[s]; ys <- ys[s]; f1 <- f1[s]
  }
  out <- data.frame(x = xs, y = ys)
  # regular lattice over the tissue: the fields are standardized and made
  # orthogonal over SPACE (area-weighted), which is what rho measures in the
  # model; two independent smooth fields can be correlated by chance on a
  # finite window, so w2 is orthogonalized against w1 on this lattice
  lat <- expand.grid(x = seq(0, width, length.out = 150),
                     y = seq(0, height, length.out = 100))
  if (tissue == "irregular") lat <- lat[.in_tissue(lat$x, lat$y, width, height), ]
  l1 <- w1(lat$x, lat$y)
  m1 <- mean(l1); s1 <- stats::sd(l1)
  z1 <- (f1 - m1) / s1; zl1 <- (l1 - m1) / s1
  for (g in names(rho)) {
    w2f <- .rff_field(ell2)
    l2 <- w2f(lat$x, lat$y)
    b <- stats::cov(l2, zl1) / stats::var(zl1)
    r2 <- l2 - b * zl1
    w2 <- (w2f(xs, ys) - b * z1 - mean(r2)) / stats::sd(r2)
    eta <- beta1 + sigma2 * (rho[[g]] * z1 + sqrt(1 - rho[[g]]^2) * w2)
    out[[g]] <- stats::rpois(length(xs), exp(eta))
  }
  if (n_other > 0) {
    s <- sample.int(length(xs_o), min(n_other, length(xs_o)))
    other <- data.frame(x = xs_o[s], y = ys_o[s])
    for (g in names(rho)) other[[g]] <- stats::rpois(nrow(other), 0.3)
    out$celltype <- celltype; other$celltype <- "Other"
    out <- rbind(out, other)[, c("x", "y", "celltype", names(rho))]
  }
  outline <- NULL
  if (tissue == "irregular") {
    th <- seq(-pi, pi, length.out = 361)
    r <- .tissue_radius(th)
    outline <- data.frame(x = width / 2 * (1 + r * cos(th)),
                          y = height / 2 * (1 + r * sin(th)))
  }
  attr(out, "truth") <- list(rho = rho, ell1 = ell1, ell2 = ell2, a11 = a11,
                             sigma2 = sigma2, beta1 = beta1,
                             window = c(width = width, height = height),
                             tissue = tissue, outline = outline)
  out
}
