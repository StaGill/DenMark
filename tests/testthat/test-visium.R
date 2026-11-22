# A hexagonal Visium-like lattice, pitch 100: rows alternate by half a pitch.
hex_spots <- function(nx = 14, ny = 12, seed = 1) {
  xy <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  set.seed(seed)
  data.frame(x = 100 * (xy$i + 0.5 * (xy$j %% 2)),
             y = 100 * sqrt(3) / 2 * xy$j,
             Y = stats::rpois(nrow(xy), 8),
             UMI = stats::rpois(nrow(xy), 5000),
             geneA = stats::rpois(nrow(xy), 3),
             geneB = stats::rpois(nrow(xy), 1))
}

# 3 x 2 grid on [0, 3] x [0, 2]: grid cell k holds k cells, each with 2 counts
# of the gene and 10 total counts.
umi_grid_df <- function() {
  k <- rep(1:6, times = 1:6)
  data.frame(x = (k - 1) %% 3 + 0.5, y = (k - 1) %/% 3 + 0.5, g = 2, u = 10)
}

test_that("buildgridpp sums the total counts per grid cell", {
  g <- buildgridpp(c(3, 2), umi_grid_df(), gene = "g", bbox = c(0, 3, 0, 2), umi = "u")
  expect_equal(g$grid$T, 10 * (1:6))
  expect_equal(g$grid$M, 2 * (1:6))
  g0 <- buildgridpp(c(3, 2), umi_grid_df(), gene = "g", bbox = c(0, 3, 0, 2))
  expect_null(g0$grid$T)
})

test_that("the umi column is not taken as the default expression column", {
  df <- umi_grid_df()[, c("x", "y", "u", "g")]
  g <- buildgridpp(c(3, 2), df, bbox = c(0, 3, 0, 2), umi = "u")
  expect_equal(g$gene, "g")
})

test_that("offset = 'umi' gives log(T + 0.05) and a per-grid area offset", {
  df <- umi_grid_df()
  p <- denmark_prepare(df, gene = "g", grid_size = c(3, 2), basis = 3,
                       bbox = c(0, 3, 0, 2), offset = "umi", umi = "u")
  sd <- p$stan_data
  expect_equal(sd$log_umi, log(10 * (1:6) + 0.05))
  expect_length(sd$log_grid_area, sd$n)
  expect_equal(sd$log_grid_area, rep(log(1), 6))
  expect_identical(p$info$offset, "umi")
  expect_false(p$info$units)

  p0 <- denmark_prepare(df, gene = "g", grid_size = c(3, 2), basis = 3, bbox = c(0, 3, 0, 2))
  expect_null(p0$stan_data$log_umi)
  expect_length(p0$stan_data$log_grid_area, 1)
  expect_identical(p0$info$offset, "cells")
})

test_that("offset = 'umi' needs the umi column", {
  df <- umi_grid_df()
  expect_error(denmark_prepare(df, gene = "g", grid_size = 3, offset = "umi"), "umi =")
  expect_error(denmark_prepare(df, gene = "g", grid_size = 3, offset = "umi", umi = "nope"),
               "not a column")
  df$u[1] <- -1
  expect_error(denmark_prepare(df, gene = "g", grid_size = 3, offset = "umi", umi = "u"),
               "non-negative")
})

test_that("denmark_init starts beta1 at the counts per unit of total count", {
  p <- denmark_prepare(umi_grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                       bbox = c(0, 3, 0, 2), offset = "umi", umi = "u")
  i <- denmark_init(p$stan_data)()
  expect_equal(i$beta1, log(sum(p$stan_data$y2) / sum(10 * (1:6) + 0.05)))
})

test_that("the umi offset selects DenMark-vararea-umi.stan", {
  old <- options(denmark.stan_dir = test_path("..", "..", "inst", "stan"))
  on.exit(options(old))
  expect_equal(basename(.denmark_stan_path(umi = TRUE)), "DenMark-vararea-umi.stan")
  expect_equal(basename(.denmark_stan_path(vararea = TRUE)), "DenMark-vararea.stan")
  expect_equal(basename(.denmark_stan_path()), "DenMark.stan")
})

test_that("DenMark refuses a custom Stan file without log_umi for offset = 'umi'", {
  f <- tempfile(fileext = ".stan")
  writeLines(c("data {", "  int n;", "  vector[n] log_grid_area;", "}"), f)
  expect_error(DenMark(umi_grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                       bbox = c(0, 3, 0, 2), offset = "umi", umi = "u", stan_file = f),
               "log_umi")
})

test_that("denmark_pool_visium pools seven spots and keeps every count", {
  spots <- hex_spots()
  u <- denmark_pool_visium(spots)
  sid <- attr(u, "spot_unit")
  expect_length(sid, nrow(spots))
  expect_setequal(unique(sid), u$unit)
  expect_equal(sum(u$Y), sum(spots$Y))
  expect_equal(sum(u$UMI), sum(spots$UMI))
  expect_equal(sum(u$geneA), sum(spots$geneA))
  expect_equal(sum(u$geneB), sum(spots$geneB))
  expect_equal(sum(u$n_spots), nrow(spots))
  expect_equal(u$area, u$n_spots * pi * 27.5^2)
  expect_true(all(u$n_spots >= 3))
})

test_that("before folding, pools hold at most seven spots, and a full pool is a spot and its six neighbours", {
  spots <- hex_spots()
  u <- denmark_pool_visium(spots, min_spots = 1)
  sid <- attr(u, "spot_unit")
  expect_true(all(u$n_spots <= 7))
  full <- u$unit[u$n_spots == 7]
  expect_gte(length(full), 5)
  for (k in full) {
    s <- spots[sid == k, c("x", "y")]
    d <- sqrt((s$x - mean(s$x))^2 + (s$y - mean(s$y))^2)
    expect_equal(sort(d), c(0, rep(100, 6)), tolerance = 1e-8)
  }
})

test_that("pools smaller than min_spots are folded into a larger pool", {
  spots <- hex_spots()
  u3 <- denmark_pool_visium(spots, min_spots = 3)
  u1 <- denmark_pool_visium(spots, min_spots = 1)
  expect_true(any(u1$n_spots < 3))
  expect_lt(nrow(u3), nrow(u1))
  expect_equal(sum(u3$n_spots), sum(u1$n_spots))
})

test_that("denmark_pool_visium checks its input", {
  spots <- hex_spots()
  expect_error(denmark_pool_visium(spots[, -3]), "Missing column")
  spots$geneA[1] <- 0.5
  expect_error(denmark_pool_visium(spots), "integer counts")
  expect_error(denmark_pool_visium(hex_spots(), genes = "nope"), "not in")
  u <- denmark_pool_visium(hex_spots()[, c("x", "y", "Y", "geneA")], umi = NULL)
  expect_named(u, c("unit", "x", "y", "Y", "n_spots", "area", "geneA"))
})

test_that("units = TRUE builds the Stan data of the spot-level analysis", {
  u <- denmark_pool_visium(hex_spots())
  p <- denmark_prepare(u, gene = "geneA", basis = 4, boundfactor = 1.5,
                       units = TRUE, offset = "umi", umi = "UMI")
  sd <- p$stan_data
  expect_equal(sd$n, nrow(u))
  expect_equal(sd$y1, u$Y)
  expect_equal(sd$y2, u$geneA)
  expect_equal(sd$log_grid_area, log(u$area))
  expect_equal(sd$log_umi, log(u$UMI + 0.05))
  cx <- u$x - mean(range(u$x)); cy <- u$y - mean(range(u$y))
  expect_equal(sd$coords, cbind(cx, cy), ignore_attr = TRUE)
  expect_equal(sd$L, 1.5 * c(max(abs(cx)), max(abs(cy))))
  expect_equal(sd$mrange, max(stats::dist(cbind(cx, cy))) / (2 * 2.75))
  expect_equal(sd$mstar, 16)
  expect_true(p$info$units)
  expect_true(p$info$vararea)
  expect_true(all(is.na(p$grid$i)))

  p_cells <- denmark_prepare(u, gene = "geneA", basis = 4, units = TRUE)
  expect_null(p_cells$stan_data$log_umi)
  expect_error(denmark_prepare(u, basis = 4, units = TRUE), "gene =")
  expect_error(denmark_prepare(u, gene = "geneA", basis = 4, units = TRUE, area = "nope"),
               "Missing column")
})

test_that("denmark_scan keeps the offset and unit columns", {
  seen <- list()
  local_mocked_bindings(
    DenMark = function(dataset, gene = NULL, ...) {
      seen[[gene]] <<- names(dataset)
      list(rho = seq(-0.2, 0.4, length.out = 50), stan_data = list(y2 = dataset[[gene]]),
           diagnostics = list(max_rhat = 1, min_ess_bulk = 400, num_divergent = 0),
           time = 0.01)
    },
    denmark_draws = function(obj) list(rho = obj$rho))
  cells <- data.frame(x = c(0, 1, 2), y = c(0, 1, 0), UMI = c(10, 12, 9),
                      geneA = c(1, 2, 0), geneB = c(0, 1, 1))
  res <- denmark_scan(cells, offset = "umi", umi = "UMI", verbose = FALSE)
  expect_setequal(res$gene, c("geneA", "geneB"))
  expect_equal(seen$geneA, c("x", "y", "UMI", "geneA"))
  units <- data.frame(unit = 1:3, x = c(0, 1, 2), y = c(0, 1, 0), Y = c(3, 5, 2),
                      UMI = c(10, 12, 9), n_spots = c(7, 7, 4), area = c(1, 1, 0.6),
                      geneA = c(1, 2, 0))
  res <- denmark_scan(units, units = TRUE, offset = "umi", umi = "UMI", verbose = FALSE)
  expect_equal(res$gene, "geneA")
  expect_equal(seen$geneA, c("x", "y", "UMI", "Y", "area", "geneA"))
  expect_error(denmark_scan(cells, units = TRUE, verbose = FALSE), "Missing column")
})

test_that("denmark_fields uses the total-count offset of an offset = 'umi' fit", {
  w <- function(v, x) posterior::as_draws_matrix(
    matrix(x, nrow = 2, byrow = TRUE, dimnames = list(NULL, paste0(v, "[", 1:3, "]"))))
  W <- list(w1 = w("w1", c(0.1, 0.2, 0.3, -0.1, 0, 0.1)),
            w2 = w("w2", c(0, 0.1, 0.2, 0.2, 0.1, 0)))
  local_mocked_bindings(denmark_draws = function(obj) data.frame(
    beta0 = c(0, 0.1), beta1 = c(-1, -0.9), a11 = c(1, 1.1), a21 = c(0.3, 0.2), a22 = c(0.5, 0.6)))
  obj <- structure(list(
    fit = list(draws = function(v) W[[v]]),
    stan_data = list(n = 3, log_grid_area = rep(0, 3), log_umi = log(c(100, 200, 300) + 0.05)),
    grid = data.frame(grid_index = 1:3, i = NA_integer_, j = NA_integer_, x = c(0, 1, 2),
                      y = c(0, 1, 0), Y = c(3, 5, 2), M = c(1, 2, 0))), class = "denmark")
  f <- denmark_fields(obj)
  expect_equal(f$log_lambda2 - f$eta2, obj$stan_data$log_umi)
  obj$stan_data$log_umi <- NULL
  f <- denmark_fields(obj)
  expect_equal(f$log_lambda2 - f$eta2, log(obj$grid$Y + 0.05))
})

test_that("maps of pre-aggregated units draw one point per unit", {
  u <- data.frame(grid_index = 1:3, i = NA_integer_, j = NA_integer_, x = c(0, 2.5, 1),
                  y = c(0, 0.9, 2.6), Y = c(3, 5, 2), M = c(1, 2, 0), type = "unit")
  expect_s3_class(.tile_plot(u, "Y", "cells")$layers[[1]]$geom, "GeomPoint")
  g <- transform(u, i = 1:3, j = 1L)
  expect_s3_class(.tile_plot(g, "Y", "cells")$layers[[1]]$geom, "GeomTile")
  p <- plot_denmark_grid(list(grid = u, gene = "geneA"))
  expect_no_error(ggplot2::ggplot_build(p[[1]]))
})
