# Grid cell (i, j) of a 3 x 2 grid on [0, 3] x [0, 2] holds i + 3 (j - 1) cells,
# each with 2 counts, so counts and ordering can be checked exactly.
grid_df <- function() {
  k <- rep(1:6, times = 1:6)
  data.frame(x = (k - 1) %% 3 + 0.5, y = (k - 1) %/% 3 + 0.5, g = 2)
}

prep_df <- function() {
  set.seed(42)
  data.frame(x = stats::runif(400, 0, 6), y = stats::runif(400, 1, 5),
             geneA = stats::rpois(400, 2), geneB = stats::rpois(400, 1))
}

test_that("denmark_filter_genes drops low-variance and rarely detected genes", {
  set.seed(1)
  expr <- data.frame(good1 = stats::rpois(200, 5), good2 = stats::rpois(200, 2),
                     rare = c(rep(5, 4), rep(0, 196)),
                     flat = c(2, rep(1, 199)))
  expect_equal(denmark_filter_genes(expr), c("good1", "good2"))
  expect_equal(denmark_filter_genes(expr, var_quantile = 0, min_nonzero = 0),
               names(expr))
  expect_error(denmark_filter_genes(matrix(1:4, 2)), "gene names")
})

test_that("buildgridpp counts cells and orders grids with x varying fastest", {
  df <- grid_df()
  g <- buildgridpp(c(3, 2), df, gene = "g", bbox = c(0, 3, 0, 2))
  expect_equal(g$grid_size, c(nx = 3L, ny = 2L))
  expect_equal(g$grid$Y, 1:6)
  expect_equal(g$grid$M, 2 * (1:6))
  expect_equal(g$grid$i, rep(1:3, 2))
  expect_equal(g$grid$j, rep(1:2, each = 3))
  expect_equal(g$grid$grid_index, g$grid$i + (g$grid$j - 1) * 3)
  expect_equal(dim(g$Ymat), c(3L, 2L))
  expect_equal(as.vector(g$Ymat), g$grid$Y)
  expect_equal(as.vector(g$Mmat), g$grid$M)
  expect_equal(unname(g$loc_mark_mat[, "point"]), g$grid$Y)
  expect_equal(unname(g$loc_mark_mat[, "mark"]), g$grid$M)
  expect_equal(sum(g$grid$Y), nrow(df))
  expect_equal(g$grid$x, rep(c(0.5, 1.5, 2.5), 2))
  expect_equal(g$grid$y, rep(c(0.5, 1.5), each = 3))
  expect_equal(g$grid_area, 1)
})

test_that("buildgridpp drops cells outside the bounding box", {
  df <- grid_df()
  g <- buildgridpp(c(2, 2), df, gene = "g", bbox = c(0, 2, 0, 2))
  expect_equal(sum(g$grid$Y), sum(df$x <= 2))
  expect_equal(g$grid$Y, c(1, 2, 4, 5))
})

test_that("buildgridpp keeps cells on the outer boundary of the data range", {
  df <- data.frame(x = c(0, 1, 2, 3), y = c(0, 3, 1, 2), g = c(1, 0, 2, 5))
  g <- buildgridpp(3, df, gene = "g")
  expect_equal(sum(g$grid$Y), 4)
  expect_equal(sum(g$grid$M), 8)
})

test_that("buildgridpp checks its input", {
  df <- grid_df()
  expect_error(buildgridpp(3, df, gene = "nope"), "not a column")
  df$g[1] <- -1
  expect_error(buildgridpp(3, df, gene = "g"), "non-negative")
  expect_error(buildgridpp(1, grid_df(), gene = "g"), "grid_size")
})

test_that("grid_pts_centroid matches the buildgridpp centroids", {
  df <- grid_df()
  gp <- grid_pts_centroid(c(3, 2), df, bbox = c(0, 3, 0, 2))
  g <- buildgridpp(c(3, 2), df, gene = "g", bbox = c(0, 3, 0, 2))
  expect_equal(gp$N, 6)
  expect_equal(gp$grid_area, 1)
  expect_equal(gp$step, c(dx = 1, dy = 1))
  expect_equal(gp$grid_points$x, g$grid$x)
  expect_equal(gp$grid_points$y, g$grid$y)
})

test_that("denmark_hsgp_indices lists all basis pairs as an integer matrix", {
  S <- denmark_hsgp_indices(c(3, 2))
  expect_true(is.matrix(S))
  expect_type(S, "integer")
  expect_equal(dim(S), c(6L, 2L))
  expect_null(dimnames(S))
  expect_equal(S[, 1], rep(1:3, 2))
  expect_equal(S[, 2], rep(1:2, each = 3))
  expect_equal(dim(denmark_hsgp_indices(10)), c(100L, 2L))
})

test_that("denmark_prepare builds the Stan data list", {
  df <- prep_df()
  prep <- denmark_prepare(df, gene = "geneA", grid_size = 10, basis = 10,
                          boundfactor = 1.5)
  sd <- prep$stan_data
  expect_named(sd, c("n", "y1", "y2", "d", "mstar", "coords", "log_grid_area",
                     "mrange", "indices", "L", "is_centerted_PHI"))
  expect_equal(sd$n, 100)
  expect_equal(sd$mstar, 100)
  expect_identical(sd$d, 2L)
  expect_identical(sd$is_centerted_PHI, 0L)
  expect_type(sd$y1, "integer")
  expect_type(sd$y2, "integer")
  expect_length(sd$y1, 100)
  expect_type(sd$coords, "double")
  expect_equal(dim(sd$coords), c(100L, 2L))
  expect_type(sd$indices, "integer")
  expect_equal(dim(sd$indices), c(100L, 2L))
  expect_length(sd$L, 2)
  expect_equal(sum(sd$y1), nrow(df))
  expect_equal(sum(sd$y2), sum(df$geneA))
  expect_equal(sd$log_grid_area, log(prep$info$grid_area))
  expect_equal(sd$mrange, prep$info$max_dist / 5.5)
  expect_equal(prep$info$gene, "geneA")
  expect_equal(nrow(prep$grid), 100)
})

test_that("denmark_prepare centers the grid and sets L from the bounding box", {
  df <- prep_df()
  prep <- denmark_prepare(df, gene = "geneA", grid_size = c(12, 8), basis = 8,
                          boundfactor = 1.5)
  sd <- prep$stan_data
  half <- c(diff(range(df$x)), diff(range(df$y))) / 2
  expect_equal(unname(colMeans(sd$coords)), c(0, 0), tolerance = 1e-10)
  expect_equal(range(sd$coords[, 1]), c(-1, 1) * (half[1] - prep$info$step[["dx"]] / 2))
  expect_equal(range(sd$coords[, 2]), c(-1, 1) * (half[2] - prep$info$step[["dy"]] / 2))
  expect_equal(sd$L, 1.5 * half)
  expect_equal(c(prep$info$L1, prep$info$L2), half)
  expect_equal(unname(prep$info$centre), c(mean(range(df$x)), mean(range(df$y))))
  expect_equal(sd$mstar, 64)
  expect_equal(prep$grid$cx, sd$coords[, 1])
})

test_that("denmark_prepare checks half-ranges and the boundary factor", {
  df <- prep_df()
  expect_warning(denmark_prepare(df, gene = "geneA", basis = 5, L1 = 100),
                 "supplied but the data give")
  expect_error(denmark_prepare(df, gene = "geneA", boundfactor = 0.9),
               "boundfactor")
})

# A tissue polygon covering x in [0, 1.5] of the 3 x 2 grid on [0, 3] x [0, 2]:
# column 1 is inside, column 2 half inside, column 3 outside.
half_tissue <- function() {
  sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1.5, 0), c(1.5, 2), c(0, 2), c(0, 0)))))
}

test_that("buildgridpp clips grid areas to the tissue and classifies grid cells", {
  skip_if_not_installed("sf")
  g <- buildgridpp(c(3, 2), grid_df(), gene = "g", bbox = c(0, 3, 0, 2),
                   boundary = half_tissue())
  expect_equal(g$grid$type, rep(c("interior", "boundary", "exterior"), 2))
  expect_equal(g$grid$area, rep(c(1, 0.5, 0), 2))
  expect_equal(g$grid$Y, 1:6)
  expect_s3_class(g$boundary, "sfc")

  plain <- buildgridpp(c(3, 2), grid_df(), gene = "g", bbox = c(0, 3, 0, 2))
  expect_equal(plain$grid$area, rep(1, 6))
  expect_equal(plain$grid$type, rep("interior", 6))
  expect_null(plain$boundary)
})

test_that("buildgridpp spans the extent of the boundary by default", {
  skip_if_not_installed("sf")
  expect_message(g <- buildgridpp(c(3, 2), grid_df(), gene = "g", boundary = half_tissue()),
                 "9 cell(s) lie outside the extent of the tissue boundary", fixed = TRUE)
  expect_equal(unname(g$bbox), c(0, 1.5, 0, 2))
  expect_true(all(g$grid$type == "interior"))
  far <- transform(grid_df(), x = x + 100)
  expect_error(suppressMessages(buildgridpp(3, far, gene = "g", boundary = half_tissue())),
               "same units")
})

test_that("buildgridpp refuses a boundary made of several features", {
  skip_if_not_installed("sf")
  square <- sf::st_sfc(sf::st_polygon(list(rbind(c(2, 0), c(3, 0), c(3, 1), c(2, 1), c(2, 0)))))
  expect_error(buildgridpp(c(3, 2), grid_df(), gene = "g", bbox = c(0, 3, 0, 2),
                           boundary = c(half_tissue(), square)),
               "has 2 features")
  merged <- buildgridpp(c(3, 2), grid_df(), gene = "g", bbox = c(0, 3, 0, 2),
                        boundary = sf::st_union(c(half_tissue(), square)))
  expect_equal(merged$grid$area, c(1, 0.5, 1, 1, 0.5, 0))
})

test_that("a non-numeric column is refused as the expression", {
  df <- data.frame(x = 1:3, y = 1:3, celltype = "T", g = 1:3)
  expect_error(buildgridpp(3, df), "Column 'celltype' is not numeric")
  expect_equal(sum(buildgridpp(3, df, gene = "g")$grid$M), 6)
})

test_that("denmark_prepare drops exterior grid cells and uses per-grid areas", {
  skip_if_not_installed("sf")
  expect_message(
    prep <- denmark_prepare(grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                            bbox = c(0, 3, 0, 2), boundary = half_tissue()),
    "9 cell(s) fall in grid cells outside the tissue", fixed = TRUE)
  sd <- prep$stan_data
  expect_equal(sd$n, 4)
  expect_equal(sd$y1, c(1L, 2L, 4L, 5L))
  expect_equal(sd$log_grid_area, log(c(1, 0.5, 1, 0.5)))
  expect_equal(sd$coords, cbind(c(-0.5, 0.5, -0.5, 0.5), c(-0.5, -0.5, 0.5, 0.5)))
  expect_equal(sd$L, 1.5 * c(1, 1))
  expect_true(prep$info$vararea)
  expect_equal(c(prep$info$n_grid_total, prep$info$n_interior, prep$info$n_boundary),
               c(6, 2, 2))
  expect_equal(prep$info$n_cells, 12)

  small <- suppressMessages(
    denmark_prepare(grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                    bbox = c(0, 3, 0, 2), boundary = half_tissue(), min_area_frac = 0.6))
  expect_equal(small$stan_data$n, 2)
  expect_equal(small$stan_data$log_grid_area, c(0, 0))
})

test_that("denmark_init starts beta0 at the number of cells per unit area", {
  skip_if_not_installed("sf")
  prep <- suppressMessages(
    denmark_prepare(grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                    bbox = c(0, 3, 0, 2), boundary = half_tissue()))
  expect_equal(denmark_init(prep$stan_data)()$beta0, log(12 / 3))
  plain <- denmark_prepare(grid_df(), gene = "g", grid_size = c(3, 2), basis = 3,
                           bbox = c(0, 3, 0, 2))
  expect_equal(denmark_init(plain$stan_data)()$beta0, log(21 / 6))
})

test_that("denmark_alpha_shape follows a concave tissue and records alpha", {
  skip_if_not_installed("sf")
  p <- expand.grid(x = seq(0, 1, by = 0.05), y = seq(0, 1, by = 0.05))
  L_shape <- p[!(p$x > 0.5 & p$y > 0.5), ]
  # default alpha: 3 x the median Delaunay edge, which is 0.05 on this lattice
  expect_equal(attr(denmark_alpha_shape(L_shape$x, L_shape$y), "alpha"), 0.15)
  # alpha = 0.04 keeps the lattice triangles (circumradius 0.035) and the one
  # triangle that cuts the concave corner (half a lattice cell): the notch stays open
  tissue <- denmark_alpha_shape(L_shape$x, L_shape$y, alpha = 0.04)
  expect_s3_class(tissue, "sfc")
  expect_equal(attr(tissue, "alpha"), 0.04)
  corner <- 0.05^2 / 2
  expect_equal(as.numeric(sf::st_area(tissue)), 0.75 + corner, tolerance = 1e-10)

  blob <- data.frame(x = p$x / 4 + 3, y = p$y / 4)     # a separate 0.25 x 0.25 piece
  both <- rbind(L_shape, blob)
  expect_equal(as.numeric(sf::st_area(denmark_alpha_shape(both$x, both$y, alpha = 0.04))),
               0.75 + corner, tolerance = 1e-10)
  expect_equal(as.numeric(sf::st_area(denmark_alpha_shape(both$x, both$y, alpha = 0.04,
                                                          keep = "all"))),
               0.75 + corner + 1 / 16, tolerance = 1e-10)
  expect_error(denmark_alpha_shape(c(0, 1), c(0, 1)), "At least 3")
})
