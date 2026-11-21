test_that(".denmark_max_dist is the largest pairwise distance", {
  set.seed(3)
  pts <- cbind(stats::runif(300, -2, 5), stats::rnorm(300, 0, 3))
  expect_equal(.denmark_max_dist(pts), max(stats::dist(pts)))
  line <- cbind(seq(0, 9, by = 1), 2 * seq(0, 9, by = 1))           # collinear
  expect_equal(.denmark_max_dist(line), sqrt(9^2 + 18^2))
  expect_equal(.denmark_max_dist(matrix(c(1, 2), 1)), 0)
})

test_that("the prior scale uses the largest distance, not the bounding-box diagonal", {
  # a disc of points: the largest distance is the diameter, while the
  # diagonal of the bounding box is sqrt(2) times longer
  th <- seq(0, 2 * pi, length.out = 361)[-361]
  disc <- cbind(cos(th), sin(th))
  expect_equal(.denmark_max_dist(disc), 2, tolerance = 1e-12)
  expect_gt(sqrt(diff(range(disc[, 1]))^2 + diff(range(disc[, 2]))^2), 2.8)
})

test_that("denmark_prepare scales the length-scale prior by the largest centroid distance", {
  set.seed(42)
  df <- data.frame(x = stats::runif(400, 0, 6), y = stats::runif(400, 1, 5),
                   geneA = stats::rpois(400, 2))
  p <- denmark_prepare(df, gene = "geneA", grid_size = c(6, 4), basis = 5)
  expect_equal(p$info$max_dist, max(stats::dist(p$stan_data$coords)))
  expect_equal(p$stan_data$mrange, p$info$max_dist / (2 * 2.75))
})
