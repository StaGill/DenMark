test_that("matern32 matches the Matern-3/2 covariance", {
  d <- c(0, 0.5, 1, 2)
  expect_equal(matern32(d, sigma = 2, rho = 1.5),
               4 * (1 + sqrt(3) * d / 1.5) * exp(-sqrt(3) * d / 1.5))
  expect_equal(matern32(0, sigma = 1.3, rho = 1), 1.69)
  expect_equal(matern32(1, rho = 1), (1 + sqrt(3)) * exp(-sqrt(3)))
})

test_that("the practical range factor at 0.05 is about 2.7389", {
  k <- .matern32_range_factor(0.05)
  expect_equal(k, 2.7389, tolerance = 1e-4)
  expect_equal(matern32(k, rho = 1), 0.05, tolerance = 1e-8)
})

test_that("quantile helpers return unnamed sample quantiles", {
  x <- 0:100
  expect_equal(c(q2.5(x), q25(x), q50(x), q75(x), q97.5(x)),
               c(2.5, 25, 50, 75, 97.5))
  expect_null(names(q50(x)))
  expect_equal(q97.5(c(3, 1, 2)), unname(stats::quantile(c(3, 1, 2), 0.975)))
})

test_that("the Stan models ship with the package", {
  stan_dir <- system.file("stan", package = "DenMark")
  expect_true(nzchar(stan_dir))
  expect_setequal(list.files(stan_dir, pattern = "\\.stan$"),
                  c("DenMark.stan", "DenMark-vararea.stan", "DenMark-vararea-umi.stan"))
})

test_that(".denmark_stan_path finds the shipped model and honors overrides", {
  expect_match(.denmark_stan_path(), "DenMark\\.stan$")
  expect_true(file.exists(.denmark_stan_path()))
  expect_match(.denmark_stan_path(vararea = TRUE), "DenMark-vararea\\.stan$")
  expect_true(file.exists(.denmark_stan_path(vararea = TRUE)))

  tmp <- tempfile("stan_dir")
  dir.create(tmp)
  file.copy(.denmark_stan_path(), tmp)
  old <- options(denmark.stan_dir = tmp)
  p <- tryCatch(.denmark_stan_path(), finally = options(old))
  expect_equal(normalizePath(p), normalizePath(file.path(tmp, "DenMark.stan")))

  expect_error(.denmark_stan_path(stan_file = file.path(tmp, "missing.stan")),
               "not found")
})
