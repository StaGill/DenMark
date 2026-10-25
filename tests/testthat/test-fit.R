# A small fit on the toy data (Tumor cells inside the alpha-shape tissue
# boundary of all cells), computed once and shared by the tests.
toy_tumor <- function() subset(denmark_toy_data(), celltype == "Tumor", select = -celltype)
toy_tissue <- function() {
  toy <- denmark_toy_data()
  denmark_alpha_shape(toy$x, toy$y)
}
toy_fits <- new.env()
toy_fit <- function() {
  if (is.null(toy_fits$fit))
    toy_fits$fit <- DenMark(toy_tumor(), gene = "geneA", boundary = toy_tissue(),
                            grid_size = 10, basis = 10, chains = 2,
                            iter_warmup = 300, iter_sampling = 300,
                            seed = 2024, quiet = TRUE)
  toy_fits$fit
}

test_that("DenMarkmain stops when asked for M1", {
  cells <- data.frame(x = 1:3, y = 1:3, g = 1:3)
  expect_error(DenMarkmain(grid_size = 10, dataset = cells, basis = 10,
                           boundfactor = 1.5, denmark_version = "M1"),
               "was removed in DenMark 0.9.0", fixed = TRUE)
})

test_that("arguments of the removed models fail before any fitting", {
  cells <- data.frame(x = 1:3, y = 1:3, g = 1:3)
  expect_error(DenMark(cells, gene = "g", model = "M2"), "`model` was removed", fixed = TRUE)
  expect_error(DenMark(cells, gene = "g", waic = TRUE), "`waic` was removed", fixed = TRUE)
  expect_error(denmark_scan(cells, fit_simple = TRUE), "`fit_simple` was removed", fixed = TRUE)
  expect_error(denmark_scan(cells, workers = 0), "positive whole number")
  expect_error(denmark_scan(cells, workers = 1.5), "positive whole number")
})

test_that("a custom Stan file without per-grid areas is refused with a boundary", {
  skip_if_not_installed("sf")
  skip_if_no_toy_data()
  expect_error(DenMark(toy_tumor(), gene = "geneA", boundary = toy_tissue(),
                       grid_size = 10, basis = 5, stan_file = .denmark_stan_path()),
               "needs `vector[n] log_grid_area`", fixed = TRUE)
})

test_that("the packaged model compiles into a per-version user cache directory", {
  skip_if_no_cmdstan()
  mod <- denmark_model()
  cache <- file.path(tools::R_user_dir("DenMark", which = "cache"),
                     as.character(utils::packageVersion("DenMark")))
  expect_equal(normalizePath(dirname(mod$exe_file())), normalizePath(cache))
  expect_identical(denmark_model(), mod)
})

test_that("DenMark fits the toy data inside the tissue boundary", {
  skip_if_no_cmdstan()
  skip_if_not_installed("sf")
  skip_if_no_toy_data()
  fit <- toy_fit()
  expect_s3_class(fit, "denmark")
  expect_true(fit$info$vararea)
  expect_equal(fit$stan_data$n, 82)
  expect_length(fit$stan_data$log_grid_area, fit$stan_data$n)
  rho <- denmark_draws(fit)$rho
  expect_length(rho, 600)
  expect_true(all(is.finite(rho)))
  expect_true(all(abs(rho) < 1))
  s <- denmark_summary(fit)
  expect_equal(s$variable, c("beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]",
                             "sigma[1]", "sigma[2]"))
  expect_true(all(c("mean", "sd", "q2.5", "q97.5", "rhat", "ess_bulk") %in% names(s)))
  expect_true(all(is.finite(s$mean)))
  expect_gt(mean(rho), 0)
  expect_output(print(fit), "rho")
})

test_that("DenMark without a boundary fits the bounding box with a common grid area", {
  skip_if_no_cmdstan()
  skip_if_no_toy_data()
  fit <- DenMark(toy_tumor(), gene = "geneA", grid_size = 10, basis = 5,
                 chains = 1, iter_warmup = 150, iter_sampling = 150,
                 seed = 2024, quiet = TRUE)
  expect_false(fit$info$vararea)
  expect_equal(fit$stan_data$n, 100)
  expect_length(fit$stan_data$log_grid_area, 1)
  expect_length(denmark_draws(fit)$rho, 150)
})

test_that("post-processing works on a DenMark fit", {
  skip_if_no_cmdstan()
  skip_if_not_installed("sf")
  skip_if_no_toy_data()
  fit <- toy_fit()
  fl <- denmark_fields(fit)
  expect_equal(nrow(fl), fit$stan_data$n)
  expect_true(all(c("eta1", "a11w1", "eta2", "a21w1", "a22w2", "log_lambda1",
                    "log_lambda2") %in% names(fl)))
  dc <- denmark_decay(fit)
  expect_equal(dc$ranges$field, c("density", "expression"))
  expect_true(all(dc$ranges$mean > 0))
  expect_equal(nrow(dc$range_draws), 600)
  expect_equal(nrow(denmark_check_hsgp(fit, warn = FALSE)), 4)

  f <- tempfile(fileext = ".rds")
  denmark_save(fit, f)
  expect_equal(denmark_draws(denmark_load(f))$rho, denmark_draws(fit)$rho)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(print(plot_denmark_fields(fit)))
  expect_no_error(print(plot_denmark_grid(fit)))
})

test_that("denmark_scan fits, classifies and resumes a two-gene scan", {
  skip_if_no_cmdstan()
  skip_if_not_installed("sf")
  skip_if_no_toy_data()
  tumor <- toy_tumor()
  tissue <- toy_tissue()
  dir <- tempfile("scan")
  run <- function() {
    denmark_scan(tumor, genes = c("geneA", "geneC"), boundary = tissue, save_dir = dir,
                 grid_size = 10, basis = 10, chains = 2, iter_warmup = 300,
                 iter_sampling = 300, seed = 2024, verbose = FALSE)
  }
  res <- run()
  expect_equal(nrow(res), 2)
  expect_setequal(res$gene, c("geneA", "geneC"))
  expect_true(all(is.na(res$error)))
  expect_true(all(res$class %in% c("PCG", "NCG", "not DCG")))
  expect_true(all(res$rho_lower <= res$rho_mean & res$rho_mean <= res$rho_upper))
  expect_gt(res$rho_mean[res$gene == "geneA"], res$rho_mean[res$gene == "geneC"])
  expect_equal(attr(res, "cred_level"), 0.95)
  expect_true(file.exists(file.path(dir, "denmark_geneA.rds")))

  again <- run()
  expect_equal(again$rho_mean, res$rho_mean)
})

test_that("a parallel scan gives the same results as the sequential scan", {
  skip_if_no_cmdstan()
  skip_if_not_installed("sf")
  skip_if_not_installed("future.apply")
  skip_if(length(find.package("DenMark", lib.loc = .libPaths(), quiet = TRUE)) == 0,
          "parallel workers need an installed DenMark")
  tumor <- toy_tumor()
  tissue <- toy_tissue()
  run <- function(workers) {
    denmark_scan(tumor, genes = c("geneA", "geneC"), boundary = tissue, workers = workers,
                 grid_size = 8, basis = 5, chains = 1, iter_warmup = 150,
                 iter_sampling = 150, seed = 2024, verbose = FALSE)
  }
  cols <- c("gene", "rho_mean", "rho_lower", "rho_upper", "max_rhat", "divergences", "class")
  expect_equal(run(2)[cols], run(1)[cols])
})
