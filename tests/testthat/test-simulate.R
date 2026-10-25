test_that("denmark_simulate returns cells, counts and the truth", {
  sim <- denmark_simulate(rho = c(up = 0.5, down = -0.5), width = 3, height = 2,
                          n_cells = 300, seed = 3)
  expect_s3_class(sim, "data.frame")
  expect_named(sim, c("x", "y", "up", "down"))
  expect_gt(nrow(sim), 0)
  expect_lte(nrow(sim), 300)
  expect_true(all(sim$x >= 0 & sim$x <= 3 & sim$y >= 0 & sim$y <= 2))
  counts <- unlist(sim[c("up", "down")])
  expect_true(all(counts >= 0 & counts == round(counts)))

  truth <- attr(sim, "truth")
  expect_equal(truth$rho, c(up = 0.5, down = -0.5))
  expect_equal(truth$window, c(width = 3, height = 2))
  expect_equal(truth$sigma2, 0.8)

  again <- denmark_simulate(rho = c(up = 0.5, down = -0.5), width = 3, height = 2,
                            n_cells = 300, seed = 3)
  expect_identical(again, sim)
})

test_that("denmark_simulate names unnamed genes", {
  sim <- denmark_simulate(rho = c(0.3, -0.3), n_cells = 100, seed = 4)
  expect_named(sim, c("x", "y", "gene1", "gene2"))
  expect_equal(names(attr(sim, "truth")$rho), c("gene1", "gene2"))
})

test_that("denmark_simulate places cells in an irregular tissue with other cell types", {
  sim <- denmark_simulate(n_cells = 500, tissue = "irregular", n_other = 300,
                          celltype = "T", seed = 2)
  expect_named(sim, c("x", "y", "celltype", "geneA", "geneB", "geneC"))
  expect_equal(sum(sim$celltype == "Other"), 300)
  expect_lte(sum(sim$celltype == "T"), 500)
  expect_true(all(.in_tissue(sim$x, sim$y, 6, 4)))
  truth <- attr(sim, "truth")
  expect_equal(truth$tissue, "irregular")
  expect_equal(dim(truth$outline), c(361L, 2L))
  expect_null(attr(denmark_simulate(n_cells = 100, seed = 4), "truth")$outline)
})
