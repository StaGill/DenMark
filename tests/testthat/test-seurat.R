# A Seurat object with centroids, built from the first cells of the toy data.
toy_seurat <- function(n = 300) {
  toy <- denmark_toy_data()[seq_len(n), ]
  cells <- paste0("cell", seq_len(n))
  counts <- Matrix::Matrix(t(as.matrix(toy[, c("geneA", "geneB", "geneC")])), sparse = TRUE)
  colnames(counts) <- cells
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  obj$celltype <- toy$celltype
  obj[["fov"]] <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(data.frame(x = toy$x, y = toy$y, cell = cells)),
    type = "centroids", assay = "RNA")
  list(object = obj, toy = toy)
}

test_that("denmark_from_seurat returns the cell table of the object", {
  skip_if_not_installed("SeuratObject")
  skip_if_no_toy_data()
  s <- toy_seurat()
  cells <- denmark_from_seurat(s$object, celltype = "celltype")
  expect_named(cells, c("x", "y", "celltype", "geneA", "geneB", "geneC"))
  expect_equal(rownames(cells), paste0("cell", 1:300))
  rownames(cells) <- NULL
  expect_equal(cells, s$toy[, names(cells)], ignore_attr = TRUE)

  two <- denmark_from_seurat(s$object, celltype = "celltype", genes = c("geneC", "geneA"))
  expect_named(two, c("x", "y", "celltype", "geneC", "geneA"))
  expect_error(denmark_from_seurat(s$object, celltype = "nope"), "not a column")
  expect_error(denmark_from_seurat(s$object, celltype = "celltype", genes = "geneZ"),
               "Not features")
})

test_that("denmark_annotate writes the scan into the feature metadata", {
  skip_if_not_installed("SeuratObject")
  skip_if_no_toy_data()
  obj <- toy_seurat(50)$object
  scan <- data.frame(gene = c("geneC", "geneA"), rho_mean = c(-0.6, 0.6),
                     rho_lower = c(-0.8, 0.3), rho_upper = c(-0.3, 0.8), class = c("NCG", "PCG"))
  md <- denmark_annotate(obj, scan)[["RNA"]][[]]
  expect_equal(md[c("geneA", "geneB", "geneC"), "denmark_rho"], c(0.6, NA, -0.6))
  expect_equal(md[c("geneA", "geneB", "geneC"), "denmark_class"], c("PCG", NA, "NCG"))
  expect_true(all(c("denmark_rho_lower", "denmark_rho_upper") %in% names(md)))
})
