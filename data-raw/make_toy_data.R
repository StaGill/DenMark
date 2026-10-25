# Toy data with a known answer, and the R/Python parity fixtures.
#
# Run from the repository root:
#   Rscript data-raw/make_toy_data.R
#
# Writes
#   inst/extdata/denmark_toy.csv.gz   simulated tissue section: x, y, celltype,
#                                     geneA, geneB, geneC
#   python/tests/data/toy_stan_data.json
#       denmark_prepare()$stan_data for the Tumor cells and geneA, without a
#       tissue boundary
#   python/tests/data/toy_stan_data_boundary.json
#       the same with the alpha-shape tissue boundary of all cells, with the
#       alpha, the tissue area and the grid areas and types
# Python's results are compared with these in the parity tests.
#
# The analyzed cell type ("Tumor") carries three genes with known
# density-expression correlations: geneA 0.6, geneB 0, geneC -0.6. The "Other"
# cells only help define the tissue boundary. This is the simulation of
# Tutorial_DenMark.Rmd, where the analyzed type is labeled "Invasive_Tumor".
# Needs pkgload, jsonlite and sf, which the package itself does not depend on.

stopifnot(file.exists("DESCRIPTION"), file.exists("data-raw/make_toy_data.R"))
pkgload::load_all(".", quiet = TRUE)

toy <- denmark_simulate(rho = c(geneA = 0.6, geneB = 0, geneC = -0.6),
                        width = 6, height = 4, n_cells = 5000, seed = 1,
                        tissue = "irregular", n_other = 4000, celltype = "Tumor")
toy$x <- round(toy$x, 6)
toy$y <- round(toy$y, 6)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
con <- gzfile("inst/extdata/denmark_toy.csv.gz", "wb")   # binary: LF line endings
utils::write.csv(toy, con, row.names = FALSE)
close(con)

# Stan data built from the file just written, as both packages read it
toy <- denmark_toy_data()
tumor <- subset(toy, celltype == "Tumor", select = -celltype)
dir.create("python/tests/data", recursive = TRUE, showWarnings = FALSE)
write_fixture <- function(x, file)
  jsonlite::write_json(x, file, auto_unbox = TRUE, digits = NA,
                       matrix = "rowmajor", pretty = FALSE)

prep <- denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10,
                        boundfactor = 1.5)
write_fixture(prep$stan_data, "python/tests/data/toy_stan_data.json")

tissue <- denmark_alpha_shape(toy$x, toy$y)
prep_b <- denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10,
                          boundfactor = 1.5, boundary = tissue)
write_fixture(list(alpha = attr(tissue, "alpha"),
                   tissue_area = as.numeric(sf::st_area(tissue)),
                   n_grid_total = prep_b$info$n_grid_total,
                   grid_area = prep_b$grid$area, grid_type = prep_b$grid$type,
                   grid_index = prep_b$grid$grid_index,
                   stan_data = prep_b$stan_data),
              "python/tests/data/toy_stan_data_boundary.json")

message("Wrote ", nrow(toy), " cells (", nrow(tumor), " Tumor) to ",
        "inst/extdata/denmark_toy.csv.gz and the geneA Stan data, without and ",
        "with the tissue boundary, to python/tests/data/")
