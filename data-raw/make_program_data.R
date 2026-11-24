# Example section for the article "Gene programs on Visium HD and Xenium".
#
# Run from the repository root:
#   Rscript data-raw/make_program_data.R
#
# Writes inst/extdata/denmark_programs.csv.gz: one row per cell with x, y,
# celltype, the total count of the cell (total) and the counts of the genes of
# three programs.
#
# The Tumor cells carry three programs with known density-expression
# correlations: proliferation 0.7, hypoxia -0.7 and epithelial 0. Each program
# count is drawn by denmark_simulate() and split among its genes with fixed
# shares, so the program has the correlation of the draw. HYP0, listed in the
# hypoxia set, follows a pattern of its own (correlation 0.5) and carries more
# than half of the counts of the set. PROL14 has too few counts to be used.
# The capture depth of a cell falls where cells are crowded and varies from
# cell to cell, as in sequencing data; it thins all counts of the cell,
# including its total.

stopifnot(file.exists("DESCRIPTION"), file.exists("data-raw/make_program_data.R"))
pkgload::load_all(".", quiet = TRUE)

rho <- c(proliferation = 0.7, hypoxia = -0.7, epithelial = 0, HYP0 = 0.5)
sim <- denmark_simulate(rho = rho, width = 6, height = 4, n_cells = 5000,
                        sigma2 = 0.8, beta1 = log(40), seed = 11,
                        tissue = "irregular", n_other = 4000, celltype = "Tumor")
set.seed(12)
n <- nrow(sim)
is_tumor <- sim$celltype == "Tumor"

# capture depth in (0, 1]: lower where cells are crowded (the 30 nearest
# cells of any type), with cell-to-cell variation
knn <- function(x, y, k) vapply(seq_along(x), function(i)
  sort((x - x[i])^2 + (y - y[i])^2, partial = k + 1)[k + 1], numeric(1))
density <- -log(knn(sim$x, sim$y, 30))
density <- (density - mean(density)) / stats::sd(density)
depth <- 0.8 * stats::plogis(0.4 - 0.5 * density) * exp(stats::rnorm(n, 0, 0.2))
depth <- pmin(depth, 1)

# genes of each program: share of the program count, in decreasing order
shares <- list(
  proliferation = setNames(1 / (2 + 1:13)^0.8, paste0("PROL", 1:13)),
  hypoxia = setNames(1 / (1 + 1:8)^0.7, paste0("HYP", 1:8)),
  epithelial = setNames(c(0.42, 0.27, 0.19, 0.12), paste0("EPI", 1:4)))
level <- c(proliferation = 0.06, hypoxia = 0.05, epithelial = 0.06)

genes <- list()
for (p in names(shares)) {
  thinned <- stats::rbinom(n, sim[[p]], level[[p]] * depth)
  pr <- shares[[p]] / sum(shares[[p]])
  counts <- vapply(seq_len(n), function(i) stats::rmultinom(1, thinned[i], pr)[, 1],
                   numeric(length(pr)))
  genes[[p]] <- as.data.frame(t(counts))
  names(genes[[p]]) <- names(pr)
}
genes$extra <- data.frame(PROL14 = stats::rbinom(n, 1, ifelse(is_tumor, 0.006 * depth, 0)),
                          HYP0 = stats::rbinom(n, sim$HYP0, 0.04 * depth))
total <- stats::rpois(n, 3000 * depth)

cells <- cbind(sim[, c("x", "y", "celltype")], total = total,
               genes$proliferation, genes$extra["PROL14"], genes$extra["HYP0"],
               genes$hypoxia, genes$epithelial)
cells$x <- round(cells$x, 6)
cells$y <- round(cells$y, 6)
rownames(cells) <- NULL

con <- gzfile("inst/extdata/denmark_programs.csv.gz", "wb")   # binary: LF line endings
utils::write.csv(cells, con, row.names = FALSE)
close(con)
