# DenMark gene scan on a SLURM cluster, one gene per array task.
#
#   Rscript denmark-scan.R prepare     once: tissue boundary, gene list, compiled model
#   sbatch --array=1-$(cat denmark_scan/n_genes.txt) denmark-scan-array.sh
#   Rscript denmark-scan.R collect     after the array: denmark_scan/denmark_scan.csv
#
# Each task fits one gene and saves its result in SAVE_DIR. A task whose gene is
# already saved ends at once, so an interrupted array can simply be resubmitted.
# "collect" runs denmark_scan() with the same SAVE_DIR, which loads the saved
# results instead of refitting and returns the usual table.

library(DenMark)

# ---- CONFIGURE ME -------------------------------------------------------------
CELLS_CSV <- "cells.csv"      # one row per cell: x, y, celltype and raw counts per gene
CELLTYPE  <- "Tumor"          # the population analyzed
SAVE_DIR  <- "denmark_scan"   # shared by all tasks
SETTINGS  <- list(grid_size = 30, basis = 25, boundfactor = 1.5,
                  chains = 4, parallel_chains = 4, seed = 2024)
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "prepare"
input <- file.path(SAVE_DIR, "scan_input.rds")
saved <- function(genes) file.exists(file.path(SAVE_DIR, paste0("denmark_", make.names(genes), ".rds")))

if (mode == "prepare") {
  dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
  cells <- read.csv(CELLS_CSV)
  tissue <- denmark_alpha_shape(cells$x, cells$y)      # tissue boundary from all cells
  focal <- cells[cells$celltype == CELLTYPE, setdiff(names(cells), "celltype")]
  genes <- denmark_filter_genes(focal[, setdiff(names(focal), c("x", "y"))])
  saveRDS(list(focal = focal, genes = genes, tissue = tissue), input)
  writeLines(as.character(length(genes)), file.path(SAVE_DIR, "n_genes.txt"))
  denmark_model(vararea = TRUE)                        # compile once, before the tasks start
  message(length(genes), " genes: submit with --array=1-", length(genes))
} else if (mode == "gene") {
  s <- readRDS(input)
  task <- as.integer(args[2])
  if (task > length(s$genes)) quit(save = "no")
  g <- s$genes[task]
  invisible(do.call(denmark_scan, c(list(s$focal[, c("x", "y", g)], genes = g,
                                         boundary = s$tissue, save_dir = SAVE_DIR), SETTINGS)))
} else if (mode == "collect") {
  s <- readRDS(input)
  done <- saved(s$genes)
  if (any(!done))
    message(sum(!done), " gene(s) have no saved result (failed or not run): ",
            paste(utils::head(s$genes[!done], 10), collapse = ", "),
            if (sum(!done) > 10) ", ..." else "", ". Resubmit the array to fit them.")
  scan <- do.call(denmark_scan, c(list(s$focal, genes = s$genes[done], boundary = s$tissue,
                                       save_dir = SAVE_DIR), SETTINGS))
  utils::write.csv(scan, file.path(SAVE_DIR, "denmark_scan.csv"), row.names = FALSE)
  print(table(scan$class, useNA = "ifany"))
} else {
  stop("The first argument must be prepare, gene or collect.")
}
