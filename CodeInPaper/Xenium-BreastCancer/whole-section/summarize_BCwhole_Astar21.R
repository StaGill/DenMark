## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # not used: the inputs are the per-gene fits of DenMark-BC-RA-cluster-whole.R
OUTPUT_DIR <- "output"   # holds Revision-R1/Xenium-BC-celltype-divisions/cluster-results/ (see ../README.md)
REPO_ROOT  <- "."        # root of the DenMark checkout
## -----------------------------------------------------------------------------

###############################################################################
##  Summarize single-cell (BC-whole) Astar[2,1] -> ONE csv per cell type.
##  Two folders (BC-whole-DCIS, BC-whole-tumor) -> two csvs, each with:
##    gene_name, posterior mean, and significance (95% CI of Astar[2,1] vs 0).
##
##  Per-gene files look like <GENE>_M2_<n>.csv (e.g. ABCC11_M2_1.csv); the gene
##  name comes from the FILENAME. Each file has an Astar[2,1] row with mean,
##  2.5% and 97.5% columns (first unnamed column is a row index, ignored).
###############################################################################

## ---- EDIT: base dir that CONTAINS the two folders --------------------------
BASE_DIR <- file.path(OUTPUT_DIR, "BC-whole")   # where whole-section/DenMark-BC-RA-cluster-whole.R writes

## folder -> output csv name  (EDIT folder names if they differ on disk)
jobs <- c(
  "BC-whole-DCIS"  = "BC-whole-DCIS_Astar21_summary.csv",
  "BC-whole-Tumor" = "BC-whole-tumor_Astar21_summary.csv"
)

## filename pattern: <GENE>_M2_<number>.csv  -> strip suffix to get the gene
GENE_SUFFIX <- "_M2_[0-9]+\\.csv$"

###############################################################################
gene_from_file <- function(fname) sub(GENE_SUFFIX, "", fname)

## pull the Astar[2,1] row from one file -> gene, mean, CI, significance
grab_A21 <- function(f) {
  d <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d)) { warning("unreadable: ", f); return(NULL) }
  i <- which(d$variable == "Astar[2,1]")
  if (length(i) != 1) { warning("no single Astar[2,1] row: ", f); return(NULL) }
  lo <- as.numeric(d[["2.5%"]][i]); hi <- as.numeric(d[["97.5%"]][i])
  data.frame(
    gene_name   = gene_from_file(basename(f)),
    A21_mean    = as.numeric(d$mean[i]),          # posterior mean
    A21_q2.5    = lo,                             # 95% CI lower
    A21_q97.5   = hi,                             # 95% CI upper
    significant = (lo > 0) | (hi < 0),           # CI excludes 0
    direction   = ifelse(lo > 0, "positive", ifelse(hi < 0, "negative", "n.s.")),
    stringsAsFactors = FALSE)
}

## ---- process each folder, write its own csv --------------------------------
for (fo in names(jobs)) {
  dirp <- file.path(BASE_DIR, fo)
  if (!dir.exists(dirp)) { warning("folder not found: ", dirp); next }

  files <- list.files(dirp, pattern = GENE_SUFFIX, full.names = TRUE)
  res   <- do.call(rbind, lapply(files, grab_A21))
  if (is.null(res) || nrow(res) == 0) { warning("no gene files matched in ", dirp); next }

  res <- res[order(-res$A21_mean), ]            # strongest positive coupling first
  out <- file.path(BASE_DIR, jobs[[fo]])
  write.csv(res, out, row.names = FALSE)

  cat(sprintf("%-16s: %d genes, %d significant (%d pos / %d neg) -> %s\n",
              fo, nrow(res), sum(res$significant),
              sum(res$direction == "positive"), sum(res$direction == "negative"),
              basename(out)))
}
cat("\nDone.\n")
###############################################################################
