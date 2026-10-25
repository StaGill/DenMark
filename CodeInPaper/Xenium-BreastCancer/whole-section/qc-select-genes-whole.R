## ---- CONFIGURE ME ------------------------------------------------------------
## Run once on the cluster before the array jobs; the shell may export these
## three variables, as cluster-long-whole-DCIS.sh does, to override the defaults.
DATA_DIR   <- Sys.getenv("DATA_DIR", "data")       # holds Dataset/xenium-breast-cancer/rep1/
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", "output")   # the QC gene lists go to qc-whole/ below it
REPO_ROOT  <- Sys.getenv("REPO_ROOT", ".")         # root of the DenMark checkout (not used here)
## -----------------------------------------------------------------------------

######################################################
# Per-group gene QC on WHOLE TISSUE (divisions B & C)
# ----------------------------------------------------
# Run ONCE before the array jobs. For each cell-type
# group, keep genes passing BOTH:
#   (1) across-cell variance in the TOP 95%
#   (2) non-zero expression in >= 5% of that group's cells
#
# Groups (whole tissue, no ROI box) - the two NEW merged groups:
#   DCIS  = DCIS_1 + DCIS_2                     (division B)
#   Tumor = DCIS_1 + DCIS_2 + Invasive_Tumor    (division C)
# (Immune / Invasive_Tumor are reused from figure A - not run here.)
######################################################

data_dir <- file.path(DATA_DIR, 'Dataset/xenium-breast-cancer/rep1/')
qc_dir   <- file.path(OUTPUT_DIR, 'qc-whole/')
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

cell_gene    <- read.csv(paste0(data_dir, 'cell_by_gene_expression.csv'))
cell_spatial <- read.csv(paste0(data_dir, 'cell_spatial.csv'))
cell_spatial$x <- cell_spatial$x / 1000
cell_spatial$y <- cell_spatial$y / 1000

# --- the two NEW merged groups -----------------------
groups <- list(
  DCIS  = c("DCIS_1", "DCIS_2"),
  Tumor = c("DCIS_1", "DCIS_2", "Invasive_Tumor")
)

# --- gene columns & control exclusion ----------------
gene_cols <- 2:ncol(cell_gene)
G <- length(gene_cols)
ctrl_pat <- "^(BLANK_|NegControlCodeword_|NegControlProbe_|antisense_|DeprecatedCodeword_|UnassignedCodeword_)"
is_gene  <- !grepl(ctrl_pat, colnames(cell_gene)[gene_cols])
n_real   <- sum(is_gene)
cat(sprintf("Columns: %d total | %d real genes | %d controls excluded\n",
            G, n_real, G - n_real))

expected_n_genes <- 313
if (n_real != expected_n_genes)
  stop(sprintf("Expected %d real genes but found %d - wrong file / offset?",
               expected_n_genes, n_real), call. = FALSE)

# --- QC thresholds -----------------------------------
var_top_frac <- 0.95
det_min      <- 0.05

# --- run QC per group --------------------------------
for (nm in names(groups)) {
  members <- groups[[nm]]
  idx <- which(cell_spatial$celltype %in% members)     # whole tissue

  expr <- as.matrix(cell_gene[idx, gene_cols, drop = FALSE])

  gene_var <- apply(expr, 2, var)
  var_thr  <- quantile(gene_var, probs = 1 - var_top_frac, na.rm = TRUE)
  pass_var <- gene_var >= var_thr

  detect   <- colMeans(expr > 0)
  pass_det <- detect >= det_min

  keep <- pass_var & pass_det & is_gene
  keep[is.na(keep)] <- FALSE

  gene_ids   <- which(keep)
  gene_names <- colnames(cell_gene)[gene_cols][keep]

  saveRDS(list(gene_ids = gene_ids, gene_names = gene_names,
               n_pass = length(gene_ids), group = nm, members = members),
          file = paste0(qc_dir, nm, '_qc_genes.rds'))
  write.csv(data.frame(gene_id = gene_ids, gene = gene_names),
            paste0(qc_dir, nm, '_qc_genes.csv'), row.names = FALSE)
  writeLines(as.character(length(gene_ids)),
             paste0(qc_dir, nm, '_ngenes.txt'))

  cat(sprintf("%s: %d cells | %d / %d real genes pass QC (var>=%.4g, detect>=%.0f%%)\n",
              nm, length(idx), length(gene_ids), n_real, var_thr, 100 * det_min))
}
