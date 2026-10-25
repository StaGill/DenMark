## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds Revision-R1/Reanalysis-MERFISH/ (annotation_out/ and the cluster summaries)
OUTPUT_DIR <- "output"   # fits, tables and figures go to Revision-R1/Reanalysis-MERFISH/ below it
REPO_ROOT  <- "."        # root of the DenMark checkout (Stan models and the grid files)
## -----------------------------------------------------------------------------

###############################################################################
# S3R1 — ASTROCYTE DenMark analysis on 5 selected genes
#   Genes: Aqp4, Sox9, Fgfr3, Ddr2, Celsr2
#   Model: DenMark-vararea.stan (UNCHANGED)  +  nested a21 = 0 variant
#
# Produces, per gene:
#   (1) cell-level expression map; grid-level cell-count & expression maps
#   (a) posterior summary of rho = Astar[2,1]
#   (b) latent-field decomposition, 2 x 4 spatial maps
#   (c) spatial correlation-decay curve for the EXPRESSION field + ranges r1, r2
#   (d) WAIC comparison: full model vs nested (a21 = 0)
#
# PARAMETERISATION (read off the Stan file):
#   A = diag_pre_multiply(sigma, Astar), Astar = cholesky_factor_corr[2]
#     => a11 = sigma[1]
#        a21 = sigma[2] * Astar[2,1]        (Astar[2,1] = rho)
#        a22 = sigma[2] * Astar[2,2] = sigma[2] * sqrt(1 - rho^2)
#   log(lambda1) = log_area + beta0 + a11*w1
#   log(lambda2) = beta1 + log(Y + 0.05) + a21*w1 + a22*w2
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(fields); library(cmdstanr); library(posterior); library(loo)
})

set.seed(2026)

## ---------------------------------------------------------------------------
## 0. Paths & settings  (EDIT)
## ---------------------------------------------------------------------------
BASE      <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/")
STAN_FULL <- file.path(REPO_ROOT, "inst", "stan", "DenMark-vararea.stan")
STAN_NEST <- file.path(REPO_ROOT, "CodeInPaper", "MERFISH-MouseBrain", "DenMark-M2-vararea-a21zero.stan")   # nested model (a21 = 0)
GRID_CSV  <- file.path(REPO_ROOT, "CodeInPaper", "MERFISH-MouseBrain", "grid_area", "grid_area_alpha_S3R1.csv")
CBG_CSV   <- file.path(BASE, "tissue-boundary-9-slices", "annotation_out",
                       "cell_by_gene_with_celltype_S3R1.csv")
META_CSV  <- file.path(BASE, "tissue-boundary-9-slices", "annotation_out",
                       "cell_metadata_with_celltype_S3R1.csv")
OUT_DIR   <- file.path(OUTPUT_DIR, "Revision-R1/Reanalysis-MERFISH", "reanalysis-S3R1/S3R1-5genes")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

TAG    <- "S3R1"
GENES  <- c("Aqp4", "Sox9", "Fgfr3", "Ddr2", "Celsr2")

basis_choice    <- 25
boundary_choice <- 1.5
coord_scale     <- 1000            # um -> mm
area_scale      <- coord_scale^2
ITER_WARMUP     <- 1000
ITER_SAMPLING   <- 1000

SCRATCH <- file.path(tempdir(), "denmark_S3R1")
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
##  theme + palettes
## ---------------------------------------------------------------------------
theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 1, colour = "grey30"),
      axis.line     = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks    = element_line(linewidth = 0.3, colour = "black"),
      axis.text     = element_text(colour = "black"),
      legend.key.size = unit(0.35, "cm"),
      strip.background = element_blank(),
      strip.text    = element_text(face = "bold", size = base_size)
    )
}
# sequential for counts/intensities; diverging for signed latent fields
scale_seq <- function(name) scale_fill_viridis_c(option = "D", name = name)
scale_div <- function(name) scale_fill_gradient2(low = "#2166AC", mid = "white",
                                                 high = "#B2182B", midpoint = 0, name = name)
okabe_ito <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00")   # 5 genes




## ---------------------------------------------------------------------------
## 0. Slice annotation, adding boundary of slice 
## ---------------------------------------------------------------------------
library(sf)
library(ggplot2)
library(patchwork)

## ---------------------------------------------------------------------------
## S3R1 — cell-type annotation map WITH the tissue boundary
## Assumes the annotation has already been run, i.e. you have:
##   results[["S3R1"]]$meta   (or the saved cell_metadata_with_celltype_S3R1.csv)
## and the grid file grid_area_alpha_S3R1.csv for the boundary.
## ---------------------------------------------------------------------------
TAG        <- "S3R1"
GRID_CSV   <- file.path(REPO_ROOT, "CodeInPaper", "MERFISH-MouseBrain", "grid_area",
                        sprintf("grid_area_alpha_%s.csv", TAG))
META_CSV   <- file.path(BASE, "tissue-boundary-9-slices", "annotation_out",
                        sprintf("cell_metadata_with_celltype_%s.csv", TAG))
coord_scale <- 1000     # um -> mm

## ---- 1. Cell-type table ----------------------------------------------------
# either reuse the in-memory object ...
m <- if (exists("results") && !is.null(results[[TAG]])) results[[TAG]]$meta else {
  # ... or read the saved file (cell IDs kept as character)
  h  <- read.csv(META_CSV, nrows = 1, check.names = FALSE)
  cc <- rep(NA, ncol(h)); cc[1] <- "character"
  read.csv(META_CSV, colClasses = cc, check.names = FALSE)
}
m$cell_type <- factor(m$cell_type,
                      levels = c("Neuron", "Astrocyte", "Microglia", "Other"))
m$x <- m$center_x / coord_scale
m$y <- m$center_y / coord_scale

## ---- 2. Tissue boundary from the grid (exterior ring only) -----------------
grid <- read.csv(GRID_CSV, stringsAsFactors = FALSE)
grid$type <- tolower(trimws(grid$type))
grid <- grid[grid$type %in% c("interior", "boundary"), , drop = FALSE]

# rebuild the lattice rectangles (mm)
ux <- tapply(grid$centroid_x, grid$i, function(v) v[1])
uy <- tapply(grid$centroid_y, grid$j, function(v) v[1])
ii <- as.integer(names(ux)); xi <- as.numeric(ux)
jj <- as.integer(names(uy)); yj <- as.numeric(uy)
step_x <- median(diff(xi[order(ii)])); step_y <- median(diff(yj[order(jj)]))

grid$xmin <- (grid$centroid_x - step_x/2) / coord_scale
grid$xmax <- (grid$centroid_x + step_x/2) / coord_scale
grid$ymin <- (grid$centroid_y - step_y/2) / coord_scale
grid$ymax <- (grid$centroid_y + step_y/2) / coord_scale

make_rect <- function(x1, x2, y1, y2)
  sf::st_polygon(list(cbind(c(x1, x2, x2, x1, x1), c(y1, y1, y2, y2, y1))))

polys  <- sf::st_sfc(mapply(make_rect, grid$xmin, grid$xmax, grid$ymin, grid$ymax,
                            SIMPLIFY = FALSE))
polys  <- sf::st_set_precision(polys, 1e6)          # snap so interior edges dissolve
tissue <- sf::st_buffer(sf::st_union(polys), 0)

outer_only <- function(geom) {                       # drop interior rings
  ps <- sf::st_cast(sf::st_geometry(geom), "POLYGON")
  sf::st_sfc(lapply(ps, function(p) sf::st_polygon(list(p[[1]]))))
}
bnd_df <- as.data.frame(sf::st_coordinates(sf::st_boundary(outer_only(tissue))))
grp    <- intersect(c("L1", "L2", "L3"), names(bnd_df))
bnd_df$grp <- interaction(bnd_df[grp], drop = TRUE)

## ---- 3. Plot ---------------------------------------------------------------
pal <- c(Neuron    = "#FDB863",   # warm  (magma-family)
         Astrocyte = "#B2182B",   # deep red — the focal cell type
         Microglia = "#2166AC",   # blue
         Other     = "grey85")

theme_panel <- function()
  theme_classic(base_size = 12) +
  theme(
    axis.line         = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.text         = element_text(size = 15, colour = "black"),
    axis.title        = element_text(size = 20, colour = "black"),
    legend.text       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.key.size   = unit(0.6, "cm"),
    plot.title        = element_text(size = 17, face = "bold", hjust = 0.5),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 14)
  )

# (a) all cell types on one panel, boundary attached
p_all <- ggplot() +
  geom_point(data = m, aes(x, y, colour = cell_type), size = 0.2, alpha = 0.8) +
  geom_path(data = bnd_df, aes(X, Y, group = grp),
            colour = "grey20", linewidth = 0.4) +
  scale_colour_manual(values = pal, name = "Cell Type") +
  coord_fixed() +
  labs(title = '', x = "x (mm)", y = "y (mm)") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_panel()

# (b) one panel per cell type (each type highlighted against the others in grey)
p_facet <- ggplot() +
  geom_point(data = dplyr::select(m, x, y), aes(x, y),
             colour = "grey92", size = 0.2) +                    # all cells, faint
  geom_point(data = m, aes(x, y, colour = cell_type), size = 0.2) +
  geom_path(data = bnd_df, aes(X, Y, group = grp),
            colour = "grey20", linewidth = 0.35) +
  scale_colour_manual(values = pal, guide = "none") +
  facet_wrap(~ cell_type, nrow = 1) +
  coord_fixed() +
  labs(title = "", x = "x (mm)", y = NULL) +
  theme_panel()

p_ct <- p_all + p_facet + patchwork::plot_layout(widths = c(1.3, 4))

ggsave(file.path(FIG_DIR, sprintf("00_celltype_map_%s.png", TAG)), p_ct,
       width = 18, height = 4.2, dpi = 300)

p_ct


## ---- Only keep astrocytes  --------------------------------------------------


p_all


# (b) astrocytes only, against a faint grey backdrop of all cells
p_astro <- ggplot() +
  geom_point(data = dplyr::select(m, x, y), aes(x, y),
             colour = "grey92", size = 0.2) +                        # all cells, faint
  geom_point(data = subset(m, cell_type == "Astrocyte"),
             aes(x, y), colour = pal["Astrocyte"], size = 0.2) +     # astrocytes highlighted
  geom_path(data = bnd_df, aes(X, Y, group = grp),
            colour = "grey20", linewidth = 0.35) +
  coord_fixed() +
  labs(title = " ", x = "x (mm)", y = NULL) +
  theme_panel()

p_astro






## ---- 4. Composition table --------------------------------------------------
tab <- table(m$cell_type)
print(data.frame(cell_type = names(tab),
                 n = as.integer(tab),
                 pct = round(100 * as.numeric(tab) / sum(tab), 1)))



library(ggplot2)
library(patchwork)

## ---------------------------------------------------------------------------
## S3R1 — cell-type composition figure
## Uses `m` (annotated metadata), `pal`, and `theme_panel()` from the map script.
## ---------------------------------------------------------------------------
comp <- as.data.frame(table(m$cell_type), stringsAsFactors = FALSE)
names(comp) <- c("cell_type", "n")
comp$cell_type <- factor(comp$cell_type,
                         levels = c("Neuron", "Astrocyte", "Microglia", "Other"))
comp$pct   <- 100 * comp$n / sum(comp$n)
comp$label <- sprintf("%s\n(%.1f%%)", format(comp$n, big.mark = ","), comp$pct)

## ---- (a) horizontal bar chart: counts, with n and % labelled ---------------
p_bar <- ggplot(comp, aes(x = n, y = factor(cell_type, levels = rev(levels(cell_type))),
                          fill = cell_type)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = label), hjust = -0.08, size = 8, lineheight = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.22))) +
  labs(title = 'Cell-type Composition',
       x = "Number of cells", y = NULL) +
  theme_panel() +
  theme(axis.text.x  = element_blank(),           
        axis.ticks.x = element_line())+
  theme(
    axis.text.x  = element_blank(),           
    axis.ticks.x = element_line(),
    #axis.text.x  = element_text(size = 25, colour = "black"),
    axis.text.y  = element_text(size = 25, colour = "black"),
    axis.title.x = element_text(size = 24, colour = "black"),
    plot.title   = element_text(size = 24, face = "bold", hjust = 0.5)
  )
 
p_bar

ggsave(file.path(FIG_DIR, sprintf("00b_celltype_composition_%s.png", TAG)),
       p_bar, width = 9, height = 5.5, dpi = 300)



## --------------------------------------------------------------------------
## Add genome-wide summary -----------
## --------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## 0. Config  (EDIT PATHS)
## ---------------------------------------------------------------------------
BASE        <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH")
RESULTS_DIR <- file.path(BASE, "cluster-for-8-rest-slices", "results-S3R1")
CBG_CSV     <- file.path(BASE, "tissue-boundary-9-slices", "annotation_out",
                         "cell_by_gene_with_celltype_S3R1.csv")
META_CSV    <- file.path(BASE, "tissue-boundary-9-slices", "annotation_out",
                         "cell_metadata_with_celltype_S3R1.csv")
GRID_CSV    <- file.path(REPO_ROOT, "CodeInPaper", "MERFISH-MouseBrain", "grid_area", "grid_area_alpha_S3R1.csv")
OUT_DIR     <- file.path(OUTPUT_DIR, "Revision-R1/Reanalysis-MERFISH", "cluster-for-8-rest-slices", "summary-S3R1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FOCAL_GENES <- c("Aqp4", "Sox9", "Fgfr3", "Ddr2", "Celsr2")

# --- QC tuning knobs --------------------------------------------------------
VAR_ON_LOG   <- FALSE   # QC1: FALSE = raw-count variance (as specified);
#      TRUE  = variance of log1p(counts) (less abundance-driven)
VAR_KEEP_TOP <- 0.95    # keep top 95% by variance -> drop the lowest 5%
DETECT_UNIT  <- "grid"  # QC2 spatial unit: "grid" (model units) or "cell"
DETECT_MIN   <- 0.05    # require non-zero in >= 5% of those locations

# --- convergence flag (reported, NOT one of your two QC gates) --------------
RHAT_MAX     <- 1.01
DIV_MAX      <- 0       # set higher (e.g. 10) if you want to tolerate a few

## ===========================================================================
## 1. Read the per-gene DenMark summaries  ->  density-expression correlation
## ===========================================================================
csv_files <- list.files(RESULTS_DIR, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(csv_files) > 0)

read_one <- function(f) {
  d <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)                 # check.names=FALSE keeps "2.5%"/"97.5%"
  if (is.null(d) || !("variable" %in% names(d))) return(NULL)
  a <- d[d$variable == "Astar[2,1]", , drop = FALSE]      # the density-expression correlation
  if (nrow(a) != 1) return(NULL)
  data.frame(
    gene_name      = a$gene_name,
    gene_id        = a$gene_id,
    slice          = a$slice,
    dens_expr_corr = a$mean,                              # <- Astar[2,1] posterior mean
    corr_sd        = a$sd,
    corr_lo        = a[["2.5%"]],
    corr_hi        = a[["97.5%"]],
    corr_rhat      = a$rhat,
    corr_ess_bulk  = a$ess_bulk,
    n_grids        = a$n_grids,
    max_rhat       = a$max_rhat,
    num_divergent  = a$num_divergent,
    waic           = a$waic,
    waic_se        = a$waic_se,
    stringsAsFactors = FALSE
  )
}

res <- dplyr::bind_rows(lapply(csv_files, read_one))
message(sprintf("Read %d gene summaries (of %d CSVs).", nrow(res), length(csv_files)))

# credible sign of the correlation + convergence flag
res <- res %>%
  mutate(
    corr_credible = dplyr::case_when(
      corr_lo > 0 ~ "Positive",
      corr_hi < 0 ~ "Negative",
      TRUE        ~ "Uninformative"
    ),
    conv_ok = max_rhat < RHAT_MAX & num_divergent <= DIV_MAX
  )

genes_all <- res$gene_name                                # the panel we have results for

## ===========================================================================
## 2. QC gates from the raw astrocyte data
## ===========================================================================

## ---- 2a. astrocyte cell-by-gene matrix -------------------------------------
meta <- read.csv(META_CSV, check.names = FALSE, stringsAsFactors = FALSE)
id_meta <- names(meta)[1]
astro_ids <- meta[[id_meta]][meta$cell_type == "Astrocyte"]

cbg <- read.csv(CBG_CSV, check.names = FALSE, stringsAsFactors = FALSE)
id_cbg <- names(cbg)[1]
cbg_a  <- cbg[cbg[[id_cbg]] %in% astro_ids, , drop = FALSE]

genes_use <- intersect(genes_all, names(cbg_a))           # genes present in both
missing   <- setdiff(genes_all, genes_use)
if (length(missing))
  warning(sprintf("%d result genes not found as columns in the cell-by-gene file: %s",
                  length(missing), paste(head(missing, 10), collapse = ", ")))

mat <- as.matrix(cbg_a[, genes_use, drop = FALSE])        # astro cells x genes (counts)
storage.mode(mat) <- "double"

## ---- 2b. QC1: variance across astrocyte cells ------------------------------
mat_for_var <- if (VAR_ON_LOG) log1p(mat) else mat
var_cells   <- apply(mat_for_var, 2, var)
var_cut     <- quantile(var_cells, 1 - VAR_KEEP_TOP, na.rm = TRUE)  # 5th pct of variance
pass_var    <- var_cells > var_cut                        # top 95%

## ---- 2c. QC2: detection fraction (grids by default, cells optional) --------
frac_cells_nonzero <- colMeans(mat > 0)                   # per-cell detection rate

# grid-level detection: bin astrocyte cells onto the model lattice
grid <- read.csv(GRID_CSV, stringsAsFactors = FALSE)
grid$type <- tolower(trimws(grid$type))
grid <- grid[grid$type %in% c("interior", "boundary"), , drop = FALSE]   # model grids only
n_model_grids <- nrow(grid)

# regular lattice: unique centroids + half-step edges
xs <- sort(unique(grid$centroid_x)); ys <- sort(unique(grid$centroid_y))
step_x <- median(diff(xs)); step_y <- median(diff(ys))
ex <- c(xs - step_x/2, xs[length(xs)] + step_x/2)
ey <- c(ys - step_y/2, ys[length(ys)] + step_y/2)

# astrocyte cell coordinates (from metadata), aligned to cbg_a row order
coord <- meta[match(cbg_a[[id_cbg]], meta[[id_meta]]), c("center_x", "center_y")]
ci <- findInterval(coord$center_x, ex)                    # 1..length(xs) if inside
cj <- findInterval(coord$center_y, ey)
inside <- ci >= 1 & ci <= length(xs) & cj >= 1 & cj <= length(ys)

# keep only cells that land in a REAL model grid (centroid pair exists in `grid`)
grid_keys  <- paste(grid$centroid_x, grid$centroid_y, sep = "_")
cell_keys  <- rep(NA_character_, nrow(cbg_a))
cell_keys[inside] <- paste(xs[ci[inside]], ys[cj[inside]], sep = "_")
valid <- !is.na(cell_keys) & cell_keys %in% grid_keys

gsum <- rowsum(mat[valid, , drop = FALSE], group = cell_keys[valid])   # occupied grids x genes
frac_grids_nonzero <- colSums(gsum > 0) / n_model_grids   # denom = ALL model grids

frac_detect <- if (DETECT_UNIT == "grid") frac_grids_nonzero else frac_cells_nonzero
pass_detect <- frac_detect >= DETECT_MIN

## ---- 2d. assemble the QC table ---------------------------------------------
qc <- data.frame(
  gene_name          = genes_use,
  var_cells          = as.numeric(var_cells),
  pass_var           = as.logical(pass_var),
  frac_cells_nonzero = as.numeric(frac_cells_nonzero),
  frac_grids_nonzero = as.numeric(frac_grids_nonzero),
  pass_detect        = as.logical(pass_detect),
  stringsAsFactors   = FALSE
) %>%
  mutate(qc_pass = pass_var & pass_detect)

## ===========================================================================
## 3. Merge results + QC, save, and report
## ===========================================================================
final <- res %>%
  left_join(qc, by = "gene_name") %>%
  mutate(is_focal = gene_name %in% FOCAL_GENES) %>%
  arrange(desc(dens_expr_corr))

write.csv(final, file.path(OUT_DIR, "panelwide_denmark_S3R1.csv"), row.names = FALSE)
saveRDS(final,   file.path(OUT_DIR, "panelwide_denmark_S3R1.rds"))

## ---- console report --------------------------------------------------------
kept <- final %>% filter(qc_pass)
cat("\n================ QC FILTERING REPORT (S3R1) ================\n")
cat(sprintf("Genes with fitted results        : %d\n", nrow(final)))
cat(sprintf("QC1 pass (variance top %.0f%%)     : %d\n", 100*VAR_KEEP_TOP, sum(final$pass_var, na.rm = TRUE)))
cat(sprintf("QC2 pass (detect >= %.0f%% of %ss) : %d\n", 100*DETECT_MIN, DETECT_UNIT, sum(final$pass_detect, na.rm = TRUE)))
cat(sprintf("Passing BOTH gates               : %d\n", nrow(kept)))
cat(sprintf("  ...also converged (rhat<%.2f,div<=%d): %d\n", RHAT_MAX, DIV_MAX, sum(kept$conv_ok, na.rm = TRUE)))
cat("\n---- Density-expression correlation, QC-passing genes ----\n")
cat(sprintf("median = %.3f | mean = %.3f | range [%.3f, %.3f]\n",
            median(kept$dens_expr_corr), mean(kept$dens_expr_corr),
            min(kept$dens_expr_corr), max(kept$dens_expr_corr)))
cat(sprintf("credibly positive (95%% CI > 0)   : %d\n", sum(kept$corr_credible == "Positive")))
cat(sprintf("credibly negative (95%% CI < 0)   : %d\n", sum(kept$corr_credible == "Negative")))
cat(sprintf("spans zero                       : %d\n", sum(kept$corr_credible == "Uninformative")))
cat("\n---- Focal genes ----\n")
print(final %>% filter(is_focal) %>%
        select(gene_name, dens_expr_corr, corr_lo, corr_hi,
               corr_credible, qc_pass, conv_ok))
cat("============================================================\n")

## ===========================================================================
## 4. Panel-wide summary figure: density-expression correlation, all genes
## ===========================================================================
plot_df <- kept %>%
  arrange(dens_expr_corr) %>%
  mutate(rank = row_number())

sign_cols <- c(Positive = "#B2182B", 
               Negative = "#2166AC", 
               `Uninformative` = "grey70")

p <- ggplot(plot_df, aes(rank, dens_expr_corr)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_linerange(aes(ymin = corr_lo, ymax = corr_hi, colour = corr_credible),
                 linewidth = 0.25, alpha = 0.7) +
  geom_point(aes(colour = corr_credible), size = 0.5) +
  geom_point(data = filter(plot_df, is_focal), colour = "black", size = 1.6) +
  ggrepel::geom_text_repel(
    data = filter(plot_df, is_focal),
    aes(label = gene_name), size = 9, fontface = "italic",
    box.padding = 0.5, min.segment.length = 0, max.overlaps = Inf) +
  scale_colour_manual(values = sign_cols, name = "95% CI") +
  labs(
    title = '',
    x = 'Gene Index',
    y = "Density-expression correlation") +
  theme_classic(base_size = 25) +
  theme(plot.title      = element_text(face = "bold", size = 11),
        axis.text       = element_text(colour = "black"),
        legend.position = c(0.02, 0.98),          # x, y in [0,1]; top-left corner
        legend.justification = c(0, 1),           # anchor the legend's top-left to that point
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        legend.title    = element_text(size = 24),
        legend.text     = element_text(size = 22),
        legend.key.size = unit(0.8, "cm"))

ggsave(file.path(OUT_DIR, "panelwide_denmark_corr_S3R1.png"),
       p, width = 8, height = 4.5, dpi = 300)
p



## ---------------------------------------------------------
## Astrocyte enrichment map (contour map) ------
## ---------------------------------------------------------
read_idchar <- function(p) {
  h  <- read.csv(p, nrows = 1, check.names = FALSE)
  cc <- rep(NA, ncol(h)); cc[1] <- "character"
  read.csv(p, colClasses = cc, check.names = FALSE)
}
cell_by_gene  <- read_idchar(CBG_CSV);  names(cell_by_gene)[1]  <- "cell_id"
cell_location <- read_idchar(META_CSV); names(cell_location)[1] <- "cell_id"

## ---------------------------------------------------------------------------
## Build the objects the bridge figure needs:  cells, grid_summ,
## EX_POS / EX_NEG (example genes) and their correlations
## ---------------------------------------------------------------------------

## ---- 1. pick the example genes: strongest +ve and -ve coupling (QC-passing) -
ex_pos_row <- kept %>% dplyr::slice_max(dens_expr_corr, n = 1)
ex_neg_row <- kept %>% dplyr::slice_min(dens_expr_corr, n = 1)

EX_POS   <- ex_pos_row$gene_name;  corr_pos <- ex_pos_row$dens_expr_corr
EX_NEG   <- ex_neg_row$gene_name;  corr_neg <- ex_neg_row$dens_expr_corr
cat(sprintf("EX_POS = %s (corr = %+.2f)\nEX_NEG = %s (corr = %+.2f)\n",
            EX_POS, corr_pos, EX_NEG, corr_neg))

## ---- 2. `cells`: astrocyte coordinates in um -------------------------------
# `coord` is already aligned to cbg_a row order from the QC step
cells <- data.frame(x = coord$center_x, y = coord$center_y)
cells <- cells[valid, , drop = FALSE]          # keep cells that fall in a model grid
cat(sprintf("cells: %d astrocytes\n", nrow(cells)))

## ---- 3. `grid_summ`: per-grid summed expression for the two example genes --
# `gsum` (from the QC step) is occupied-grids x genes, rownames = "centroidx_centroidy"
key_xy <- do.call(rbind, strsplit(rownames(gsum), "_", fixed = TRUE))
occ <- data.frame(gx  = as.numeric(key_xy[, 1]),
                  gy  = as.numeric(key_xy[, 2]),
                  pos = as.numeric(gsum[, EX_POS]),
                  neg = as.numeric(gsum[, EX_NEG]),
                  row.names = NULL)

# fill in the model grids that contain no astrocytes, so the raster is complete
grid_summ <- merge(data.frame(gx = grid$centroid_x, gy = grid$centroid_y),
                   occ, by = c("gx", "gy"), all.x = TRUE)
grid_summ$pos[is.na(grid_summ$pos)] <- 0
grid_summ$neg[is.na(grid_summ$neg)] <- 0

cat(sprintf("grid_summ: %d grids (%d with astrocytes)\n",
            nrow(grid_summ), nrow(occ)))

## ---- rebuild `cells` (astrocyte-only, coords in um) ------------------------
cbg_a  <- cell_by_gene[cell_by_gene$cell_type  == "Astrocyte", , drop = FALSE]
meta_a <- cell_location[cell_location$cell_type == "Astrocyte", , drop = FALSE]
meta_a <- meta_a[match(cbg_a$cell_id, meta_a$cell_id), , drop = FALSE]
stopifnot(all(cbg_a$cell_id == meta_a$cell_id))

cells <- data.frame(cell_id = meta_a$cell_id,
                    x = meta_a$center_x,
                    y = meta_a$center_y,
                    cbg_a[, GENES, drop = FALSE],
                    check.names = FALSE)
cat(sprintf("cells: %d astrocytes, %d gene columns\n", nrow(cells), length(GENES)))



cells_mm <- transform(cells,     xm = x/1000, ym = y/1000)   # astrocyte coords, mm
grid_mm  <- transform(grid_summ, gx = gx/1000, gy = gy/1000) # gridded values, mm



p_dens <- ggplot(cells_mm, aes(xm, ym)) +
  stat_density_2d(aes(fill = after_stat(ndensity)),   # ndensity = 0..1, more readable than raw
                  geom = "raster", contour = FALSE) +
  geom_density_2d(colour = "white", linewidth = 0.2, bins = 8) +  # optional: contour lines on top
  coord_fixed(expand = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Scaled \nDensity") +
  labs(title = "", x = "x (mm)", y = "y (mm)") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(face = "bold", hjust = .5),
        legend.key.width = unit(0.6, "cm"))+
  guides(fill = guide_colourbar(
    barwidth  = unit(0.8, "cm"),    # wider
    barheight = unit(6,   "cm"),    # longer
    frame.colour = "black",
    frame.linewidth = 0.3,
    ticks.colour = "black"
  ))
p_dens




## ---------------------------------------------------------------------------
## 1. Load data; keep ASTROCYTES only
## ---------------------------------------------------------------------------
# keep the long cell IDs as character (they are ~39-digit integers)
read_idchar <- function(p) {
  h  <- read.csv(p, nrows = 1, check.names = FALSE)
  cc <- rep(NA, ncol(h)); cc[1] <- "character"
  read.csv(p, colClasses = cc, check.names = FALSE)
}
cell_by_gene  <- read_idchar(CBG_CSV)
cell_location <- read_idchar(META_CSV)
names(cell_by_gene)[1]  <- "cell_id"
names(cell_location)[1] <- "cell_id"

cbg_a  <- cell_by_gene [cell_by_gene$cell_type  == "Astrocyte", , drop = FALSE]
meta_a <- cell_location[cell_location$cell_type == "Astrocyte", , drop = FALSE]
meta_a <- meta_a[match(cbg_a$cell_id, meta_a$cell_id), , drop = FALSE]   # align
stopifnot(all(cbg_a$cell_id == meta_a$cell_id))

missing_genes <- setdiff(GENES, names(cbg_a))
if (length(missing_genes))
  stop("Gene(s) not in the panel: ", paste(missing_genes, collapse = ", "))

cells <- data.frame(cell_id = meta_a$cell_id,
                    x = meta_a$center_x, y = meta_a$center_y,
                    cbg_a[, GENES, drop = FALSE],
                    check.names = FALSE)
cat(sprintf("[%s] astrocytes = %d\n", TAG, nrow(cells)))

## ---------------------------------------------------------------------------
## 2. Grid: load, lattice geometry, bin astrocytes
## ---------------------------------------------------------------------------
grid <- read.csv(GRID_CSV, stringsAsFactors = FALSE)
grid$type <- tolower(trimws(grid$type))
grid <- grid[grid$type %in% c("interior", "boundary"), , drop = FALSE]
stopifnot(nrow(grid) > 0, all(c("i", "j") %in% names(grid)))

# lattice step (rectangles for plotting; boundary cells are clipped in AREA only)
ux <- tapply(grid$centroid_x, grid$i, function(v) v[1])
uy <- tapply(grid$centroid_y, grid$j, function(v) v[1])
ii <- as.integer(names(ux)); xi <- as.numeric(ux)
jj <- as.integer(names(uy)); yj <- as.numeric(uy)
step_x <- median(diff(xi[order(ii)])); step_y <- median(diff(yj[order(jj)]))
x0 <- xi[which.min(ii)] - min(ii) * step_x
y0 <- yj[which.min(jj)] - min(jj) * step_y

# assign each astrocyte to a grid cell
ci <- round((cells$x - x0) / step_x); cj <- round((cells$y - y0) / step_y)
inside <- ci >= min(ii) & ci <= max(ii) & cj >= min(jj) & cj <= max(jj)
gkey <- grid$i * 100000L + grid$j
ckey <- rep(NA_real_, nrow(cells)); ckey[inside] <- ci[inside] * 100000L + cj[inside]
cells$grid_row <- match(ckey, gkey)
valid <- !is.na(cells$grid_row)
cat(sprintf("binned %d/%d astrocytes into %d tissue grids (%d interior, %d boundary)\n",
            sum(valid), nrow(cells), nrow(grid),
            sum(grid$type == "interior"), sum(grid$type == "boundary")))

# Y = astrocyte count per grid; M_<gene> = summed expression per grid
grid$Y <- tabulate(cells$grid_row[valid], nbins = nrow(grid))
for (g in GENES) {
  mk <- tapply(cells[[g]][valid], cells$grid_row[valid], sum)
  v  <- rep(0, nrow(grid)); v[as.integer(names(mk))] <- as.numeric(mk)
  grid[[paste0("M_", g)]] <- v
}

# per-grid exposure: interior = full cell area; boundary = clipped area (from CSV)
full_area     <- median(grid$area[grid$type == "interior"])
grid$area_use <- ifelse(grid$type == "interior", full_area, grid$area)
stopifnot(all(grid$area_use > 0))

# plotting rectangles (mm)
grid$xmin <- (grid$centroid_x - step_x/2) / coord_scale
grid$xmax <- (grid$centroid_x + step_x/2) / coord_scale
grid$ymin <- (grid$centroid_y - step_y/2) / coord_scale
grid$ymax <- (grid$centroid_y + step_y/2) / coord_scale
grid$gx_mm <- grid$centroid_x / coord_scale
grid$gy_mm <- grid$centroid_y / coord_scale

# helper: draw a grid-level field, with boundary grids outlined (context)
plot_grid_field <- function(values, title, diverging = FALSE, legend = "") {
  d <- grid; d$val <- values
  p <- ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = val)) +
    # boundary grids: dashed outline so clipped-area context is visible
    geom_rect(data = subset(d, type == "boundary"),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, colour = "grey25", linewidth = 0.15, linetype = "22") +
    coord_fixed() +
    labs(title = title, x = "x (mm)", y = "y (mm)") +
    theme_nature()
  p + if (diverging) scale_div(legend) else scale_seq(legend)
}
  
## ---------------------------------------------------------------------------
## (1) VISUALISATION — cell-level expression, magma, tissue outline only
##   - colour = log(expression + 0.05)   (matches the model's log treatment)
##   - zero-count cells drawn separately in grey  -> palette spans only the
##     expressing cells, which is what gives the contrast
##   - scale stretched over the 2%-98% range of non-zero values (outliers squished)
##   - single tissue boundary: grid rectangles unioned with coordinate snapping
##     so interior edges dissolve; only the exterior ring is drawn
## ---------------------------------------------------------------------------
library(sf)
library(scales)

## ---- 1. Tissue outline (exterior ring only) --------------------------------
make_rect <- function(x1, x2, y1, y2)
  sf::st_polygon(list(cbind(c(x1, x2, x2, x1, x1),
                            c(y1, y1, y2, y2, y1))))

polys <- sf::st_sfc(mapply(make_rect, grid$xmin, grid$xmax, grid$ymin, grid$ymax,
                           SIMPLIFY = FALSE))
polys  <- sf::st_set_precision(polys, 1e6)   # snap so shared edges are identical
tissue <- sf::st_union(polys)
tissue <- sf::st_buffer(tissue, 0)           # heal hairline seams

# keep ONLY the exterior ring of each polygon -> no lines inside the tissue
outer_only <- function(geom) {
  ps <- sf::st_cast(sf::st_geometry(geom), "POLYGON")
  sf::st_sfc(lapply(ps, function(p) sf::st_polygon(list(p[[1]]))))
}
tissue_outer <- outer_only(tissue)

bnd_df <- as.data.frame(sf::st_coordinates(sf::st_boundary(tissue_outer)))
grp    <- intersect(c("L1", "L2", "L3"), names(bnd_df))
bnd_df$grp <- interaction(bnd_df[grp], drop = TRUE)

## ---- 2. Cell data: log transform, split zeros from expressing cells --------
cells_long <- cells %>%
  dplyr::filter(!is.na(grid_row)) %>%
  dplyr::select(x, y, dplyr::all_of(GENES)) %>%
  tidyr::pivot_longer(dplyr::all_of(GENES), names_to = "gene", values_to = "expr") %>%
  dplyr::mutate(gene     = factor(gene, levels = GENES),
                log_expr = log(expr + 0.05),
                x = x / coord_scale,
                y = y / coord_scale)

cells_zero <- dplyr::filter(cells_long, expr == 0)
cells_pos  <- dplyr::filter(cells_long, expr >  0) %>%
  dplyr::arrange(log_expr)                  # bright cells drawn last (on top)

# stretch the palette across the middle 96% of the NON-ZERO range
qlim <- stats::quantile(cells_pos$log_expr, c(0.02, 0.98), na.rm = TRUE)

## ---- 3. Plot ---------------------------------------------------------------
p_cells <- ggplot() +
  geom_path(data = bnd_df, aes(X, Y, group = grp),
            colour = "grey30", linewidth = 0.3) +
  geom_point(data = cells_zero, aes(x, y),
             colour = "grey92", size = 0.6, alpha = 0.6) +
  geom_point(data = cells_pos, aes(x, y, colour = log_expr),
             size = 0.6) +
  scale_colour_viridis_c(option = "A", direction = 1,      # magma; dark = high
                         limits = qlim, oob = scales::squish,
                         begin = 0.05, end = 0.98,
                         name = '') +
  facet_wrap(~ gene, nrow = 1) +
  coord_fixed() +
  #labs(title    = sprintf("%s — astrocyte expression per cell", TAG),
  #     subtitle = "colour = log(expression + 0.05); grey = zero-count cells; line = tissue boundary",
  #     x = "x (mm)", y = "y (mm)") +
  labs(title    = '',
       x = "x (mm)", y = "y (mm)") +
  theme_nature()+
  theme(
    axis.text = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 15),                
    legend.title = element_text(size = 11),                
    legend.key.height = unit(0.9, "cm")                   
  )




ggsave(file.path(FIG_DIR, "01_cell_level_expression.png"), p_cells,
       width = 14, height = 3.6, dpi = 300)

p_cells
 


# 1b. grid-level astrocyte counts + per-gene summed expression -----------------

## ---- rebuild the binning columns on `grid` ---------------------------------
# lattice geometry (from the grid's own i/j indices)
ux <- tapply(grid$centroid_x, grid$i, function(v) v[1])
uy <- tapply(grid$centroid_y, grid$j, function(v) v[1])
ii <- as.integer(names(ux)); xi <- as.numeric(ux)
jj <- as.integer(names(uy)); yj <- as.numeric(uy)
step_x <- median(diff(xi[order(ii)])); step_y <- median(diff(yj[order(jj)]))
x0 <- xi[which.min(ii)] - min(ii) * step_x
y0 <- yj[which.min(jj)] - min(jj) * step_y

# assign each astrocyte to a grid cell (recompute — do NOT reuse a stale grid_row)
ci <- round((cells$x - x0) / step_x)
cj <- round((cells$y - y0) / step_y)
inside <- ci >= min(ii) & ci <= max(ii) & cj >= min(jj) & cj <= max(jj)
gkey <- grid$i * 100000L + grid$j
ckey <- rep(NA_real_, nrow(cells))
ckey[inside] <- ci[inside] * 100000L + cj[inside]
cells$grid_row <- match(ckey, gkey)
valid <- !is.na(cells$grid_row)

# Y = astrocyte count per grid
grid$Y <- tabulate(cells$grid_row[valid], nbins = nrow(grid))

# M_<gene> = summed expression per grid
for (g in GENES) {
  mk <- tapply(cells[[g]][valid], cells$grid_row[valid], sum)
  v  <- rep(0, nrow(grid)); v[as.integer(names(mk))] <- as.numeric(mk)
  grid[[paste0("M_", g)]] <- v
}

# per-grid exposure: interior = full area; boundary = clipped area from the CSV
full_area      <- median(grid$area[grid$type == "interior"])
grid$area_use  <- ifelse(grid$type == "interior", full_area, grid$area)
grid$area_frac <- grid$area_use / full_area          # bounded 0-1 exposure

# verify
stopifnot(length(setdiff(c("Y","area_use", paste0("M_", GENES)), names(grid))) == 0)
cat(sprintf("binned %d/%d astrocytes | grids=%d | sum(Y)=%d\n",
            sum(valid), nrow(cells), nrow(grid), sum(grid$Y)))

## ---------------------------------------------------------------------------
## Grid-level: area, astrocyte count, per-gene summed expression (ONE row)
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## Grid-level: astrocyte count + per-gene summed expression (ONE row)
##   - grid area (exposure) panel removed
##   - larger x / y axis fonts
##   - ONE y-axis label (on the count panel); gene panels share it
##   - shared magma scale + a single collected colourbar for the 5 genes
## ---------------------------------------------------------------------------
library(scales)
library(patchwork)

theme_row <- function()
  theme_nature(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text         = element_text(size = 14, colour = "black"),
    axis.title        = element_text(size = 18, colour = "black"),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(0.30, "cm"),
    legend.key.height = unit(0.55,  "cm")
  )

## ---- magma field plotter ---------------------------------------------------
plot_grid_magma <- function(values, title, legend = "", log_scale = FALSE,
                            limits = NULL, show_y = FALSE,
                            direction = 1, trans_pad = 0.05) {
  d <- grid
  if (length(values) != nrow(d))
    stop(sprintf("plot_grid_magma: %d values but grid has %d rows (binning not run?)",
                 length(values), nrow(d)))
  d$val <- if (log_scale) log(values + trans_pad) else values
  
  qlim <- if (is.null(limits))
    stats::quantile(d$val, c(0.02, 0.98), na.rm = TRUE) else limits
  
  ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = val)) +
    geom_path(data = bnd_df, aes(X, Y, group = grp),
              colour = "grey20", linewidth = 0.35) +
    scale_fill_viridis_c(option = "A", direction = direction,
                         limits = qlim, oob = scales::squish,
                         begin = 0.05, end = 0.98, name = legend,
                         breaks = scales::pretty_breaks(n = 3)) +
    coord_fixed() +
    labs(title = title,
         x = "x (mm)",
         y = if (show_y) "y (mm)" else NULL) +
    theme_row() +
    { if (!show_y) theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank()) else NULL }
}

## ---- (i) astrocyte count per grid  (carries the single y label) ------------
p_count <- plot_grid_magma(grid$Y, "Astrocytes",
                           legend = "cells", show_y = TRUE)

## ---- (ii) per-gene summed expression, SHARED log scale ---------------------
qlim_all <- stats::quantile(
  log(unlist(grid[paste0("M_", GENES)]) + 0.05), c(0.02, 0.98), na.rm = TRUE)

# first gene panel keeps the y axis; the rest drop it
p_expr_grid <- lapply(seq_along(GENES), function(k) {
  g <- GENES[k]
  plot_grid_magma(grid[[paste0("M_", g)]], g,
                  legend = NULL, log_scale = TRUE,
                  limits = qlim_all, show_y = (k == 1))
})

## ---- combine: one row; the 5 gene panels share ONE colourbar ---------------
p_genes <- patchwork::wrap_plots(p_expr_grid, nrow = 1) +
  patchwork::plot_layout(guides = "collect")

p_grid_all <- patchwork::wrap_plots(list(p_count, p_genes),
                                    nrow = 1, widths = c(1, 5))

ggsave(file.path(FIG_DIR, "02_grid_level_counts_expression.png"), p_grid_all,
       width = 24, height = 5.5, dpi = 300)

p_grid_all


# scaled by the cell counts: 

library(scales)
library(patchwork)

theme_row <- function()
  theme_nature(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text         = element_text(size = 14, colour = "black"),
    axis.title        = element_text(size = 18, colour = "black"),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(0.30, "cm"),
    legend.key.height = unit(0.55, "cm")
  )

## ---- magma field plotter ---------------------------------------------------
plot_grid_magma <- function(values, title, legend = "", log_scale = FALSE,
                            limits = NULL, show_y = FALSE,
                            direction = 1, trans_pad = 0.05) {
  d <- grid
  if (length(values) != nrow(d))
    stop(sprintf("plot_grid_magma: %d values but grid has %d rows (binning not run?)",
                 length(values), nrow(d)))
  d$val <- if (log_scale) log(values + trans_pad) else values
  
  qlim <- if (is.null(limits))
    stats::quantile(d$val, c(0.02, 0.98), na.rm = TRUE) else limits
  
  ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = val)) +
    geom_path(data = bnd_df, aes(X, Y, group = grp),
              colour = "grey20", linewidth = 0.35) +
    scale_fill_viridis_c(option = "A", direction = direction,
                         limits = qlim, oob = scales::squish,
                         begin = 0.05, end = 0.98, name = legend,
                         na.value = "grey90",
                         breaks = scales::pretty_breaks(n = 3)) +
    coord_fixed() +
    labs(title = title,
         x = "x (mm)",
         y = if (show_y) "y (mm)" else NULL) +
    theme_row() +
    { if (!show_y) theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank()) else NULL }
}

## ---- (i) astrocyte count per grid  (carries the single y label) ------------
p_count <- plot_grid_magma(grid$Y, "Astrocytes",
                           legend = "cells", show_y = TRUE)

## ---- (ii) per-gene expression PER CELL, SHARED log scale ------------------
PAD <- 0.05                                   # same offset as the model

## expression density = summed counts / cells in that grid; NA where no cells
expr_per_cell <- function(g) {
  m <- grid[[paste0("M_", g)]]
  ifelse(grid$Y > 0, m / grid$Y, NA_real_)
}

ratio_list <- lapply(GENES, expr_per_cell)
names(ratio_list) <- GENES

qlim_all <- stats::quantile(
  log(unlist(ratio_list) + PAD), c(0.02, 0.98), na.rm = TRUE)

LEG_LAB <- expression(log(counts/cell + 0.05))

LEG_LAB <- ''


# first gene panel keeps the y axis; the rest drop it
p_expr_grid <- lapply(seq_along(GENES), function(k) {
  plot_grid_magma(ratio_list[[k]], GENES[k],
                  legend = LEG_LAB, log_scale = TRUE, trans_pad = PAD,
                  limits = qlim_all, show_y = (k == 1))
})

## ---- combine: one row; the 5 gene panels share ONE colourbar ---------------
p_genes <- patchwork::wrap_plots(p_expr_grid, nrow = 1) +
  patchwork::plot_layout(guides = "collect")

p_grid_all <- patchwork::wrap_plots(list(p_count, p_genes),
                                    nrow = 1, widths = c(1, 5))

ggsave(file.path(FIG_DIR, "02_grid_level_counts_expression.png"), p_grid_all,
       width = 24, height = 5.5, dpi = 300)

p_grid_all




## ---------------------------------------------------------------------------
## 3. Shared Stan inputs (identical for both models)
## ---------------------------------------------------------------------------
cx <- (min(grid$centroid_x) + max(grid$centroid_x)) / 2
cy <- (min(grid$centroid_y) + max(grid$centroid_y)) / 2
gx <- (grid$centroid_x - cx) / coord_scale
gy <- (grid$centroid_y - cy) / coord_scale
coords <- cbind(gx, gy); colnames(coords) <- NULL

N  <- nrow(grid)
log_grid_area <- log(grid$area_use / area_scale)
L      <- c(max(abs(gx)), max(abs(gy))) * boundary_choice
mrange <- max(fields::rdist(coords)) / (2 * 2.75)
m1 <- basis_choice; mstar <- m1 * m1
S  <- as.matrix(expand.grid(1:m1, 1:m1)); colnames(S) <- NULL

make_input <- function(Y, M) list(
  n = N, y1 = Y, y2 = M, d = 2, mstar = mstar, coords = coords,
  log_grid_area = log_grid_area, mrange = mrange,
  indices = S, L = L, is_centerted_PHI = 0)

init_full <- function() list(Astar = matrix(c(1,0,0.5,1),2,2), sigma = c(1,1),
                             beta0 = 0, beta1 = -1, ell = c(1,1),
                             betab1 = rep(0,mstar), betab2 = rep(0,mstar))
init_nest <- function() list(sigma = c(1,1), beta0 = 0, beta1 = -1, ell = c(1,1),
                             betab1 = rep(0,mstar), betab2 = rep(0,mstar))

mod_full <- cmdstan_model(STAN_FULL, compile = TRUE)
mod_nest <- cmdstan_model(STAN_NEST, compile = TRUE)

## ---------------------------------------------------------------------------
## Matern-3/2 correlation + practical ranges
##   rho_k(r)   = (1 + sqrt(3) r/ell_k) * exp(-sqrt(3) r/ell_k)
##   expression field decay:
##     C(r) = ( a21^2 rho_1(r) + a22^2 rho_2(r) ) / ( a21^2 + a22^2 )
##   ranges: solve rho_1(r1) = 0.05  and  C(r2) = 0.05
## ---------------------------------------------------------------------------
matern32 <- function(r, ell) (1 + sqrt(3) * r / ell) * exp(-sqrt(3) * r / ell)

solve_range <- function(f, upper) {
  if (f(1e-8) <= 0.05) return(0)
  if (f(upper)  > 0.05) return(NA_real_)          # decay never reaches 0.05 in range
  stats::uniroot(function(r) f(r) - 0.05, c(1e-8, upper))$root
}


###############################################################################
# S3R1 — FIT ONLY. Runs DenMark (full + nested a21=0) for each of the 5 genes
# and saves one RDS per gene containing everything the figures need.
# No plotting here — see S3R1-figures-from-rds.R for that.
#
# Prereqs already in the session (from the setup part of S3R1-DenMark-5genes.R):
#   cells, grid (WITH Y / M_* / area_use), GENES, TAG, N, coords,
#   log_grid_area, make_input(), init_full(), init_nest(),
#   mod_full, mod_nest, matern32(), solve_range(),
#   ITER_WARMUP, ITER_SAMPLING, SCRATCH, OUT_DIR
###############################################################################

stopifnot(length(setdiff(c("Y","area_use", paste0("M_", GENES)), names(grid))) == 0)
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

fit_gene <- function(g) {
  cat(sprintf("\n=== %s / %s ===\n", TAG, g))
  M     <- grid[[paste0("M_", g)]]
  input <- make_input(grid$Y, M)
  
  # cmdstanr does NOT create output_dir — make the per-gene folders first
  dir_f <- file.path(SCRATCH, paste0("full_", g))
  dir_n <- file.path(SCRATCH, paste0("nest_", g))
  dir.create(dir_f, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_n, showWarnings = FALSE, recursive = TRUE)
  
  ## ---- fit both models ----------------------------------------------------
  fit_f <- mod_full$sample(data = input, chains = 4, parallel_chains = 4,
                           iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                           adapt_delta = 0.95, init = init_full, max_treedepth = 10,
                           seed = 2026, refresh = 200, output_dir = dir_f)
  fit_n <- mod_nest$sample(data = input, chains = 4, parallel_chains = 4,
                           iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                           adapt_delta = 0.95, init = init_nest, max_treedepth = 10,
                           seed = 2026, refresh = 200, output_dir = dir_n)
  
  ## ---- (a) posterior summary ----------------------------------------------
  s_rho <- fit_f$summary(variables = c("Astar[2,1]","sigma[1]","sigma[2]",
                                       "ell[1]","ell[2]","beta0","beta1"),
                         mean = mean, sd = sd,
                         q2.5  = ~quantile(.x, 0.025), q50 = ~quantile(.x, 0.5),
                         q97.5 = ~quantile(.x, 0.975),
                         rhat = posterior::rhat,
                         ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  s_rho$gene <- g; s_rho$slice <- TAG
  cat("--- posterior summary (full model) ---\n"); print(as.data.frame(s_rho), digits = 3)
  
  dg <- fit_f$diagnostic_summary()
  cat(sprintf("    divergences=%d  max_treedepth=%d  time=%.1f min\n",
              sum(dg$num_divergent), sum(dg$num_max_treedepth),
              fit_f$time()$total / 60))
  
  ## ---- coregionalization draws:  a11=s1, a21=s2*rho, a22=s2*sqrt(1-rho^2) --
  dr   <- fit_f$draws(format = "draws_matrix")
  rho  <- as.numeric(dr[, "Astar[2,1]"])
  sig1 <- as.numeric(dr[, "sigma[1]"]);  sig2 <- as.numeric(dr[, "sigma[2]"])
  b0   <- as.numeric(dr[, "beta0"]);     b1   <- as.numeric(dr[, "beta1"])
  ell1 <- as.numeric(dr[, "ell[1]"]);    ell2 <- as.numeric(dr[, "ell[2]"])
  a11  <- sig1
  a21  <- sig2 * rho
  a22  <- sig2 * sqrt(pmax(0, 1 - rho^2))
  W1 <- dr[, paste0("w1[", 1:N, "]")]
  W2 <- dr[, paste0("w2[", 1:N, "]")]
  
  # posterior-mean fields (what the decomposition maps plot)
  a11w1   <- colMeans(a11 * W1)
  a21w1   <- colMeans(a21 * W1)
  a22w2   <- colMeans(a22 * W2)
  loglam1 <- colMeans(sweep(a11 * W1, 1, b0, "+")) + log_grid_area
  loglam2 <- colMeans(sweep(a21 * W1 + a22 * W2, 1, b1, "+")) + log(grid$Y + 0.05)
  
  ## ---- (c) decay curve + practical ranges ---------------------------------
  dmax <- max(fields::rdist(coords))
  rseq <- seq(0, dmax, length.out = 200)
  nd   <- length(rho)
  Cmat <- matrix(NA_real_, nd, length(rseq))
  for (k in seq_len(nd)) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    Cmat[k, ] <- (w21 * matern32(rseq, ell1[k]) + w22 * matern32(rseq, ell2[k])) / (w21 + w22)
  }
  decay <- data.frame(r   = rseq,
                      med = apply(Cmat, 2, median),
                      lo  = apply(Cmat, 2, quantile, 0.025),
                      hi  = apply(Cmat, 2, quantile, 0.975))
  
  r1_draws <- vapply(seq_len(nd), function(k)
    solve_range(function(r) matern32(r, ell1[k]), upper = 10 * dmax), numeric(1))
  r2_draws <- vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    solve_range(function(r) (w21*matern32(r, ell1[k]) + w22*matern32(r, ell2[k])) / (w21 + w22),
                upper = 10 * dmax)
  }, numeric(1))
  
  ## ---- (d) WAIC: full vs nested -------------------------------------------
  w_f <- loo::waic(fit_f$draws("log_lik", format = "draws_matrix"))
  w_n <- loo::waic(fit_n$draws("log_lik", format = "draws_matrix"))
  cmp <- loo::loo_compare(list(full = w_f, nested_a21_0 = w_n))
  cat("--- WAIC (full vs nested a21 = 0) ---\n"); print(cmp)
  
  ## ---- SAVE everything the figures need -----------------------------------
  out <- list(
    slice = TAG, gene = g,
    draws  = data.frame(rho = rho, sigma1 = sig1, sigma2 = sig2,
                        beta0 = b0, beta1 = b1, ell1 = ell1, ell2 = ell2),
    fields = data.frame(Y = grid$Y, M = M,
                        loglam1 = loglam1, loglam2 = loglam2,
                        a11w1 = a11w1, a21w1 = a21w1, a22w2 = a22w2),
    geom   = grid[, c("i","j","centroid_x","centroid_y","type",
                      "xmin","xmax","ymin","ymax")],
    log_grid_area = log_grid_area,
    summary = s_rho, decay = decay, r1 = r1_draws, r2 = r2_draws,
    waic_full = w_f, waic_nested = w_n, compare = cmp,
    diagnostics = dg)
  
  saveRDS(out, file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))
  cat(sprintf("    saved -> fit_%s_%s.rds\n", TAG, g))
  invisible(out)
}

## ---------------------------------------------------------------------------
## Run. Fit ALL five genes:
## ---------------------------------------------------------------------------
res <- lapply(GENES, fit_gene)
names(res) <- GENES

## ...or fit ONE gene at a time (re-runnable; skips nothing, just refits that gene):
#fit_gene("Aqp4")
#fit_gene("Sox9")
#fit_gene("Fgfr3")
#fit_gene("Ddr2")
#fit_gene("Celsr2")
 
f <- res[["Aqp4"]]        # same as readRDS(file.path(OUT_DIR, "fit_S3R1_Aqp4.rds"))
names(res[["Aqp4"]])
# slice, gene, draws, fields, geom, log_grid_area, summary, decay,
# r1, r2, waic_full, waic_nested, compare, diagnostics
sapply(res, function(f) sum(f$diagnostics$num_divergent))

input <- make_input(grid$Y, grid[["M_Aqp4"]])
dir_c <- file.path(SCRATCH, "full_Aqp4_ad99"); dir.create(dir_c, showWarnings = FALSE, recursive = TRUE)

fit_check <- mod_full$sample(data = input, chains = 4, parallel_chains = 4,
                             iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                             adapt_delta = 0.99, init = init_full, max_treedepth = 12,
                             seed = 2026, refresh = 500, output_dir = dir_c)

cat("divergences:", sum(fit_check$diagnostic_summary()$num_divergent), "\n")
print(fit_check$summary(c("Astar[2,1]","sigma[1]","sigma[2]","ell[1]","ell[2]")), digits = 3)
# compare against: res[["Aqp4"]]$summary

library(loo)

## y1 / y2 / input — needed by the comparison code below
## (input already exists from fit_check, but define y1/y2 explicitly)
g  <- "Aqp4"
y1 <- grid$Y
y2 <- grid[[paste0("M_", g)]]

## ---- the only new fit: nested model, matched settings ---------------------
dir_n <- file.path(SCRATCH, paste0("nest_", g, "_ad99"))
dir.create(dir_n, showWarnings = FALSE, recursive = TRUE)

fit_n99 <- mod_nest$sample(data = input, chains = 4, parallel_chains = 4,
                           iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                           adapt_delta = 0.99, init = init_nest, max_treedepth = 12,
                           seed = 2026, refresh = 500, output_dir = dir_n)

cat(sprintf("nested divergences: %d\n", sum(fit_n99$diagnostic_summary()$num_divergent)))
print(fit_n99$summary(c("sigma[1]","sigma[2]","ell[1]","ell[2]","beta0","beta1")), digits = 3)


library(loo); library(ggplot2); library(patchwork)

y1 <- grid$Y; y2 <- grid[["M_Aqp4"]]

## ---- combined log_lik (what the pipeline reports) --------------------------
w_f <- loo::waic(fit_check$draws("log_lik", format = "draws_matrix"))
w_n <- loo::waic(fit_n99 $draws("log_lik", format = "draws_matrix"))
cmp_comb <- loo::loo_compare(list(full = w_f, nested = w_n))

## ---- expression-only, reconstructed from draws (a21 only touches lambda2) --
ll2_full <- function(fit) {
  d <- fit$draws(c("Astar[2,1]","sigma","beta1","w1","w2"), format = "draws_matrix")
  rho <- as.numeric(d[,"Astar[2,1]"]); s2 <- as.numeric(d[,"sigma[2]"]); b1 <- as.numeric(d[,"beta1"])
  W1 <- as.matrix(d[,paste0("w1[",1:N,"]")]); W2 <- as.matrix(d[,paste0("w2[",1:N,"]")])
  lam2 <- sweep(sweep(s2*rho*W1 + s2*sqrt(pmax(0,1-rho^2))*W2, 1, b1, "+"), 2, log(y1+0.05), "+")
  dpois(matrix(y2, nrow(lam2), ncol(lam2), byrow = TRUE), exp(lam2), log = TRUE)
}
ll2_nest <- function(fit) {
  d <- fit$draws(c("sigma","beta1","w2"), format = "draws_matrix")
  s2 <- as.numeric(d[,"sigma[2]"]); b1 <- as.numeric(d[,"beta1"]); W2 <- as.matrix(d[,paste0("w2[",1:N,"]")])
  lam2 <- sweep(sweep(s2*W2, 1, b1, "+"), 2, log(y1+0.05), "+")
  dpois(matrix(y2, nrow(lam2), ncol(lam2), byrow = TRUE), exp(lam2), log = TRUE)
}
w2_f <- loo::waic(ll2_full(fit_check))
w2_n <- loo::waic(ll2_nest(fit_n99))
cmp_y2 <- loo::loo_compare(list(full = w2_f, nested = w2_n))

## ---- table -----------------------------------------------------------------
tab <- data.frame(
  comparison = c("Combined (y1+y2)", "Expression only (y2)"),
  elpd_diff  = c(cmp_comb[2,"elpd_diff"], cmp_y2[2,"elpd_diff"]),
  se_diff    = c(cmp_comb[2,"se_diff"],   cmp_y2[2,"se_diff"]))
tab$ratio   <- abs(tab$elpd_diff) / tab$se_diff
tab$verdict <- ifelse(tab$ratio > 2, "full model favoured", "indistinguishable")
print(tab, digits = 3, row.names = FALSE)

## ---- figure ----------------------------------------------------------------
tab$comparison <- factor(tab$comparison, levels = rev(tab$comparison))
p <- ggplot(tab, aes(elpd_diff, comparison, colour = verdict)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#0072B2", alpha = 0.05) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_pointrange(aes(xmin = elpd_diff - 2*se_diff, xmax = elpd_diff + 2*se_diff),
                  linewidth = 0.7, size = 0.7) +
  scale_colour_manual(values = c("full model favoured" = "#0072B2",
                                 "indistinguishable"   = "grey55"), name = NULL) +
  labs(title = "Aqp4 — full vs nested (a21 = 0)",
       subtitle = "elpd difference (full - nested) \u00b1 2 SE; < 0 favours the full model",
       x = "elpd difference", y = NULL) +
  theme_nature(base_size = 12) + theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "waic_aqp4_combined_vs_y2.png"), p, width = 7, height = 3, dpi = 300)
p




###############################################################################
# S3R1 — compare + visualize the 5 fitted genes (reads the saved RDS files)
###############################################################################
library(ggplot2); library(patchwork); library(scales)

## ---- 1. Load the five fits -------------------------------------------------
res <- lapply(GENES, function(g)
  readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g))))
names(res) <- GENES

okabe_ito <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00")
names(okabe_ito) <- GENES

## ---- 2. Comparison tables --------------------------------------------------
rho_tab <- do.call(rbind, lapply(res, function(f) {
  s   <- as.data.frame(f$summary)
  row <- s[s$variable == "Astar[2,1]", , drop = FALSE]
  data.frame(gene     = f$gene,
             mean     = row$mean,
             sd       = row$sd,
             q2.5     = row$`2.5%`,      # <-- backticks
             q97.5    = row$`97.5%`,     # <-- backticks
             rhat     = row$rhat,
             ess_bulk = row$ess_bulk,
             row.names = NULL)
}))

range_tab <- do.call(rbind, lapply(res, function(f) data.frame(
  gene = f$gene,
  r1_median = median(f$r1, na.rm = TRUE),
  r1_lo = quantile(f$r1, .025, na.rm = TRUE), r1_hi = quantile(f$r1, .975, na.rm = TRUE),
  r2_median = median(f$r2, na.rm = TRUE),
  r2_lo = quantile(f$r2, .025, na.rm = TRUE), r2_hi = quantile(f$r2, .975, na.rm = TRUE),
  row.names = NULL)))

waic_tab <- do.call(rbind, lapply(res, function(f) data.frame(
  gene        = f$gene,
  waic_full   = f$waic_full$estimates["waic","Estimate"],
  waic_nested = f$waic_nested$estimates["waic","Estimate"],
  dwaic       = f$waic_full$estimates["waic","Estimate"] -
    f$waic_nested$estimates["waic","Estimate"],
  elpd_diff   = f$compare[2, "elpd_diff"],
  elpd_se     = f$compare[2, "se_diff"],
  better      = rownames(f$compare)[1],
  row.names   = NULL)))

# headline table: does rho exclude 0, and does WAIC agree?
comparison <- merge(rho_tab, waic_tab, by = "gene")
comparison$rho_excludes_0 <- with(comparison, (q2.5 > 0) | (q97.5 < 0))
comparison$waic_prefers   <- comparison$better
print(comparison[, c("gene","mean","q2.5","q97.5","rho_excludes_0",
                     "dwaic","elpd_diff","elpd_se","waic_prefers","rhat")],
      digits = 3, row.names = FALSE)

write.csv(comparison, file.path(OUT_DIR, sprintf("comparison_%s.csv", TAG)), row.names = FALSE)
write.csv(range_tab,  file.path(OUT_DIR, sprintf("ranges_%s.csv", TAG)),     row.names = FALSE)





library(ggplot2); library(patchwork)

## ---- WAIC comparison table -------------------------------------------------
waic_tab <- do.call(rbind, lapply(res, function(f) data.frame(
  gene        = f$gene,
  waic_full   = f$waic_full$estimates["waic","Estimate"],
  waic_nested = f$waic_nested$estimates["waic","Estimate"],
  dwaic       = f$waic_full$estimates["waic","Estimate"] -
    f$waic_nested$estimates["waic","Estimate"],
  elpd_diff   = f$compare[2, "elpd_diff"],
  elpd_se     = f$compare[2, "se_diff"],
  rho_mean    = mean(f$draws$rho),
  rho_q2.5    = unname(quantile(f$draws$rho, 0.025)),
  rho_q97.5   = unname(quantile(f$draws$rho, 0.975)),
  row.names   = NULL)))

# a difference is only meaningful if |elpd_diff| > 2 x its SE
waic_tab$ratio   <- abs(waic_tab$elpd_diff) / waic_tab$elpd_se
waic_tab$verdict <- ifelse(waic_tab$ratio < 2, "Indistinguishable",
                           ifelse(waic_tab$elpd_diff < 0, "Full (a21 free)", "Nested (a21 = 0)"))
waic_tab$rho_excl0 <- with(waic_tab, (rho_q2.5 > 0) | (rho_q97.5 < 0))

print(waic_tab[, c("gene","waic_full","waic_nested","dwaic",
                   "elpd_diff","elpd_se","ratio","verdict",
                   "rho_mean","rho_excl0")], digits = 3, row.names = FALSE)
write.csv(waic_tab, file.path(OUT_DIR, sprintf("waic_comparison_%s.csv", TAG)),
          row.names = FALSE)

## ---- (A) elpd difference with +/- 2 SE, the principled view ----------------
verdict_cols <- c("Full (a21 free)"   = "#0072B2",
                  "Nested (a21 = 0)"  = "#D55E00",
                  "Indistinguishable" = "grey60")

waic_tab$gene <- factor(waic_tab$gene, levels = rev(GENES))

p_elpd <- ggplot(waic_tab, aes(elpd_diff, gene, colour = verdict)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#0072B2", alpha = 0.04) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_pointrange(aes(xmin = elpd_diff - 2*elpd_se,
                      xmax = elpd_diff + 2*elpd_se),
                  size = 0.6, linewidth = 0.6) +
  scale_colour_manual(values = verdict_cols, name = NULL) +
  labs(title = sprintf("%s — model comparison (WAIC)", TAG),
       subtitle = "elpd difference (full - nested) \u00b1 2 SE; < 0 favours the full model",
       x = "elpd difference", y = NULL) +
  theme_nature(base_size = 12) +
  theme(legend.position = "bottom")

## ---- (B) rho posterior, the parameter-level view ---------------------------
p_rho <- ggplot(waic_tab, aes(rho_mean, gene, colour = verdict)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(xmin = rho_q2.5, xmax = rho_q97.5),
                  size = 0.6, linewidth = 0.6) +
  scale_colour_manual(values = verdict_cols, guide = "none") +
  labs(title = "Density-expression correlation",
       subtitle = "posterior mean of \u03c1 with 95% CI",
       x = expression(rho), y = NULL) +
  theme_nature(base_size = 12)

## ---- (C) raw WAIC values, side by side -------------------------------------
waic_long <- rbind(
  data.frame(gene = waic_tab$gene, model = "Full (a21 free)",  waic = waic_tab$waic_full),
  data.frame(gene = waic_tab$gene, model = "Nested (a21 = 0)", waic = waic_tab$waic_nested))

p_raw <- ggplot(waic_long, aes(waic, gene, colour = model, shape = model)) +
  geom_line(aes(group = gene), colour = "grey75", linewidth = 0.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Full (a21 free)"  = "#0072B2",
                                 "Nested (a21 = 0)" = "#D55E00"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(title = "WAIC (lower = better)", x = "WAIC", y = NULL) +
  theme_nature(base_size = 12) +
  theme(legend.position = "bottom")

p_waic_all <- (p_elpd | p_rho) / p_raw + plot_layout(heights = c(1.2, 1))
ggsave(file.path(FIG_DIR, sprintf("06_waic_comparison_%s.png", TAG)),
       p_waic_all, width = 12, height = 8, dpi = 300)
p_waic_all




library(ggplot2); library(ggridges)

## ---- prior draws for rho under LKJ(eta = 1.2), K = 2 -----------------------
NPRIOR <- 8000
# Option A: rethinking
# prior_rho <- rethinking::rlkjcorr(NPRIOR, K = 2, eta = 1.2)[, 1, 2]
# Option B: exact, no extra package. For K = 2 the LKJ marginal is
#   p(rho) ∝ (1 - rho^2)^(eta - 1)  <=>  (rho + 1)/2 ~ Beta(eta, eta)
prior_rho <- 2 * rbeta(NPRIOR, 1.2, 1.2) - 1

## ---- stack prior + the 5 posteriors ---------------------------------------
draws_df <- do.call(rbind, lapply(res, function(f)
  data.frame(gene = f$gene, rho = f$draws$rho)))
draws_df <- rbind(data.frame(gene = "Prior", rho = prior_rho), draws_df)

# Prior on top: ggridges draws the LAST factor level at the top
lev <- rev(c("Prior", GENES))
draws_df$gene <- factor(draws_df$gene, levels = lev)

## ---- per-row summaries: mean + 95% CI -------------------------------------
summ <- do.call(rbind, lapply(split(draws_df, draws_df$gene), function(d)
  data.frame(gene  = d$gene[1],
             mean  = mean(d$rho),
             lo    = unname(quantile(d$rho, 0.025)),
             hi    = unname(quantile(d$rho, 0.975)))))
summ$gene <- factor(as.character(summ$gene), levels = lev)
summ$is_prior <- summ$gene == "Prior"

OFF <- 0.18   # how far above the baseline the interval bar sits

p_rho <- ggplot(draws_df, aes(x = rho, y = gene)) +
  # densities: prior grey, posteriors salmon
  geom_density_ridges(aes(fill = gene == "Prior"),
                      colour = "grey25", linewidth = 0.35,
                      scale = 1.5, rel_min_height = 0.002) +
  scale_fill_manual(values = c(`FALSE` = "#F08080", `TRUE` = "grey85"),
                    guide = "none") +
  # zero reference
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  # 95% CI bar with end dots + mean dot, drawn just above each baseline
  geom_segment(data = summ,
               aes(x = lo, xend = hi,
                   y = as.numeric(gene) + OFF, yend = as.numeric(gene) + OFF),
               colour = "#2C5F8A", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = summ,
             aes(x = lo,   y = as.numeric(gene) + OFF),
             colour = "#2C5F8A", size = 1.6, inherit.aes = FALSE) +
  geom_point(data = summ,
             aes(x = hi,   y = as.numeric(gene) + OFF),
             colour = "#2C5F8A", size = 1.6, inherit.aes = FALSE) +
  geom_point(data = summ,
             aes(x = mean, y = as.numeric(gene) + OFF),
             colour = "#2C5F8A", size = 2.4, inherit.aes = FALSE) +
  # "Mean=0.xx" labels
  geom_text(data = summ,
            aes(x = mean, y = as.numeric(gene) + OFF + 0.30,
                label = sprintf("Mean=%.2f", mean)),
            size = 8, hjust = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-1.05, 1.05),
                     breaks = seq(-1, 1, 0.5)) +
  labs(x = "Correlation", y = NULL) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.y  = element_text(size = 25, face = "italic", colour = "black"),
    axis.text.x  = element_text(size = 25, colour = "black"),
    axis.title.x = element_text(size = 25, colour = "black"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(file.path(FIG_DIR, sprintf("07_rho_prior_posterior_%s.png", TAG)),
       p_rho, width = 5, height = 6.5, dpi = 300)
p_rho


## ---- WAIC difference arrows (full - nested), baseline at 0 -----------------
waic_arrow <- do.call(rbind, lapply(res, function(f) data.frame(
  gene  = f$gene,
  dwaic = f$waic_full$estimates["waic","Estimate"] -
    f$waic_nested$estimates["waic","Estimate"],
  row.names = NULL)))
waic_arrow$gene <- factor(waic_arrow$gene, levels = rev(GENES))

# pad the axis so arrowheads and labels aren't clipped
xr  <- range(c(0, waic_arrow$dwaic))
pad <- 0.30 * diff(xr)

p_waic_arrow <- ggplot(waic_arrow, aes(y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = dwaic, y = gene, yend = gene),
               colour = "red", linewidth = 0.7,
               arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  geom_text(aes(x = dwaic, label = sprintf("%.1f", dwaic)),
            vjust = -0.9, hjust = 0.5, size = 8, colour = "black") +
  annotate("text", x = 0, y = length(GENES) + 0.55,
           label = "Baseline", fontface = "bold", size = 8, hjust = 0.5) +
  scale_x_continuous(limits = c(xr[1] - pad, xr[2] + pad * 0.4)) +
  coord_cartesian(clip = "off") +
  labs(title = "", x = "WAIC Difference (DenMark - Baseline)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold.italic", size = 20, hjust = 0),
    axis.text.y  = element_text(size = 20, face = "italic", colour = "black"),
    axis.text.x  = element_text(size = 20, colour = "black"),
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 20, r = 10, b = 5, l = 5)
  )

ggsave(file.path(FIG_DIR, sprintf("08_waic_arrows_%s.png", TAG)),
       p_waic_arrow, width = 5.5, height = 5, dpi = 300)
p_waic_arrow


library(ggplot2)

## dwaic = full - nested (WAIC scale; lower is better);
## verdict from the 2-SE rule on elpd_diff so ties are flagged
waic_arrow <- do.call(rbind, lapply(res, function(f) data.frame(
  gene      = f$gene,
  dwaic     = f$waic_full$estimates["waic","Estimate"] -
    f$waic_nested$estimates["waic","Estimate"],
  elpd_diff = f$compare[2, "elpd_diff"],
  se_diff   = f$compare[2, "se_diff"],
  row.names = NULL)))

waic_arrow$verdict <- with(waic_arrow, ifelse(
  abs(elpd_diff) <= 2 * se_diff, "Indistinguishable",
  ifelse(dwaic < 0, "DenMark better", "Baseline better")))
waic_arrow$gene <- factor(waic_arrow$gene, levels = rev(GENES))

cols <- c("DenMark better"    = "#B2182B",
          "Baseline better"   = "#2166AC",
          "Indistinguishable" = "grey55")

xr  <- range(c(0, waic_arrow$dwaic)); pad <- 0.30 * diff(xr)

p_waic_arrow <- ggplot(waic_arrow, aes(y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = dwaic, y = gene, yend = gene, colour = verdict),
               linewidth = 0.9, arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  geom_text(aes(x = dwaic, label = sprintf("%.1f", dwaic)),
            vjust = -0.9, hjust = 0.5, size = 8, colour = "black") +
  annotate("text", x = 0, y = length(GENES) + 0.55,
           label = "Baseline", fontface = "bold", size = 8, hjust = 0.5) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(xr[1] - pad, xr[2] + pad * 0.4)) +
  coord_cartesian(clip = "off") +
  labs(x = "WAIC Difference (DenMark - Baseline)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y     = element_text(size = 20, face = "italic", colour = "black"),
    axis.text.x     = element_text(size = 20, colour = "black"),
    axis.title.x    = element_text(size = 20, colour = "black"),
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    legend.position = "bottom",
    legend.text     = element_text(size = 15),
    plot.margin     = margin(t = 20, r = 10, b = 5, l = 5)
  )

ggsave(file.path(FIG_DIR, sprintf("08_waic_arrows_%s.png", TAG)),
       p_waic_arrow, width = 6.5, height = 5.5, dpi = 300)
p_waic_arrow
# then in geom_segment: colour = waic_arrow$col instead of "red"
#   arrows for Aqp4/Sox9 become grey = "not distinguishable from baseline"

 
library(ggplot2); library(patchwork)

## ---------------------------------------------------------------------------
## Correlation decay of the EXPRESSION latent field  a21*w1 + a22*w2
##   C(r) = [ a21^2 * rho1(r) + a22^2 * rho2(r) ] / ( a21^2 + a22^2 )
##   with rho_k(r) = Matern-3/2 = (1 + sqrt(3) r/ell_k) exp(-sqrt(3) r/ell_k)
##   range r2 solves C(r2) = 0.05
## Everything is recomputed PER DRAW from the saved RDS, so the ribbon is the
## posterior distribution of the decay curve, not a plug-in curve.
## ---------------------------------------------------------------------------
library(ggplot2); library(ggrepel)

RMAX <- 4                                   # x-axis: 0 to 6 mm
rseq <- seq(0, RMAX, length.out = 400)

matern32 <- function(r, ell) (1 + sqrt(3) * r / ell) * exp(-sqrt(3) * r / ell)
solve_range <- function(f, upper) {
  if (f(1e-8) <= 0.05) return(0)
  if (f(upper)  > 0.05) return(NA_real_)
  stats::uniroot(function(r) f(r) - 0.05, c(1e-8, upper))$root
}

## ---- decay curves + practical range r2, per gene ---------------------------
decay_all <- list(); range_all <- list()

for (g in GENES) {
  f   <- readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))
  d   <- f$draws
  a21 <- d$sigma2 * d$rho
  a22 <- d$sigma2 * sqrt(pmax(0, 1 - d$rho^2))
  nd  <- nrow(d)
  
  Cmat <- t(vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    (w21 * matern32(rseq, d$ell1[k]) + w22 * matern32(rseq, d$ell2[k])) / (w21 + w22)
  }, numeric(length(rseq))))
  
  decay_all[[g]] <- data.frame(
    gene = g, r = rseq,
    med  = apply(Cmat, 2, median),
    lo   = apply(Cmat, 2, quantile, 0.025),
    hi   = apply(Cmat, 2, quantile, 0.975))
  
  # practical range: solve C(r2) = 0.05, per draw (upper well beyond RMAX)
  r2 <- vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    solve_range(function(r) (w21 * matern32(r, d$ell1[k]) +
                               w22 * matern32(r, d$ell2[k])) / (w21 + w22), upper = 50)
  }, numeric(1))
  
  range_all[[g]] <- data.frame(
    gene = g,
    r2   = median(r2, na.rm = TRUE),
    lo   = unname(quantile(r2, 0.025, na.rm = TRUE)),
    hi   = unname(quantile(r2, 0.975, na.rm = TRUE)))
}

decay_all <- do.call(rbind, decay_all)
range_all <- do.call(rbind, range_all)

decay_all$gene <- factor(decay_all$gene, levels = GENES)
range_all$gene <- factor(range_all$gene, levels = GENES)
range_all$y    <- 0.05                       # ranges sit on the 0.05 line

okabe_ito <- setNames(c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00"), GENES)

## ---- overlay ---------------------------------------------------------------
p_overlay <- ggplot(decay_all, aes(r, med, colour = gene, fill = gene)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  # 0.05 threshold
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  # drop-lines at each gene's practical range
  geom_segment(data = range_all,
               aes(x = r2, xend = r2, y = 0, yend = 0.05, colour = gene),
               linetype = "dotted", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = range_all, aes(x = r2, y = y, colour = gene),
             size = 3, inherit.aes = FALSE) +
  # labelled ranges, repelled so they don't collide
  ggrepel::geom_text_repel(
    data = range_all,
    aes(x = r2, y = y, label = sprintf("%s: %.2f mm", gene, r2), colour = gene),
    inherit.aes = FALSE, size = 6, fontface = "italic",
    direction = "y", nudge_y = 0.13, segment.size = 0.3,
    box.padding = 0.35, min.segment.length = 0, seed = 1, show.legend = FALSE) +
  annotate("text", x = RMAX, y = 0.085, hjust = 1, size = 6, colour = "grey35",
           label = "correlation = 0.05") +
  scale_colour_manual(values = okabe_ito, name = NULL) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  scale_x_continuous(limits = c(0, RMAX), breaks = seq(0, RMAX, 1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = '',
       #sprintf("%s — spatial correlation decay of the expression field", TAG),
       #subtitle = expression("median with 95% CI;  "*C(r)==(a[21]^2*rho[1](r) + a[22]^2*rho[2](r))/(a[21]^2 + a[22]^2)),
       x = "Distance (mm)",
       y = "Correlation") +
  theme_nature(base_size = 20) +
  theme(
    legend.position      = c(0.98, 0.98),   # x, y in [0,1] panel coords
    legend.justification = c(1, 1),         # anchor the legend's top-right corner there
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 20),
    legend.key.size      = unit(0.9, "cm"),
    legend.margin        = margin(0, 0, 0, 0)
  )

ggsave(file.path(FIG_DIR, sprintf("16_decay_overlay_%s.png", TAG)),
       p_overlay, width = 7, height = 5, dpi = 300)

p_overlay
print(range_all, digits = 3)



#############################
# Decay: cell dens

## ---------------------------------------------------------------------------
## Correlation decay of the EXPRESSION latent field  a21*w1 + a22*w2
##   C(r) = (a21^2 rho1(r) + a22^2 rho2(r)) / (a21^2 + a22^2)
## plus the CELL-DENSITY latent field decay  rho1(r) = Matern-3/2(ell1),
## shown as the average across the five gene fits (pooled posterior draws).
## Ribbons are posterior (recomputed per draw from the saved RDS).
## ---------------------------------------------------------------------------
library(ggplot2)

RMAX <- 4
rseq <- seq(0, RMAX, length.out = 400)

matern32 <- function(r, ell) (1 + sqrt(3) * r / ell) * exp(-sqrt(3) * r / ell)
solve_range <- function(f, upper) {
  if (f(1e-8) <= 0.05) return(0)
  if (f(upper)  > 0.05) return(NA_real_)
  stats::uniroot(function(r) f(r) - 0.05, c(1e-8, upper))$root
}

decay_all <- list(); range_all <- list()
cell_draws <- list(); ell1_pool <- c()

for (g in GENES) {
  f   <- readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))
  d   <- f$draws
  a21 <- d$sigma2 * d$rho
  a22 <- d$sigma2 * sqrt(pmax(0, 1 - d$rho^2))
  nd  <- nrow(d)
  
  ## expression-field correlation decay, per draw
  Cmat <- t(vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    (w21 * matern32(rseq, d$ell1[k]) + w22 * matern32(rseq, d$ell2[k])) / (w21 + w22)
  }, numeric(length(rseq))))
  
  decay_all[[g]] <- data.frame(
    gene = g, r = rseq,
    med  = apply(Cmat, 2, median),
    lo   = apply(Cmat, 2, quantile, 0.025),
    hi   = apply(Cmat, 2, quantile, 0.975))
  
  r2 <- vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    solve_range(function(r) (w21 * matern32(r, d$ell1[k]) +
                               w22 * matern32(r, d$ell2[k])) / (w21 + w22), upper = 50)
  }, numeric(1))
  range_all[[g]] <- data.frame(gene = g, r2 = median(r2, na.rm = TRUE))
  
  ## cell-density-field decay = Matern-3/2 with ell1, per draw
  cell_draws[[g]] <- t(vapply(seq_len(nd),
                              function(k) matern32(rseq, d$ell1[k]),
                              numeric(length(rseq))))
  ell1_pool <- c(ell1_pool, d$ell1)
}

decay_all <- do.call(rbind, decay_all)
range_all <- do.call(rbind, range_all)

## average cell-density decay across the five genes (pool draws)
cell_all  <- do.call(rbind, cell_draws)
cell_decay <- data.frame(
  gene = "Cell density (avg)", r = rseq,
  med  = apply(cell_all, 2, median),
  lo   = apply(cell_all, 2, quantile, 0.025),
  hi   = apply(cell_all, 2, quantile, 0.975))

r2_cell <- median(vapply(ell1_pool,
                         function(e) solve_range(function(r) matern32(r, e), upper = 50),
                         numeric(1)), na.rm = TRUE)

## combine, order, legend labels (gene name + practical range)
okabe_ito <- setNames(c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00"), GENES)
CELL <- "Cell density (avg)"
lvl  <- c(GENES, CELL)
combined <- rbind(decay_all, cell_decay)
combined$gene <- factor(combined$gene, levels = lvl)

r2_by_gene <- range_all$r2[match(GENES, range_all$gene)]
lab <- c(sprintf("%s (%.2f mm)", GENES, r2_by_gene),
         sprintf("Cell density, avg (%.2f mm)", r2_cell))
names(lab) <- lvl
labels_ord <- unname(lab[lvl])

col_vals  <- c(okabe_ito, setNames("black",  CELL))
lty_vals  <- c(setNames(rep("solid", length(GENES)), GENES), setNames("longdash", CELL))
fill_vals <- c(okabe_ito, setNames("grey50", CELL))

p_overlay <- ggplot(combined, aes(r, med, colour = gene, fill = gene, linetype = gene)) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  annotate("text", x = RMAX, y = 0.085, hjust = 1, size = 6, colour = "grey35",
           label = "correlation = 0.05") +
  scale_colour_manual(values = col_vals,  breaks = lvl, labels = labels_ord, name = NULL) +
  scale_linetype_manual(values = lty_vals, breaks = lvl, labels = labels_ord, name = NULL) +
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_x_continuous(limits = c(0, RMAX), breaks = seq(0, RMAX, 1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Distance (mm)", y = "Correlation") +
  theme_nature(base_size = 20) +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 20),
    legend.key.size      = unit(0.9, "cm"),
    legend.margin        = margin(0, 0, 0, 0))
p_overlay



#################################################################
# Decay figure: gene expression decay + cell density decay ----
#################################################################
## ---------------------------------------------------------------------------
## Correlation decay: EXPRESSION latent field (per gene) + POOLED DENSITY field
##   Expression:  C(r) = [a21^2 rho1(r) + a22^2 rho2(r)] / (a21^2 + a22^2)
##   Density:     rho1(r) = Matern-3/2 with length-scale ell1 (single field)
##   rho_k(r) = (1 + sqrt(3) r/ell_k) exp(-sqrt(3) r/ell_k)
##   ranges solve C(r2) = 0.05 and rho1(r1) = 0.05
##
## w1 is driven by the SAME grid$Y in every gene's fit, so the 5 ell1 posteriors
## are repeat measurements of one field. We POOL THE DRAWS (concatenate all
## 5 x 4000 = 20000 draws), not the point estimates, so the grey band reflects
## both within-fit posterior width and between-fit disagreement.
## ---------------------------------------------------------------------------
library(ggplot2); library(ggrepel)

RMAX     <- 4
rseq     <- seq(0, RMAX, length.out = 400)
DENS_LAB <- "Pooled density"          # <- single place to edit the legend label

matern32 <- function(r, ell) (1 + sqrt(3) * r / ell) * exp(-sqrt(3) * r / ell)
solve_range <- function(f, upper) {
  if (f(1e-8) <= 0.05) return(0)
  if (f(upper)  > 0.05) return(NA_real_)
  stats::uniroot(function(r) f(r) - 0.05, c(1e-8, upper))$root
}

## ===========================================================================
## 1. Expression-field decay, per gene
## ===========================================================================
decay_all <- list(); range_all <- list()

for (g in GENES) {
  f   <- readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))
  d   <- f$draws
  a21 <- d$sigma2 * d$rho
  a22 <- d$sigma2 * sqrt(pmax(0, 1 - d$rho^2))
  nd  <- nrow(d)
  
  Cmat <- t(vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    (w21 * matern32(rseq, d$ell1[k]) + w22 * matern32(rseq, d$ell2[k])) / (w21 + w22)
  }, numeric(length(rseq))))
  
  decay_all[[g]] <- data.frame(
    gene = g, r = rseq,
    med  = apply(Cmat, 2, median),
    lo   = apply(Cmat, 2, quantile, 0.025),
    hi   = apply(Cmat, 2, quantile, 0.975))
  
  r2 <- vapply(seq_len(nd), function(k) {
    w21 <- a21[k]^2; w22 <- a22[k]^2
    solve_range(function(r) (w21 * matern32(r, d$ell1[k]) +
                               w22 * matern32(r, d$ell2[k])) / (w21 + w22), upper = 50)
  }, numeric(1))
  
  range_all[[g]] <- data.frame(
    gene = g,
    r2   = median(r2, na.rm = TRUE),
    lo   = unname(quantile(r2, 0.025, na.rm = TRUE)),
    hi   = unname(quantile(r2, 0.975, na.rm = TRUE)))
}

decay_all <- do.call(rbind, decay_all)
range_all <- do.call(rbind, range_all)

## ===========================================================================
## 2. Density-field decay: POOL the ell1 draws across genes
## ===========================================================================
# diagnostic: are the 5 ell1 estimates consistent? (they should be)
ell1_by_gene <- sapply(GENES, function(g)
  median(readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))$draws$ell1))
cat("per-gene ell1 medians (same density field -> should be similar):\n")
print(round(ell1_by_gene, 3))

ell1_pool <- unlist(lapply(GENES, function(g)
  readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))$draws$ell1))
cat(sprintf("pooled ell1: n = %d draws | median = %.3f mm | 95%% CI [%.3f, %.3f]\n",
            length(ell1_pool), median(ell1_pool),
            quantile(ell1_pool, .025), quantile(ell1_pool, .975)))

D1 <- t(vapply(ell1_pool, function(e) matern32(rseq, e), numeric(length(rseq))))
density_decay <- data.frame(
  gene = DENS_LAB, r = rseq,
  med  = apply(D1, 2, median),
  lo   = apply(D1, 2, quantile, 0.025),
  hi   = apply(D1, 2, quantile, 0.975))

r1_pool <- vapply(ell1_pool, function(e)
  solve_range(function(r) matern32(r, e), upper = 50), numeric(1))
r1_med  <- median(r1_pool, na.rm = TRUE)
cat(sprintf("pooled density range r1 = %.2f mm [%.2f, %.2f]\n",
            r1_med, quantile(r1_pool, .025, na.rm = TRUE),
            quantile(r1_pool, .975, na.rm = TRUE)))

range_density <- data.frame(gene = DENS_LAB, r2 = r1_med,
                            lo = quantile(r1_pool, .025, na.rm = TRUE),
                            hi = quantile(r1_pool, .975, na.rm = TRUE),
                            row.names = NULL)

## ===========================================================================
## 3. Assemble
## ===========================================================================
decay_plot <- rbind(decay_all, density_decay)
range_plot <- rbind(range_all, range_density)

LEVELS <- c(GENES, DENS_LAB)
decay_plot$gene <- factor(decay_plot$gene, levels = LEVELS)
range_plot$gene <- factor(range_plot$gene, levels = LEVELS)
range_plot$y    <- 0.05

pal <- setNames(c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","grey30"), LEVELS)
lw  <- setNames(c(rep(0.9, length(GENES)), 1.4), LEVELS)

## ===========================================================================
## 4. Overlay figure
## ===========================================================================
p_overlay <- ggplot(decay_plot, aes(r, med, colour = gene, fill = gene)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(aes(linewidth = gene)) +
  scale_linewidth_manual(values = lw, guide = "none") +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  # dot on the 0.05 line marking each practical range
  geom_point(data = range_plot, aes(x = r2, y = 0.05, colour = gene),
             size = 3.5, inherit.aes = FALSE) +
  # labels raised well above the curves, connected by dotted leaders
  ggrepel::geom_text_repel(
    data = range_plot,
    aes(x = r2, y = 0.05, label = sprintf("%s: %.2f mm", gene, r2), colour = gene),
    inherit.aes = FALSE, size = 7.5, fontface = "italic",
    direction = "y", nudge_y = 0.22, ylim = c(0.30, 0.85),
    segment.size = 0.4, segment.linetype = "dotted",
    box.padding = 0.6, point.padding = 0.3,
    min.segment.length = 0, seed = 1, max.overlaps = Inf, show.legend = FALSE) +
  annotate("text", x = RMAX, y = 0.10, hjust = 1, size = 7, colour = "grey35",
           label = "correlation = 0.05") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(0, RMAX), breaks = seq(0, RMAX, 1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "", x = "Distance (mm)", y = "Correlation") +
  theme_nature(base_size = 20) +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 20),
    legend.key.size      = unit(0.9, "cm"),
    legend.margin        = margin(0, 0, 0, 0)
  )

ggsave(file.path(FIG_DIR, sprintf("16_decay_overlay_%s.png", TAG)),
       p_overlay, width = 7.5, height = 5.5, dpi = 300)

p_overlay
print(range_plot, digits = 3)




library(scales)
library(sf)




## ---------------------------------------------------------------------------
## Latent-field decomposition, 2 rows x 3 columns, per gene
##   Row 1 (density):     lambda1 | log(lambda1) | a11*w1
##   Row 2 (expression):  lambda2 | log(lambda2) | a21*w1 | a22*w2
##
## From the Stan parameterisation A = diag_pre_multiply(sigma, Astar):
##   a11 = sigma1,  a21 = sigma2*rho,  a22 = sigma2*sqrt(1 - rho^2)
##   log(lambda1) = log_area + beta0 + a11*w1
##   log(lambda2) = beta1 + log(Y + 0.05) + a21*w1 + a22*w2
## ---------------------------------------------------------------------------

## ---- shared magma panel maker ---------------------------------------------
theme_dec <- function()
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5,
                                   margin = margin(b = 4)),
    legend.text     = element_text(size = 10),
    legend.title    = element_blank(),
    legend.key.width  = unit(0.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    plot.margin     = margin(2, 2, 2, 2)
  )

field_panel <- function(values, title, gridg, bnd, limits = NULL, symmetric = FALSE) {
  d <- gridg; d$val <- values
  qlim <- if (!is.null(limits)) limits else stats::quantile(d$val, c(0.02, 0.98), na.rm = TRUE)
  if (is.null(limits) && symmetric) { m <- max(abs(qlim)); qlim <- c(-m, m) }
  ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = val)) +
    geom_path(data = bnd, aes(X, Y, group = grp),
              colour = "grey25", linewidth = 0.3, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "A", direction = 1,
                         limits = qlim, oob = scales::squish,
                         begin = 0.03, end = 0.98,
                         breaks = scales::pretty_breaks(n = 4)) +
    coord_fixed() + labs(title = title) + theme_dec()
}

## ---- build one gene's decomposition ---------------------------------------
decomp_gene <- function(g) {
  f  <- readRDS(file.path(OUT_DIR, sprintf("fit_%s_%s.rds", TAG, g)))
  gg <- f$geom
  fl <- f$fields
  
  lam1 <- exp(fl$loglam1)
  lam2 <- exp(fl$loglam2)
  
  # ---- SHARED SCALES for the additive components -------------------------
  # Group A (density row):    log(lambda1) and a11*w1
  grpA <- c(fl$loglam1, fl$a11w1)
  limA <- stats::quantile(grpA, c(0.02, 0.98), na.rm = TRUE)
  
  # Group B (expression row): log(lambda2), a21*w1, a22*w2
  grpB <- c(fl$loglam2, fl$a21w1, fl$a22w2)
  limB <- stats::quantile(grpB, c(0.02, 0.98), na.rm = TRUE)
  
  # ---- row 1: density ----
  p11 <- field_panel(lam1,       expression(lambda[1](g)),      gg, bnd_df)              # own scale
  p12 <- field_panel(fl$loglam1, expression(log~lambda[1](g)),  gg, bnd_df, limits = limA)
  p13 <- field_panel(fl$a11w1,   expression(a[11]*omega[1](g)), gg, bnd_df, limits = limA)
  
  # ---- row 2: expression ----
  p21 <- field_panel(lam2,       expression(lambda[2](g)),      gg, bnd_df)              # own scale
  p22 <- field_panel(fl$loglam2, expression(log~lambda[2](g)),  gg, bnd_df, limits = limB)
  p23 <- field_panel(fl$a21w1,   expression(a[21]*omega[1](g)), gg, bnd_df, limits = limB)
  p24 <- field_panel(fl$a22w2,   expression(a[22]*omega[2](g)), gg, bnd_df, limits = limB)
  
  cap <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, hjust = 0.5, vjust = 0.5,
             label = sprintf("Decomposition\nvisualization for\nstudying %s", g),
             fontface = "bold.italic", size = 5, lineheight = 1.1)
  
  (p11 | p12 | p13 | cap) /
    (p21 | p22 | p23 | p24)
}

## ---- write one figure per gene --------------------------------------------
#for (g in GENES) {
#  p <- decomp_gene(g)
#  ggsave(file.path(FIG_DIR, sprintf("17_decomposition_%s_%s.png", TAG, g)),
#         p, width = 14, height = 7, dpi = 300)
#  message("wrote decomposition for ", g)
#}

## ---- optional: all 5 genes stacked in one tall figure ---------------------
#p_all <- patchwork::wrap_plots(lapply(GENES, decomp_gene), ncol = 1)
#ggsave(file.path(FIG_DIR, sprintf("17_decomposition_ALL_%s.png", TAG)),
#       p_all, width = 14, height = 34, dpi = 300)
#p_all

# the candidate genes: 
GENES

decomp_gene('Aqp4')
decomp_gene('Sox9')
decomp_gene('Fgfr3')
decomp_gene('Ddr2')
decomp_gene('Celsr2')





library(Seurat); library(ggplot2); library(ggrepel); library(dplyr)

## ---------------------------------------------------------------------------
## 0. Paths
## ---------------------------------------------------------------------------
RES_S3R1 <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/cluster-for-8-rest-slices/results-S3R1")

## ---------------------------------------------------------------------------
## 1. Astrocyte marker specificity via FindMarkers(test.use = "roc")
##    Run on ALL cells / ALL cell types (not just astrocytes).
##    `power` = 2*|AUC - 0.5|  ->  0 = no discrimination, 1 = perfect.
## ---------------------------------------------------------------------------
drop_cols <- c("cell_id", "cell_type", "cell_type_fine", "slice")
gene_cols <- setdiff(names(cell_by_gene), drop_cols)
gene_cols <- gene_cols[!grepl("^blank", gene_cols, ignore.case = TRUE)]

counts <- t(as.matrix(cell_by_gene[, gene_cols, drop = FALSE]))   # genes x cells
colnames(counts) <- make.unique(as.character(cell_by_gene$cell_id))

so <- CreateSeuratObject(counts = counts, project = TAG)
so <- NormalizeData(so, verbose = FALSE)
Idents(so) <- factor(cell_by_gene$cell_type)

# thresholds at 0 so every gene is scored (defaults silently drop most)
roc <- FindMarkers(so, ident.1 = "Astrocyte", test.use = "roc",
                   logfc.threshold = 0, min.pct = 0,
                   only.pos = FALSE, verbose = FALSE)
roc$gene <- rownames(roc)

spec <- roc %>%
  dplyr::transmute(gene,
                   auc         = myAUC,
                   specificity = power,        # 2*|AUC - 0.5|
                   avg_diff    = avg_diff)
cat(sprintf("Seurat scored %d genes\n", nrow(spec)))

## ---------------------------------------------------------------------------
## 2. rho for ALL genes, from the cluster per-gene summaries
##    NOTE: the CSVs have columns literally named "2.5%" / "97.5%", so we read
##    with check.names = FALSE and reference them with backticks.
## ---------------------------------------------------------------------------
files <- list.files(RES_S3R1,
                    pattern = sprintf("^denmark_summary_%s_g.*\\.csv$", TAG),
                    full.names = TRUE)
cat(sprintf("found %d per-gene summary files\n", length(files)))
stopifnot(length(files) > 0)

allsum <- do.call(rbind, lapply(files, read.csv,
                                stringsAsFactors = FALSE, check.names = FALSE))

rho_raw <- subset(allsum, variable == "Astar[2,1]")
rho_df <- data.frame(
  gene      = rho_raw$gene_name,
  rho       = rho_raw$mean,
  lo        = rho_raw$`2.5%`,
  hi        = rho_raw$`97.5%`,
  rhat      = rho_raw$rhat,
  divergent = rho_raw$num_divergent,
  row.names = NULL)

# QC: drop non-converged fits (with ~483 automated fits, some will fail)
n0  <- nrow(rho_df)
rho_df <- rho_df[!is.na(rho_df$rhat) & rho_df$rhat <= 1.05, ]
rho_df <- rho_df[!duplicated(rho_df$gene), ]
cat(sprintf("rho: %d genes (dropped %d with rhat > 1.05)\n", nrow(rho_df), n0 - nrow(rho_df)))

rho_df$sig <- (rho_df$lo > 0) | (rho_df$hi < 0)     # 95% CI excludes 0

dat <- merge(rho_df, spec, by = "gene")
cat(sprintf("merged: %d genes with both rho and specificity\n", nrow(dat)))


## ---------------------------------------------------------------------------
## 3. Classes: highlight markers; flag the CIRCULAR label-defining genes
## ---------------------------------------------------------------------------
gating_genes  <- c("Aqp4", "Gfap", "Aldh1l1")     # defined the Astrocyte label
astro_markers <- c("Aqp4", "Gfap", "Aldh1l1", "Slc1a3", "Sox9",
                   "Gja1", "Aldoc", "S100b", "Fgfr3")
FOCAL <- c("Aqp4", "Sox9", "Fgfr3", "Ddr2", "Celsr2")   # <- only these get labels

dat$gating <- dat$gene %in% gating_genes
dat$class  <- factor(ifelse(dat$gene %in% astro_markers, "Astrocyte marker", "Other gene"),
                     levels = c("Astrocyte marker", "Other gene"))

d_clean <- subset(dat, !gating)       # statistics still exclude the circular genes
ct <- suppressWarnings(cor.test(d_clean$specificity, d_clean$rho, method = "spearman"))
lab <- sprintf("Spearman \u03c1 = %.2f (p = %.3g), n = %d genes (3 label-defining genes excluded)",
               unname(ct$estimate), ct$p.value, nrow(d_clean))
cat(lab, "\n")

## ---------------------------------------------------------------------------
## 4. Scatter
## ---------------------------------------------------------------------------
p_scatter <- ggplot(dat, aes(specificity, rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  # no geom_smooth: trend line and shaded band removed
  geom_point(aes(colour = class, size = class, alpha = class), shape = 16) +   # all solid circles
  ggrepel::geom_text_repel(
    data = subset(dat, gene %in% FOCAL),          # <- only the 5 focal genes labelled
    aes(label = gene), size = 12, fontface = "italic",
    box.padding = 0.7, point.padding = 0.3,
    min.segment.length = 0, segment.size = 0.4,
    max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c("Astrocyte marker" = "#B2182B",
                                 "Other gene"       = "grey55"), name = NULL) +
  scale_size_manual(values  = c("Astrocyte marker" = 3.5, "Other gene" = 1.3),
                    guide = "none") +
  scale_alpha_manual(values = c("Astrocyte marker" = 1, "Other gene" = 0.45),
                     guide = "none") +
  labs(title = '',
       x = "Astrocyte marker specificity  (ROC power)",
       y = "density-expression correlation") +
  theme_nature(base_size = 25) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, sprintf("09_rho_vs_specificity_%s.png", TAG)),
       p_scatter, width = 7.5, height = 6.5, dpi = 300)
p_scatter






#==================================
# DenMark vs. SPARK -----
#==================================

library(SPARK); library(ggplot2); library(ggrepel); library(dplyr)

ASTAR_DIR <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/cluster-for-8-rest-slices/astar_summaries")

## ---------------------------------------------------------------------------
## 1. DCG ranking from DenMark rho
##    The astar_summary_<TAG>.csv is ALREADY one row per gene, already filtered
##    to Astar[2,1], with clean q2.5 / q97.5 names. No globbing, no subset().
## ---------------------------------------------------------------------------
dcg <- read.csv(file.path(ASTAR_DIR, sprintf("astar_summary_%s.csv", TAG)),
                stringsAsFactors = FALSE)
names(dcg)[names(dcg) == "gene_name"] <- "gene"
names(dcg)[names(dcg) == "mean"]      <- "rho"
names(dcg)[names(dcg) == "q2.5"]      <- "lo"
names(dcg)[names(dcg) == "q97.5"]     <- "hi"

dcg$abs_rho  <- abs(dcg$rho)
dcg$dcg_rank <- rank(-dcg$abs_rho, ties.method = "average")   # 1 = strongest coupling
dcg$sig_rho  <- (dcg$lo > 0) | (dcg$hi < 0)                   # 95% CI excludes 0

cat(sprintf("DCG: %d genes | rho %.2f to %.2f | %d (%.0f%%) with CI excluding 0\n",
            nrow(dcg), min(dcg$rho), max(dcg$rho),
            sum(dcg$sig_rho), 100*mean(dcg$sig_rho)))

## ---------------------------------------------------------------------------
## 2. SVG detection with SPARK-X, on the ASTROCYTE cells
##    (same population DenMark was fitted on)
## ---------------------------------------------------------------------------
cbg_a  <- cell_by_gene[cell_by_gene$cell_type == "Astrocyte", , drop = FALSE]
meta_b <- cell_location[cell_location$cell_type == "Astrocyte", , drop = FALSE]
meta_b <- meta_b[match(cbg_a$cell_id, meta_b$cell_id), , drop = FALSE]
stopifnot(all(cbg_a$cell_id == meta_b$cell_id))

gene_cols <- setdiff(names(cbg_a), c("cell_id","cell_type","cell_type_fine","slice"))
gene_cols <- gene_cols[!grepl("^blank", gene_cols, ignore.case = TRUE)]

expr_mat <- t(as.matrix(cbg_a[, gene_cols, drop = FALSE]))    # genes x cells
colnames(expr_mat) <- cbg_a$cell_id
loc_mat  <- as.matrix(meta_b[, c("center_x", "center_y")])    # raw um
rownames(loc_mat) <- meta_b$cell_id
stopifnot(identical(colnames(expr_mat), rownames(loc_mat)))

sx <- sparkx(count_in = expr_mat, locus_in = loc_mat,
              option = "mixture")

svg <- data.frame(gene     = rownames(sx$res_mtest),
                  pval     = sx$res_mtest$combinedPval,
                  adjp     = sx$res_mtest$adjustedPval,
                  row.names = NULL)
svg$svg_rank <- rank(svg$pval, ties.method = "average")       # 1 = most spatially variable
svg$is_svg   <- svg$adjp < 0.05
cat(sprintf("SVG: %d genes | %d significant (FDR < 0.05)\n", nrow(svg), sum(svg$is_svg)))

## ---------------------------------------------------------------------------
## 3. Join + correlate
## ---------------------------------------------------------------------------
dat <- merge(dcg, svg, by = "gene")
cat(sprintf("merged: %d genes with both DCG and SVG\n", nrow(dat)))

pear  <- cor.test(dat$dcg_rank, dat$svg_rank, method = "pearson")
spear <- suppressWarnings(cor.test(dat$dcg_rank, dat$svg_rank, method = "spearman"))
cap <- sprintf("Pearson r = %.3f (95%% CI %.3f\u2013%.3f, p = %.2g);  Spearman \u03c1s = %.3f;  n = %d genes",
               unname(pear$estimate), pear$conf.int[1], pear$conf.int[2],
               pear$p.value, unname(spear$estimate), nrow(dat))
cat(cap, "\n")

## ---------------------------------------------------------------------------
## 4. Rank-rank scatter
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## 4. Rank-rank scatter
## ---------------------------------------------------------------------------
FOCAL <- c("Aqp4","Sox9","Fgfr3","Ddr2","Celsr2")
dat$highlight <- dat$gene %in% FOCAL
lab_dat <- subset(dat, highlight | (dcg_rank <= 8 & svg_rank <= 8))

p_rank <- ggplot(dat, aes(dcg_rank, svg_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.4) +
  # geom_smooth removed: no regression line, no shaded band
  geom_point(aes(colour = sig_rho, size = highlight, alpha = highlight)) +
  ggrepel::geom_text_repel(data = lab_dat, aes(label = gene),
                           size = 13, fontface = "italic",
                           box.padding = 0.8, point.padding = 0.35,
                           min.segment.length = 0, segment.size = 0.4,
                           max.overlaps = Inf, seed = 1,
                           show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "\u03c1 CI excludes 0", `FALSE` = "n.s."),
                      name = NULL) +
  scale_size_manual(values  = c(`TRUE` = 3.5, `FALSE` = 1.4), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5),  guide = "none") +
  scale_x_reverse() + scale_y_reverse() +          # rank 1 (best) at top-right
  labs(title = '',
       x = "DCG rank (DenMark, 1 = strongest density-correlation)",
       y = "SVG rank (SPARK-X, 1 = most spatially variable)") +
  theme_nature(base_size = 20) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, sprintf("10_dcg_vs_svg_rank_%s.png", TAG)),
       p_rank, width = 8.5, height = 7.5, dpi = 300)
p_rank


p_rank <- ggplot(dat, aes(dcg_rank, svg_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = sig_rho, size = highlight, alpha = highlight)) +
  ggrepel::geom_text_repel(data = lab_dat, aes(label = gene),
                           size = 13, fontface = "italic",
                           box.padding = 0.8, point.padding = 0.35,
                           min.segment.length = 0, segment.size = 0.4,
                           max.overlaps = Inf, seed = 1,
                           show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"),
                      labels = c(`TRUE`  = "95% CI excludes 0",
                                 `FALSE` = "Not informative"),   
                      name = NULL) +
  scale_size_manual(values  = c(`TRUE` = 3.5, `FALSE` = 1.4), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5),  guide = "none") +
  scale_x_reverse() + scale_y_reverse() +
  labs(title = '',
       x = "DCG rank (1 = strongest density-correlation)",
       y = "SVG rank (1 = most spatially variable)") +
  theme_nature(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 20),
    axis.text       = element_text(size = 22, colour = "black"),   # tick numbers
    axis.title      = element_text(size = 24, colour = "black")    # axis labels
  )

ggsave(file.path(FIG_DIR, sprintf("10_dcg_vs_svg_rank_%s.png", TAG)),
       p_rank, width = 9.5, height = 8.5, dpi = 300)
p_rank


## ---------------------------------------------------------------------------
## 5. Quadrant counts (what the Venn conveyed, but quantitative)
## ---------------------------------------------------------------------------
dat$top_dcg <- dat$dcg_rank <= 50
dat$top_svg <- dat$svg_rank <= 50
print(table(DCG_top50 = dat$top_dcg, SVG_top50 = dat$top_svg))
write.csv(dat, file.path(OUT_DIR, sprintf("dcg_vs_svg_%s.csv", TAG)), row.names = FALSE)




#=======================================
# Slice by Slice — reproducibility panel
#=======================================
library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(ggrepel)

TAGS9     <- c("S1R1","S1R2","S1R3","S2R1","S2R2","S2R3","S3R1","S3R2","S3R3")
ASTAR_DIR <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/cluster-for-8-rest-slices/astar_summaries")

## ---------------------------------------------------------------------------
## 1. Load rho for all 9 slices
##    Each astar_summary_<TAG>.csv is ALREADY one row per gene, already filtered
##    to Astar[2,1], with columns: slice, gene_id, gene_name, mean, q2.5, q97.5
## ---------------------------------------------------------------------------
rho <- dplyr::bind_rows(lapply(TAGS9, function(tg) {
  f <- file.path(ASTAR_DIR, sprintf("astar_summary_%s.csv", tg))
  if (!file.exists(f)) { message("missing: ", basename(f)); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(slice = tg,
             gene  = d$gene_name,
             rho   = d$mean,
             lo    = d$q2.5,
             hi    = d$q97.5,
             row.names = NULL)
})) %>%
  dplyr::distinct(slice, gene, .keep_all = TRUE) %>%
  dplyr::mutate(slice = factor(slice, levels = TAGS9),
                sig   = (lo > 0) | (hi < 0))       # 95% CI excludes 0

cat(sprintf("%d slices x %d genes | %.0f%% of (slice,gene) fits have CI excluding 0\n",
            dplyr::n_distinct(rho$slice), dplyr::n_distinct(rho$gene),
            100 * mean(rho$sig)))

## wide matrix: genes x slices
wide <- rho %>% dplyr::select(gene, slice, rho) %>%
  tidyr::pivot_wider(names_from = slice, values_from = rho, values_fn = mean)
mat  <- as.matrix(wide[, -1]); rownames(mat) <- wide$gene
mat  <- mat[, TAGS9[TAGS9 %in% colnames(mat)], drop = FALSE]

## ---------------------------------------------------------------------------
## 2. (A) slice x slice concordance heatmap
## ---------------------------------------------------------------------------
cmat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

cdf <- as.data.frame(as.table(cmat))
names(cdf) <- c("slice_i", "slice_j", "r")
cdf$slice_i <- factor(cdf$slice_i, levels = TAGS9)
cdf$slice_j <- factor(cdf$slice_j, levels = rev(TAGS9))

offdiag <- mean(cmat[upper.tri(cmat)], na.rm = TRUE)   # headline reproducibility number
cat(sprintf("mean off-diagonal concordance = %.2f (range %.2f - %.2f)\n",
            offdiag, min(cmat[upper.tri(cmat)]), max(cmat[upper.tri(cmat)])))


p_heat <- ggplot(cdf, aes(slice_i, slice_j, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 9, colour = "black") +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       name = 'Correlation') +
  coord_fixed() +
  labs(title = '',
         #sprintf("Agreement of Astar[2,1] between slices (mean off-diagonal = %.2f)", offdiag),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 18, colour = "black"),
    axis.ticks       = element_blank(),
    plot.title       = element_text(size = 13, hjust = 0),
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 20),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height= unit(1.2, "cm")
  )

ggsave(file.path(FIG_DIR, "11a_slice_concordance_heatmap.png"),
       p_heat, width = 8, height = 7, dpi = 300)
p_heat


## ---------------------------------------------------------------------------
## 3. (B) STABLE CORE: consistent sign AND significance across slices
## ---------------------------------------------------------------------------
K <- 7                                   # of 9 slices (7 = strict; 5 = simple majority)

core <- rho %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_slices  = dplyr::n(),
    n_pos_sig = sum(sig & rho > 0, na.rm = TRUE),
    n_neg_sig = sum(sig & rho < 0, na.rm = TRUE),
    mean_rho  = mean(rho, na.rm = TRUE),
    sd_rho    = sd(rho,   na.rm = TRUE),
    min_rho   = min(rho,  na.rm = TRUE),
    max_rho   = max(rho,  na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(core = dplyr::case_when(
    n_pos_sig >= K & n_neg_sig == 0 ~ "Stable positive",
    n_neg_sig >= K & n_pos_sig == 0 ~ "Stable negative",
    TRUE                            ~ "Not stable"))

core$core <- factor(core$core,
                    levels = c("Stable positive", "Stable negative", "Not stable"))
print(table(core$core))
write.csv(core, file.path(OUT_DIR, "stable_core_genes.csv"), row.names = FALSE)

core_cols <- c("Stable positive" = "#B2182B",
               "Stable negative" = "#2166AC",
               "Not stable"      = "grey75")

## (B1) volcano-style: mean rho vs cross-slice consistency
p_core <- ggplot(core, aes(mean_rho, n_pos_sig - n_neg_sig, colour = core)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = c(-K, K), linetype = "dotted", colour = "grey50", linewidth = 0.3) +
  geom_point(aes(size = core, alpha = core)) +
  ggrepel::geom_text_repel(
    data = core %>% dplyr::filter(core != "Not stable") %>%
      dplyr::group_by(core) %>% dplyr::slice_max(abs(mean_rho), n = 8) %>% dplyr::ungroup(),
    aes(label = gene), size = 7, fontface = "italic",
    max.overlaps = 25, box.padding = 0.4, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = core_cols, name = NULL) +
  scale_size_manual(values  = c("Stable positive" = 2.4, "Stable negative" = 2.4,
                                "Not stable" = 1.2), guide = "none") +
  scale_alpha_manual(values = c("Stable positive" = 1, "Stable negative" = 1,
                                "Not stable" = 0.45), guide = "none") +
  labs(title ='',
         #"Stable core",
       #subtitle = sprintf("significant with a consistent sign in \u2265 %d of 9 slices", K),
       x = expression(mean~rho~across~slices),
       y = "consistency  (# slices sig. +  \u2212  # slices sig. \u2212)") +
  theme_nature(base_size = 20) + theme(legend.position = "bottom")

p_core



## (B2) the core genes, one row each, all 9 slice estimates shown
core_genes <- core %>% dplyr::filter(core != "Not stable") %>%
  dplyr::arrange(core, mean_rho) %>% dplyr::pull(gene)
cat(sprintf("stable core: %d genes\n", length(core_genes)))

core_long <- rho %>% dplyr::filter(gene %in% core_genes) %>%
  dplyr::left_join(dplyr::select(core, gene, core, mean_rho), by = "gene") %>%
  dplyr::mutate(gene = factor(gene, levels = core_genes))

p_core_genes <- ggplot(core_long, aes(rho, gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = core), size = 1.4, alpha = 0.65,
             position = position_jitter(height = 0.12, width = 0)) +
  stat_summary(aes(colour = core), fun = mean, geom = "point",
               shape = 18, size = 3.4) +
  scale_colour_manual(values = core_cols, guide = "none") +
  facet_grid(core ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Stable-core genes across all 9 slices",
       #subtitle = "small dots = per-slice \u03c1; diamond = mean",
       x = expression(rho), y = NULL) +
  theme_nature(base_size = 11) +
  theme(axis.text.y = element_text(size = 7, face = "italic"))
 
p_core_genes

#----------------------------------------
# GSEA -----
#----------------------------------------



library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
library(ggplot2); library(dplyr); library(stringr)

## ---------------------------------------------------------------------------
## 1. Ranked gene list from DenMark rho (S3R1, all 483 genes)
## ---------------------------------------------------------------------------
d <- read.csv(file.path(ASTAR_DIR, sprintf("astar_summary_%s.csv", TAG)),
              stringsAsFactors = FALSE)
geneList <- setNames(d$mean, d$gene_name)
geneList <- sort(geneList[!duplicated(names(geneList))], decreasing = TRUE)

## ---------------------------------------------------------------------------
## 2. GSEA on GO — MOUSE (org.Mm.eg.db), all three ontologies like the reference
##    ("trace-amine receptor activity" is MF, "nervous system development" is BP,
##     "cellular component organization" is CC -> the reference used ont = "ALL")
## ---------------------------------------------------------------------------
set.seed(2026)
gse <- gseGO(geneList      = geneList,
             OrgDb         = org.Mm.eg.db,      # <- MOUSE
             keyType       = "SYMBOL",
             ont           = "ALL",             # BP + MF + CC
             minGSSize     = 5,                 # panel is only ~483 genes
             maxGSSize     = 300,
             pvalueCutoff  = 0.05,
             pAdjustMethod = "BH",
             eps           = 0,
             verbose       = FALSE)

cat(sprintf("enriched sets: %d\n", nrow(as.data.frame(gse))))

## ---------------------------------------------------------------------------
## 3. THE FIGURE — enrichplot::dotplot with activated/suppressed facets
##    This is exactly the reference style; split by NES sign.
## ---------------------------------------------------------------------------
p_dot <- dotplot(gse,
                 showCategory = 5,                 # per facet
                 split        = ".sign",           # activated vs suppressed
                 x            = "GeneRatio",
                 color        = "p.adjust",
                 size         = "Count") +
  facet_grid(. ~ .sign) +
  scale_colour_gradient(low = "#2166AC", high = "#B2182B", name = "p.adjust") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 28)) +
  labs(x = "GeneRatio", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.y      = element_text(size = 20, colour = "black"),
    axis.text.x      = element_text(size = 20, colour = "black"),
    axis.title.x     = element_text(size = 20),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 15)
  )

ggsave(file.path(FIG_DIR, sprintf("12_gsea_dotplot_%s.png", TAG)),
       p_dot, width = 8.5, height = 5, dpi = 300)
p_dot





## ===========================================================================
## Some other tasks: (1) 
# Mertk expression by cell type (S3R1) — confirm astrocyte/microglia, not neuron
## ===========================================================================
library(dplyr)

GENE <- "Mertk"

## cell_by_gene already loaded earlier with cell IDs as character.
## If not in the session, reload:
if (!exists("cell_by_gene")) {
  read_idchar <- function(p) {
    h  <- read.csv(p, nrows = 1, check.names = FALSE)
    cc <- rep(NA, ncol(h)); cc[1] <- "character"
    read.csv(p, colClasses = cc, check.names = FALSE)
  }
  cell_by_gene <- read_idchar(CBG_CSV); names(cell_by_gene)[1] <- "cell_id"
}
stopifnot(GENE %in% names(cell_by_gene))          # gene must be in the panel

## --- mean RAW count per cell type -------------------------------------------
tab <- cell_by_gene %>%
  group_by(cell_type) %>%
  summarise(n_cells    = n(),
            mean_expr  = mean(.data[[GENE]]),                 # mean raw count / cell
            pct_pos    = 100 * mean(.data[[GENE]] > 0),       # % of cells detecting it
            median_expr= median(.data[[GENE]]),
            .groups = "drop") %>%
  arrange(desc(mean_expr))

cat(sprintf("\n=== %s expression by cell type (S3R1) ===\n", GENE))
print(as.data.frame(tab), digits = 3, row.names = FALSE)

## --- normalized cross-check (mean of log1p, less abundance-driven) ----------
tab_norm <- cell_by_gene %>%
  group_by(cell_type) %>%
  summarise(mean_log1p = mean(log1p(.data[[GENE]])), .groups = "drop") %>%
  arrange(desc(mean_log1p))
cat("\n-- mean log1p(count) per cell type --\n")
print(as.data.frame(tab_norm), digits = 3, row.names = FALSE)

## --- explicit verdict vs neurons --------------------------------------------
mv <- setNames(tab$mean_expr, tab$cell_type)
cat(sprintf("\nAstrocyte / Neuron ratio = %.2f\nMicroglia / Neuron ratio = %.2f\n",
            mv["Astrocyte"] / mv["Neuron"], mv["Microglia"] / mv["Neuron"]))



