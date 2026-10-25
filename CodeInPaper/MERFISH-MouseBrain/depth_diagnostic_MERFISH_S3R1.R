## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds Revision-R1/Reanalysis-MERFISH/tissue-boundary-9-slices/annotation_out/
OUTPUT_DIR <- "output"   # the depth diagnostic figure and table go to Revision-R1/seq-depth-diagnostics/
REPO_ROOT  <- "."        # root of the DenMark checkout (the grid file is read from CodeInPaper/)
dir.create(file.path(OUTPUT_DIR, "Revision-R1/seq-depth-diagnostics"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

# =====================================================================
# Depth diagnostic - MERFISH mouse brain, slice S3R1
# (A) local cell density map  (B) mean depth map
# (C) astrocyte per-cell depth vs local ASTROCYTE density, styled as Fig. S12B
# density on the interior 30x30 alpha grid (constant area, no boundary leverage)
# =====================================================================

## ---- paths (edit to your machine) ----
base_dir  <- file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/tissue-boundary-9-slices/annotation_out")
gene_path <- file.path(base_dir, "cell_by_gene_with_celltype_S3R1.csv")
meta_path <- file.path(base_dir, "cell_metadata_with_celltype_S3R1.csv")
grid_path <- file.path(REPO_ROOT, "CodeInPaper", "MERFISH-MouseBrain", "grid_area", "grid_area_alpha_S3R1.csv")
out_dir   <- file.path(OUTPUT_DIR, "Revision-R1/seq-depth-diagnostics")

## cell_id values are 39-digit integers: read as CHARACTER or the merge breaks.
cg   <- read.csv(gene_path, check.names = FALSE, colClasses = c(cell_id = "character"))
meta <- read.csv(meta_path, colClasses = c(cell_id = "character"))
grid <- read.csv(grid_path)

## ---------------------------------------------------------------
## 1. per-cell detection depth = sum over gene columns
## ---------------------------------------------------------------
meta_cols <- c("cell_id", "cell_type", "cell_type_fine")
genes     <- setdiff(colnames(cg),
                     c(meta_cols, grep("^Blank[-_]", colnames(cg), value = TRUE)))
depth_df  <- data.frame(cell_id = cg$cell_id, depth = rowSums(cg[, genes]))
cat(sprintf("depth: %d genes summed | median %.0f  mean %.1f\n",
            length(genes), median(depth_df$depth), mean(depth_df$depth)))

## ---------------------------------------------------------------
## 2. assemble cells (coords in um), join depth by cell_id
## ---------------------------------------------------------------
d <- merge(data.frame(cell_id = meta$cell_id, cx = meta$center_x, cy = meta$center_y,
                      celltype = meta$cell_type_fine),
           depth_df, by = "cell_id")
d <- d[is.finite(d$depth), ]

## ---------------------------------------------------------------
## 3. assign each cell to a grid bin, matching the provided 30x30 grid
## ---------------------------------------------------------------
cxs <- sort(unique(grid$centroid_x)); cys <- sort(unique(grid$centroid_y))
wx  <- (max(cxs) - min(cxs)) / (length(cxs) - 1)
wy  <- (max(cys) - min(cys)) / (length(cys) - 1)
d$i <- pmin(pmax(round((d$cx - min(cxs)) / wx), 0), length(cxs) - 1)
d$j <- pmin(pmax(round((d$cy - min(cys)) / wy), 0), length(cys) - 1)
d$gk          <- paste(d$i, d$j, sep = "_")
grid$gk       <- paste(grid$i, grid$j, sep = "_")
grid$area_mm2 <- grid$area / 1e6

## INTERIOR grids only -> constant area (0.0614 mm^2)
gint <- grid[grid$type == "interior", ]
area <- setNames(gint$area_mm2, gint$gk)
d    <- d[d$gk %in% names(area), ]

ncell <- tapply(rep(1, nrow(d)), d$gk, sum)
dens  <- ncell[names(area)] / area                 # cells / mm^2 (all cells)
names(dens) <- names(area)
dens  <- dens[is.finite(dens)]

## ---------------------------------------------------------------
## 4. astrocytes: per-grid mean astrocyte depth vs astrocyte density
##    (x = astrocytes per quadrat, as in Fig. S12B; interior quadrats
##     all have the same area, so this is proportional to cells/mm^2)
## ---------------------------------------------------------------
ast    <- d[d$celltype == "Astrocyte", ]
ast_n  <- tapply(rep(1, nrow(ast)), ast$gk, sum)   # astrocytes / quadrat
ast_md <- tapply(ast$depth, ast$gk, mean)          # mean astrocyte depth / quadrat
q      <- intersect(names(ast_n), names(area))
A_D    <- as.numeric(ast_n[q])
A_M    <- as.numeric(ast_md[q])
ok     <- is.finite(A_D) & is.finite(A_M); A_D <- A_D[ok]; A_M <- A_M[ok]

A_mean <- mean(ast$depth)                           # mean per-cell astrocyte depth
A_r    <- cor(A_D, A_M)
A_rho  <- suppressWarnings(cor(A_D, A_M, method = "spearman"))
A_p    <- cor.test(A_D, A_M)$p.value

cat("\n===== astrocyte depth ~ astrocyte density (interior quadrats), S3R1 =====\n")
cat(sprintf("n quadrats = %d | n astrocytes = %d | mean depth = %.1f\n",
            length(A_D), nrow(ast), A_mean))
cat(sprintf("Pearson r = %+.3f | Spearman rho = %+.3f | p = %.3g\n", A_r, A_rho, A_p))

## ---------------------------------------------------------------
## 5. grid matrices for the spatial maps (interior grids)
## ---------------------------------------------------------------
nI <- length(cxs); nJ <- length(cys)
Dmap <- matrix(NA, nI, nJ); MDmap <- matrix(NA, nI, nJ)
mdall <- tapply(d$depth, d$gk, mean)
for (gk in names(area)) {
  ij <- as.integer(strsplit(gk, "_")[[1]]); i <- ij[1] + 1; j <- ij[2] + 1
  if (gk %in% names(dens))  Dmap[i, j]  <- dens[gk]
  if (gk %in% names(mdall)) MDmap[i, j] <- mdall[gk]
}
xm <- cxs / 1000; ym <- cys / 1000                 # mm

## ---------------------------------------------------------------
## 6. figure: (A) density map  (B) depth map  (C) astrocyte depth vs density
## ---------------------------------------------------------------
pal <- if (requireNamespace("viridisLite", quietly = TRUE))
  viridisLite::magma(256) else hcl.colors(256, "Inferno")

draw_bar <- function(zlim, pal, lab) {
  par(mar = c(4, 0.6, 3, 5.4))
  zz <- seq(zlim[1], zlim[2], length.out = length(pal))
  image(1, zz, matrix(zz, nrow = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 1.0, mgp = c(3, 0.6, 0), tcl = -0.3)
  mtext(lab, side = 4, line = 3.8, cex = 0.95); box()
}

png(file.path(out_dir, "MERFISH_S3R1_depth_diagnostic.png"),
    width = 2400, height = 900, res = 200)
layout(matrix(1:5, nrow = 1), widths = c(1, 0.32, 1, 0.32, 1.3))

## (A) density
par(mar = c(4, 4, 3, 1))
image(xm, ym, Dmap, col = pal, asp = 1, xlab = "x (mm)", ylab = "y (mm)",
      main = "Local cell density")
mtext("A", side = 3, adj = -0.18, line = 1, font = 2, cex = 1.2)
draw_bar(range(Dmap, na.rm = TRUE), pal, expression(cells~mm^-2))

## (B) mean depth
par(mar = c(4, 4, 3, 1))
image(xm, ym, MDmap, col = pal, asp = 1, xlab = "x (mm)", ylab = "y (mm)",
      main = "Mean per-cell detection depth")
mtext("B", side = 3, adj = -0.18, line = 1, font = 2, cex = 1.2)
draw_bar(range(MDmap, na.rm = TRUE), pal, expression(transcripts~cell^-1))

## (C) astrocyte depth vs astrocyte density, Fig. S12B style
pt_col <- adjustcolor("#4B5A8A", alpha.f = 0.70)
par(mar = c(4, 4.2, 3, 1))
plot(A_D, A_M, pch = 16, cex = 0.5, col = pt_col,
     xlim = c(0, max(A_D)), ylim = c(0, max(A_M) * 1.05),
     xlab = "Local astrocyte density (cells / quadrat)",
     ylab = "Mean per-cell depth (transcripts)", main = "")
abline(h = A_mean, col = "#E64B35", lty = 2, lwd = 2)

## grey title strip, as in S12B
usr <- par("usr"); hh <- 0.09 * (usr[4] - usr[3])
rect(usr[1], usr[4], usr[2], usr[4] + hh, col = "grey85", border = "grey60", xpd = NA)
text(mean(usr[1:2]), usr[4] + hh / 2, "Astrocytes", font = 2, cex = 1.15, xpd = NA)
mtext("C", side = 3, adj = -0.20, line = 1, font = 2, cex = 1.2)

legend("topright", bg = "white", box.col = "grey70", box.lwd = 0.8, cex = 0.9,
       pch = c(16, NA), lty = c(NA, 2), lwd = c(NA, 2),
       col = c(pt_col, "#E64B35"), seg.len = 1.6, y.intersp = 1.15,
       legend = c(sprintf("Quadrat mean (r = %+.2f, n = %d)", A_r, length(A_D)),
                  sprintf("Mean depth = %.0f", A_mean)))
dev.off()
cat("\nWrote MERFISH_S3R1_depth_diagnostic.png\n")

## optional: per-quadrat data behind panel C
write.csv(data.frame(cell_type = "Astrocyte",
                     astrocytes_per_quadrat = A_D,
                     mean_depth = A_M),
          file.path(out_dir, "MERFISH_S3R1_astrocyte_depth_density_perquadrat.csv"),
          row.names = FALSE)
cat("Wrote MERFISH_S3R1_astrocyte_depth_density_perquadrat.csv\n")





