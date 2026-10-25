## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds the Xenium breast tables in Revision-R1/Visium-Breast-Cancer/xenium-breast-cancer/rep1/ (see README.md)
OUTPUT_DIR <- "output"   # CSV tables and figures are written here
REPO_ROOT  <- "."        # root of the DenMark checkout (holds inst/stan/)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

###############################################################################
##  DenMark pseudo-Visium  ->  Astar[2,1] ONLY
##  -------------------------------------------------------------------------
##  GOAL: report the density-expression coupling rho = Astar[2,1] from
##  SPOT-AGGREGATED (pseudo-Visium) data, and compare it to the single-cell
##  (Xenium) result gene by gene.
##
##  PROVENANCE: single cells (fixed to one cell type as an oracle) are degraded
##  to a 100um pseudo-Visium grid, then pooled with the SAME aperture-7 (7-spot)
##  design used previously (Section 3). The pooled super-cells (`agg`) are the
##  ONLY input to the fit, so every Astar[2,1] below is a pseudo-Visium quantity.
##
##  This version reports ONE thing: the posterior of Astar[2,1] (mean + 95% CI)
##  per gene. No cor(w1,w2), no cor(z1,z2), no sign-flipping, no theta, no raw
##  Spearman cross-checks. Astar[2,1] is the coregionalization correlation
##  (R[2,1] = Astar[2,1], since Astar[1,1] = 1), i.e. the model coupling itself.
##
##  Run once per CELLTYPE. Colour: NPG / magma (Nature style).
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(fields)
  library(cmdstanr); library(posterior)
  library(ggplot2); library(ggsci); library(viridis); library(patchwork); library(dplyr)
})
set.seed(2024)

## ================== CHOOSE ONE CELL TYPE ==================
CELLTYPE <- "Invasive_Tumor"     # "Invasive_Tumor" | "DCIS" | "CD8+_T_Cells"
## =========================================================

## ---- paths (EDIT for your machine) ------------------------------------------
DIR       <- file.path(DATA_DIR, "Revision-R1/Visium-Breast-Cancer/xenium-breast-cancer/rep1")                       # cell_spatial.csv etc.
stan_file <- file.path(REPO_ROOT, "inst", "stan", "DenMark-vararea.stan")
## existing single-cell (Xenium) couplings for the comparison overlay (Section 6):
##   a CSV with columns:  celltype, gene, sc_coupling   (sc_coupling = single-cell Astar[2,1])
SC_FILE   <- file.path(DIR, "xenium_singlecell_coupling.csv")  # set NA to skip

UM_PER_PX <- as.numeric(readLines(file.path(DIR,"pixel-size-raw.txt")))     # 0.5
SPOT_R_UM <- as.numeric(readLines(file.path(DIR,"radius-raw.txt")))*UM_PER_PX # 28 um
SPOT_PITCH<- 100                  # pseudo-Visium pitch (um), standard Visium
basis <- 8; boundfactor <- 1.5; SPOT_AREA <- pi*SPOT_R_UM^2

## label groups + gene panel
label_set <- switch(CELLTYPE,
  "Invasive_Tumor" = "Invasive_Tumor",
  "DCIS"           = c("DCIS_1","DCIS_2"),
  "CD8+_T_Cells"   = "CD8+_T_Cells")
if (CELLTYPE=="CD8+_T_Cells") {
  GENES  <- c("CD8A","GZMB","PTPRC","MALAT1")
} else {                              # tumour (DCIS / Invasive)
  GENES  <- c("ERBB2","EPCAM","KRT8","MALAT1")
}

theme_nature <- function(b=7) theme_classic(base_size=b, base_family="sans") +
  theme(axis.line=element_line(linewidth=.3), axis.ticks=element_line(linewidth=.3),
        axis.text=element_text(colour="black"), plot.title=element_text(size=b+1, hjust=0),
        legend.key.height=unit(3,"mm"), legend.key.width=unit(2.5,"mm"),
        legend.title=element_text(size=b), legend.text=element_text(size=b-1))
theme_map <- function(b=7) theme_nature(b)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text=element_blank(),axis.title=element_blank())
NPG <- pal_npg("nrc")(10)


## ---- 1. Load cells (+labels) and cell x gene; SELECT the cell type ----------
cells <- fread(file.path(DIR,"cell_spatial.csv"))
setnames(cells, c("cell_row","x","y","celltype")[seq_len(ncol(cells))])
cbg <- fread(file.path(DIR,"cell_by_gene_expression.csv"))
gene_cols <- setdiff(names(cbg), names(cbg)[1])
present <- intersect(GENES, gene_cols)
cat("cell type:", CELLTYPE, "| genes found:", paste(present, collapse=", "), "\n")
missing <- setdiff(GENES, present); if(length(missing)) cat("  (missing:", paste(missing,collapse=", "),")\n")

sel <- which(cells$celltype %in% label_set)
cat(sprintf("selected %d %s cells (of %d total)\n", length(sel), CELLTYPE, nrow(cells)))
cs  <- cells[sel, ]
cs$x_um <- cs$x
cs$y_um <- cs$y
gmat <- as.matrix(cbg[sel, ..present]); storage.mode(gmat) <- "double"


## ---- 2. Build pseudo-Visium spot grid; Y (cells/spot) and expr/spot ---------
xr <- range(cs$x_um); yr <- range(cs$y_um)
dy <- SPOT_PITCH*sqrt(3)/2
rows <- seq(yr[1], yr[2], by=dy)
grid <- do.call(rbind, lapply(seq_along(rows), function(i){
  off <- if(i%%2==0) SPOT_PITCH/2 else 0
  xs <- seq(xr[1]+off, xr[2], by=SPOT_PITCH); data.frame(x_um=xs, y_um=rows[i], row=i-1) }))
grid$spot_id <- seq_len(nrow(grid))

g_xy <- as.matrix(grid[,c("x_um","y_um")]); c_xy <- as.matrix(cs[,c("x_um","y_um")])
nn <- apply(c_xy, 1, function(p){ d<-sqrt((g_xy[,1]-p[1])^2+(g_xy[,2]-p[2])^2)
  j<-which.min(d); if(d[j]<=SPOT_R_UM) j else NA_integer_ })
cs$spot_id <- nn
cat(sprintf("cells landing in a pseudo-spot: %.0f%%  (Visium gap coverage ~25-30%% is normal)\n",
            100*mean(!is.na(nn))))

inb <- !is.na(cs$spot_id)
Ytab <- tapply(rep(1L,sum(inb)), cs$spot_id[inb], sum)
spots <- data.frame(spot_id=as.integer(names(Ytab)), Y=as.integer(Ytab))
for(g in present){ mv<-tapply(gmat[inb,g], cs$spot_id[inb], sum)
  spots[[g]] <- as.integer(mv[match(spots$spot_id, as.integer(names(mv)))]); spots[[g]][is.na(spots[[g]])]<-0L }
spots <- merge(spots, grid[,c("spot_id","x_um","y_um")], by="spot_id")
spots <- spots[spots$Y>0, ]
cat("pseudo-spots with cells:", nrow(spots), "\n")



# check the number of genes----
cbg <- fread(file.path(DIR,"cell_by_gene_expression.csv"))
gene_cols <- setdiff(names(cbg), names(cbg)[1])   # every column except the cell-id column

length(gene_cols)        # total genes in the dataset

a<-read.csv(file.path(DIR,"cell_by_gene_expression.csv"))
colnames(a)


gvec <- cbg[[gene_name]]                           # pick the column by gene name


# which genes pass the preprocessing ----
expr <- as.matrix(cbg[, ..gene_cols]); storage.mode(expr) <- "double"   # cells x genes, ROI

## criterion 1: variance across cells in the top 95% (drop the bottom 5% by variance)
v      <- apply(expr, 2, var)
keep_v <- v >= quantile(v, 0.05, na.rm = TRUE)

## criterion 2: non-zero expression in >= 5% of locations (cells)
keep_f <- colMeans(expr > 0) >= 0.05

## genes kept (= J)
genes_kept <- gene_cols[keep_v & keep_f]

length(genes_kept)   # J = number of genes passing
genes_kept           # the gene names that pass


## ---- 3. Aperture-7 (7-spot) pooling  ->  super-cells `agg` -------------------
##  This is the 7-spot design carried over unchanged; `agg` is the pseudo-Visium
##  input to DenMark, so Astar[2,1] below is a pseudo-Visium quantity.
pool_aperture7_xy <- function(x_um, y_um, pitch){
  ax<-x_um/pitch; ay<-y_um/(pitch*sqrt(3)/2); r<-round(ay); q<-round(ax-r/2)
  ex<-q+r/2; ey<-(sqrt(3)/2)*r
  v1<-c(2.5,sqrt(3)/2); v2<-c(0.5,1.5*sqrt(3)); B<-cbind(v1,v2); P<-cbind(ex,ey)
  anchor<-P[which.min(rowSums(sweep(P,2,colMeans(P))^2)),]
  AB<-t(solve(B,t(sweep(P,2,anchor))))
  G<-as.matrix(expand.grid(A=(floor(min(AB[,1]))-1):(ceiling(max(AB[,1]))+1),
                           B=(floor(min(AB[,2]))-1):(ceiling(max(AB[,2]))+1)))
  list(sid=apply(fields::rdist(P, sweep(G%*%t(B),2,anchor,"+")),1,which.min), ex=ex, ey=ey)
}
ap <- pool_aperture7_xy(spots$x_um, spots$y_um, SPOT_PITCH)
spots$sid<-as.integer(factor(ap$sid)); spots$ex<-ap$ex; spots$ey<-ap$ey
sz<-table(spots$sid); large<-as.integer(names(sz[sz>=3])); sm<-!spots$sid%in%large
if(any(sm)&&length(large)>0){ lc<-spots[!sm,]|>group_by(sid)|>summarise(ex=mean(ex),ey=mean(ey))
  j<-apply(fields::rdist(as.matrix(spots[sm,c("ex","ey")]),as.matrix(lc[,c("ex","ey")])),1,which.min)
  spots$sid[sm]<-lc$sid[j] }
spots$sidf<-factor(spots$sid); lev<-levels(spots$sidf)
agg<-data.frame(sid=lev,
  x=as.numeric(tapply(spots$x_um,spots$sidf,mean)[lev]),
  y=as.numeric(tapply(spots$y_um,spots$sidf,mean)[lev]),
  Y=as.integer(tapply(spots$Y,spots$sidf,sum)[lev]),
  n_spots=as.integer(table(spots$sidf)[lev]))
for(g in present) agg[[g]]<-as.integer(tapply(spots[[g]],spots$sidf,sum)[lev])
agg$dens<-agg$Y/agg$n_spots
x0<-mean(range(agg$x));y0<-mean(range(agg$y));s<-max(diff(range(agg$x))/2,diff(range(agg$y))/2)
agg$sx<-(agg$x-x0)/s; agg$sy<-(agg$y-y0)/s
cat(sprintf("super-cells: %d | pool sizes: %s\n", nrow(agg),
            paste(as.integer(table(table(spots$sid))), collapse="/")))





## ---- 4. DenMark fit  ->  Astar[2,1] posterior (mean + 95% CI) per gene -------
mod   <- cmdstan_model(stan_file)
m1 <- m2 <- basis; mstar <- m1*m2
S <- unname(as.matrix(expand.grid(S2=1:m1,S1=1:m2)[,2:1])); storage.mode(S) <- "integer"

fit_gene <- function(agg, gene){
  coords <- as.matrix(agg[,c("sx","sy")]); n <- nrow(coords)
  L  <- boundfactor*c(max(abs(coords[,1])),max(abs(coords[,2])))
  la <- log(agg$n_spots*SPOT_AREA); la <- la - mean(la)      # centered per-grid area offset
  input <- list(d=2, L=L, mstar=mstar, n=n, coords=coords, log_grid_area=la,
                y1=as.integer(agg$Y), y2=as.integer(agg[[gene]]),
                mrange=max(fields::rdist(coords))/(2*2.75), indices=S, is_centerted_PHI=0L)
  init <- function() list(Astar=matrix(c(1,0.3,0,sqrt(1-0.09)),2,2), sigma=c(1,1),
                          beta0=0, beta1=-1, ell=c(1,1), betab1=rep(0,mstar), betab2=rep(0,mstar))
  fit <- mod$sample(data=input, chains=4, parallel_chains=4,
                    iter_warmup=1500, iter_sampling=1500,
                    adapt_delta=0.97, init=init, refresh=600, show_messages=TRUE)

  ## the coupling: rho = Astar[2,1]  (= R[2,1], the coregionalization correlation)
  A21 <- as.numeric(as_draws_matrix(fit$draws("Astar[2,1]")))
  sm  <- fit$summary("Astar[2,1]")                          # rhat / ess for convergence only

  data.frame(gene       = gene,
             A21_mean   = mean(A21),
             A21_q2.5   = unname(quantile(A21, .025)),
             A21_q50    = median(A21),
             A21_q97.5  = unname(quantile(A21, .975)),
             rhat       = sm$rhat,
             ess_bulk   = sm$ess_bulk)
}

res <- do.call(rbind, lapply(present, function(g){ cat("fitting", g, "\n"); fit_gene(agg, g) }))
res <- cbind(celltype = CELLTYPE, res)

cat("\n===== pseudo-Visium  Astar[2,1]  (mean + 95% CI),", CELLTYPE, "=====\n")
print(res[,c("gene","A21_mean","A21_q2.5","A21_q50","A21_q97.5","rhat","ess_bulk")],
      row.names=FALSE, digits=3)
write.csv(res, file.path(OUTPUT_DIR, sprintf("pseudoVisium_%s_Astar21.csv", gsub("[^A-Za-z0-9]","",CELLTYPE))),
          row.names=FALSE)


## ---- 5. Figures -------------------------------------------------------------
##  Provenance (a-c: the 7-spot pseudo-Visium construction) + result
##  (d: Astar[2,1] per gene, posterior mean + 95% CI).



figA <- ggplot(cs, aes(x_um/1000, y_um/1000)) +       # positive y, no minus
  geom_point(size = .15, colour = "#4E7B9C", alpha = .5) +
  geom_point(data = grid, aes(x_um/1000, y_um/1000),
             shape = 21, size = 1.8, stroke = 0.5, colour = "grey20", fill = "grey85") +
  coord_equal() +                                            # <- no scale_y_reverse()
  labs(x = "x (mm)", y = "y (mm)")+
       #title = "Invasive Tumour Cells: Generated Visium spots"
  theme_nature() +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15))
figA



figB <- ggplot(spots, aes(x_um, -y_um, colour = Y)) + geom_point(size = 1.2) +
  scale_colour_viridis_c(option = "magma", name = "cells/spot", trans = "log1p",
                         breaks = c(10, 20, 40),
                         guide = guide_colourbar(barheight = unit(30, "mm"),
                                                 barwidth  = unit(3, "mm"))) +
  scale_y_reverse() + coord_equal() +
  labs(x = expression("x ("*mu*"m)"),
       y = expression("y ("*mu*"m)")) +
  theme_nature() +
  theme(legend.text  = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 12)) +
  ggtitle("")

figB



figB <- ggplot(spots, aes(x_um/1000, y_um/1000, colour = Y)) + geom_point(size = 1.2) +
  scale_colour_viridis_c(option = "magma", name = "cells/spot", trans = "log1p",
                         breaks = c(10, 20, 40),
                         guide = guide_colourbar(barheight = unit(30, "mm"),
                                                 barwidth  = unit(3, "mm"))) +
  coord_equal() +                                           # <- no scale_y_reverse()
  labs(x = "x (mm)", y = "y (mm)") +
  theme_nature() +
  theme(legend.text  = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 12)) +
  ggtitle("")

figB



figC <- ggplot(agg, aes(x/1000, y/1000, colour = dens)) + geom_point(size = 1.2) +
  scale_colour_viridis_c(option = "magma", name = "cells/spot", trans = "log1p") +
  coord_equal() +                                          # <- no scale_y_reverse()
  labs(x = "x (mm)", y = "y (mm)") +
  theme_nature() +                                         # <- see note below
  theme(legend.text  = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 12)) +
  ggtitle("Pooled density (7-spot super-cell)")

figC




ggplot(agg, aes(x, -y, colour = dens)) + geom_point(size = 2) +
  scale_colour_viridis_c(option = "magma", name = "cells/grid", trans = "log1p",
                         breaks = c(10, 20, 40),
                         guide = guide_colourbar(barheight = unit(30, "mm"),
                                                 barwidth  = unit(3, "mm"))) +
  scale_y_reverse() + coord_equal() +
  labs(x = expression("x ("*mu*"m)"),
       y = expression("y ("*mu*"m)")) +
  theme_nature() +
  theme(legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 10)) +
  ggtitle("")



## Astar[2,1] per gene: posterior mean (point) + 95% CI (bar)
resF <- res; resF$gene <- factor(resF$gene, levels=rev(present))
figD <- ggplot(resF, aes(A21_mean, gene)) +
  geom_vline(xintercept=0, colour="grey60", linetype=2, linewidth=.4) +
  geom_errorbarh(aes(xmin=A21_q2.5, xmax=A21_q97.5), height=.22, linewidth=.5, colour=NPG[1]) +
  geom_point(size=2, colour=NPG[1]) +
  labs(x=expression(A[21]*" (density–expression coupling): mean + 95% CI"),
       y=NULL, title=sprintf("%s: pseudo-Visium Astar[2,1]", CELLTYPE)) +
  theme_nature()

fig <- (figA | figB) / (figC | figD) +
  plot_annotation(title=sprintf("Pseudo-Visium Astar[2,1] — %s", CELLTYPE), tag_levels="a") &
  theme(plot.tag=element_text(size=10, face="bold"))
ggsave(file.path(OUTPUT_DIR, sprintf("Fig_pseudoVisium_%s.pdf", gsub("[^A-Za-z0-9]","",CELLTYPE))), fig,
       width=180, height=140, units="mm")
ggsave(file.path(OUTPUT_DIR, sprintf("Fig_pseudoVisium_%s.png", gsub("[^A-Za-z0-9]","",CELLTYPE))), fig,
       width=180, height=140, units="mm", dpi=320)


## ---- 6. Compare to single-cell (Xenium) Astar[2,1] --------------------------
##  SC_FILE columns:  celltype, gene, sc_coupling   (sc_coupling = single-cell Astar[2,1])
##  Recovery = pseudo-Visium Astar[2,1] (mean + 95% CI) plotted against single-cell,
##  with the y = x line: points on the line mean the 7-spot aggregation recovers
##  the single-cell coupling.
if (!is.na(SC_FILE) && file.exists(SC_FILE)) {
  sc  <- read.csv(SC_FILE)
  cmp <- merge(res, sc, by=c("celltype","gene"))
  if (nrow(cmp) > 0) {
    figV <- ggplot(cmp, aes(sc_coupling, A21_mean)) +
      geom_abline(slope=1, intercept=0, linetype=2, colour="grey60") +
      geom_errorbar(aes(ymin=A21_q2.5, ymax=A21_q97.5), width=0, linewidth=.4, colour=NPG[1]) +
      geom_point(size=2, colour=NPG[1]) +
      ggrepel::geom_text_repel(aes(label=gene), size=2.4) +
      labs(x="single-cell (Xenium) Astar[2,1]",
           y="pseudo-Visium Astar[2,1] (mean + 95% CI)",
           title=sprintf("Recovery: pseudo-Visium vs single-cell — %s", CELLTYPE)) +
      coord_equal() + theme_nature()
    ggsave(file.path(OUTPUT_DIR, sprintf("Fig_recovery_%s.pdf", gsub("[^A-Za-z0-9]","",CELLTYPE))), figV,
           width=90, height=90, units="mm")
    ggsave(file.path(OUTPUT_DIR, sprintf("Fig_recovery_%s.png", gsub("[^A-Za-z0-9]","",CELLTYPE))), figV,
           width=90, height=90, units="mm", dpi=320)
    cat("\n===== recovery: pseudo-Visium vs single-cell Astar[2,1],", CELLTYPE, "=====\n")
    print(cmp[,c("gene","sc_coupling","A21_mean","A21_q2.5","A21_q97.5")],
          row.names=FALSE, digits=3)
    write.csv(cmp[,c("celltype","gene","sc_coupling","A21_mean","A21_q2.5","A21_q50","A21_q97.5")],
              file.path(OUTPUT_DIR, sprintf("recovery_%s_Astar21.csv", gsub("[^A-Za-z0-9]","",CELLTYPE))), row.names=FALSE)
  } else {
    cat("\n[no gene overlap between pseudo-Visium result and SC_FILE]\n")
  }
} else {
  cat("\n[single-cell comparison skipped: set SC_FILE to a CSV with columns\n",
      " celltype, gene, sc_coupling (single-cell Astar[2,1]) to produce the recovery figure]\n")
}

cat("\nDone.", CELLTYPE, "— Astar[2,1] CSV + figures written.\n")
###############################################################################



## ---- X. Spot occupancy per cell type (56-um pseudo-Visium spots) -----------
LABELS <- list("CD8+_T_Cells"   = "CD8+_T_Cells",
               "Invasive_Tumor" = "Invasive_Tumor",
               "DCIS"           = c("DCIS_1", "DCIS_2"))

## Same rule as Section 2 (a cell belongs to the nearest spot centre if within SPOT_R_UM),
## vectorised so it runs on all cells in seconds. Spots are 100 um apart and 56 um wide,
## so at most one spot can contain a given cell.
spot_key <- function(x, y, xr, yr, pitch = SPOT_PITCH, R = SPOT_R_UM) {
  dy     <- pitch * sqrt(3) / 2
  nrow_g <- floor((yr[2] - yr[1]) / dy) + 1
  key    <- rep(NA_character_, length(x))
  r0     <- floor((y - yr[1]) / dy)
  for (r in list(r0, r0 + 1)) {                           # the two nearest grid rows
    off  <- ifelse(r %% 2 == 1, pitch / 2, 0)             # odd rows shifted by half a pitch (as in Section 2)
    cmax <- floor((xr[2] - xr[1] - off) / pitch)
    cc   <- pmin(pmax(round((x - xr[1] - off) / pitch), 0), cmax)
    ok   <- r >= 0 & r < nrow_g & cmax >= 0
    d    <- sqrt((x - (xr[1] + off + cc * pitch))^2 + (y - (yr[1] + r * dy))^2)
    hit  <- ok & d <= R & is.na(key)
    key[hit] <- paste(r[hit], cc[hit], sep = "_")
  }
  key
}


## common tissue-wide grid -> denominator = spots holding >= 1 cell of ANY type
key_all  <- spot_key(cells$x, cells$y, range(cells$x), range(cells$y))
n_tissue <- length(unique(na.omit(key_all)))

occ <- do.call(rbind, lapply(names(LABELS), function(ct) {
  i  <- which(cells$celltype %in% LABELS[[ct]])
  xi <- cells$x[i]; yi <- cells$y[i]
  k_own <- spot_key(xi, yi, range(xi), range(yi))   # per-type grid, exactly what Section 2 feeds DenMark
  k_com <- key_all[i]                               # common grid, for comparing across types
  cps   <- as.integer(table(k_own[!is.na(k_own)]))  # cells per occupied spot
  n_com <- length(unique(na.omit(k_com)))
  data.frame(celltype            = ct,
             n_cells             = length(i),
             pct_cells_in_spots  = 100 * mean(!is.na(k_own)),
             n_spots_fit         = length(cps),
             n_spots             = n_com,
             pct_of_tissue_spots = 100 * n_com / n_tissue,
             cells_per_spot_med  = median(cps),
             cells_per_spot_min  = min(cps),
             cells_per_spot_max  = max(cps))
}))

cat("56-um pseudo-spots containing >= 1 cell (any type):", n_tissue, "\n")
print(occ, row.names = FALSE, digits = 3)
write.csv(occ, file.path(OUTPUT_DIR, "pseudoVisium_spot_occupancy.csv"), row.names = FALSE)

## sanity check against the Section-2 run for the current CELLTYPE
if (exists("spots")) cat(sprintf("check: Section 2 gave %d spots for %s; here n_spots_fit = %d\n",
                                 nrow(spots), CELLTYPE, occ$n_spots_fit[occ$celltype == CELLTYPE]))





