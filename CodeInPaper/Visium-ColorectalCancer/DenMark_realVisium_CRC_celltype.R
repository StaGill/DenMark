## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds the 10x Visium CRC download in Revision-R1/spot-Visium/ (see README.md)
OUTPUT_DIR <- "output"   # Revision-R1/spot-Visium/ below it holds nuclei.csv and spot_celltype.csv from segment_nuclei_settled.ipynb and receives the results
REPO_ROOT  <- "."        # root of the DenMark checkout
STAN_DIR   <- file.path(REPO_ROOT, "CodeInPaper", "Visium-ColorectalCancer")   # holds DenMark-M2-vararea-libsize.stan, not in this repository (see "Missing files" in README.md)
dir.create(file.path(OUTPUT_DIR, "Revision-R1/spot-Visium"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

###############################################################################
##  DenMark on REAL 10x Visium (55um spots) — CELL-TYPE REGION + log(UMI)
##  -------------------------------------------------------------------------
##  Verification that DenMark works at spot level and recovers a biologically
##  sensible density-expression coupling for ONE cell type.
##
##  Cell-type definition (Route A + option B): tumor-epithelial spots are read
##  from spot_celltype.csv, produced in Python by Leiden clustering + marker
##  labelling (cluster is "tumor" if mean tumor score > mean stroma score).
##  DenMark is then run ONLY on those spots:
##     density Y   = nuclei per spot (Cellpose segmentation) within tumor spots
##     expression  = spot gene counts, offset by log(total UMI) (library size)
##     coupling    = Astar[2,1] between tumor-region cell density and expression
##
##  ** Anti-circularity: the genes that DEFINE the region (in Python: EPCAM/KRT/
##     CDX2/CEACAM5) must be DISJOINT from the genes TESTED here (TEST_GENES). **
##
##  Expected biology (CRC): dense tumor glands -> higher proliferation share
##  (MKI67) -> positive Astar[2,1]. MALAT1 = abundance control, expected ~0.
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(fields); library(jsonlite)
  library(RANN); library(Seurat)
  library(cmdstanr); library(posterior)
  library(ggplot2); library(ggsci); library(viridis); library(patchwork); library(dplyr)
  library(sf)
})
set.seed(2024)

## ================== EDIT: paths ==========================================
OUTS_DIR  <- file.path(DATA_DIR, "Revision-R1/spot-Visium")
H5_FILE   <- file.path(OUTS_DIR, "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5")
SPATIAL   <- file.path(OUTS_DIR, "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial", "spatial")
NUC_FILE  <- file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", "nuclei.csv")
CT_FILE   <- file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", "spot_celltype.csv")         # <- from Python (barcode, celltype)
stan_file <- file.path(STAN_DIR, "DenMark-M2-vararea-libsize.stan")     # <- your compiled Stan model
SAMPLE    <- "CRC_tumor"
## =========================================================================

NUC_X <- "x"; NUC_Y <- "y"; NUC_SCALE <- 1

## genes TESTED for coupling (must be DISJOINT from the Python DEFINE genes)
TEST_GENES     <- c("MKI67","MYC","PPIA")
## DEFINE genes = the markers the Python step actually scored on. Kept here only
## as an anti-circularity tripwire; CEACAM5 IS a CRC defining marker, so it must
## appear here (and therefore must NOT be tested). Long-term: have Python export
## its define-gene list to a file and read it in, so the two can't silently drift.
DEFINE_GENES   <- c("EPCAM","KRT8","KRT18","KRT20","CDX2","CEACAM5")
stopifnot(length(intersect(DEFINE_GENES, TEST_GENES)) == 0)

basis <- 8; boundfactor <- 1.5

theme_nature <- function(b=7) theme_classic(base_size=b, base_family="sans") +
  theme(axis.line=element_line(linewidth=.3), axis.ticks=element_line(linewidth=.3),
        axis.text=element_text(colour="black"), plot.title=element_text(size=b+1, hjust=0),
        legend.key.height=unit(3,"mm"), legend.key.width=unit(2.5,"mm"),
        legend.title=element_text(size=b), legend.text=element_text(size=b-1))
theme_map <- function(b=7) theme_nature(b)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text=element_blank(),axis.title=element_blank())
NPG <- pal_npg("nrc")(10)


## ---- 1. Load matrix, spot coords, scalefactors; MEASURED pixel scale -------
mat <- Read10X_h5(H5_FILE)                            # genes x barcodes (in-tissue spots)

sf_json <- fromJSON(file.path(SPATIAL, "scalefactors_json.json"))
pos_all <- fread(file.path(SPATIAL, "tissue_positions.csv")); setnames(pos_all, tolower(names(pos_all)))

## scale from the measured 100um spot pitch (more reliable than drawn spot_diameter)
xy_all   <- as.matrix(pos_all[, .(pxl_col_in_fullres, pxl_row_in_fullres)])
pitch_px <- median(RANN::nn2(xy_all, xy_all, k=2)$nn.dists[,2])
UM_PER_PX <- 100 / pitch_px                           # ~0.275 um/px
SPOT_R_UM <- 27.5                                     # physical 55um spot radius
SPOT_AREA <- pi*SPOT_R_UM^2
cat(sprintf("pitch %.1f px -> %.4f um/px (spot radius %.0f px)\n",
            pitch_px, UM_PER_PX, SPOT_R_UM/UM_PER_PX))

pos <- pos_all[in_tissue==1]
common <- intersect(pos$barcode, colnames(mat))
pos <- pos[match(common, pos$barcode)]; mat <- mat[, common]
pos[, x_um := pxl_col_in_fullres * UM_PER_PX]
pos[, y_um := pxl_row_in_fullres * UM_PER_PX]
pos[, UMI := Matrix::colSums(mat)]                    # library size = sequencing depth
cat(sprintf("in-tissue spots: %d | median UMI/spot: %.0f\n", nrow(pos), median(pos$UMI)))


## ---- 2. Tumor region: read the Python cluster annotation ------------------
##  spot_celltype.csv (barcode, celltype) from Leiden clustering + marker labels.
ct  <- fread(CT_FILE); setnames(ct, tolower(names(ct)))
pos <- merge(pos, ct[, .(barcode, celltype)], by = "barcode", all.x = TRUE, sort = FALSE)
pos <- pos[match(common, pos$barcode)]                # keep matrix/order alignment
pos[, tumor := celltype == "Tumor"]
stopifnot(sum(pos$tumor, na.rm=TRUE) > 0)             # tripwire: annotation actually merged
cat(sprintf("tumor spots (from Python annotation): %d / %d\n",
            sum(pos$tumor, na.rm=TRUE), nrow(pos)))


## ---- 3. Nuclei -> cells per spot (Y) --------------------------------------
nuc <- fread(NUC_FILE)
nuc <- data.frame(x_um = nuc[[NUC_X]]*NUC_SCALE*UM_PER_PX,
                  y_um = nuc[[NUC_Y]]*NUC_SCALE*UM_PER_PX)
nn   <- RANN::nn2(as.matrix(pos[, .(x_um, y_um)]), as.matrix(nuc[, c("x_um","y_um")]), k=1)
sidx <- ifelse(nn$nn.dists[,1] <= SPOT_R_UM, nn$nn.idx[,1], NA_integer_)
pos[, Y := tabulate(sidx[!is.na(sidx)], nbins = nrow(pos))]
cat(sprintf("nuclei: %d | landing in a spot: %.0f%% | mean cells/spot: %.1f\n",
            nrow(nuc), 100*mean(!is.na(sidx)), mean(pos$Y)))

## per-spot counts for the TEST genes
present <- intersect(TEST_GENES, rownames(mat))
cat("test genes present:", paste(present, collapse=", "), "\n")
miss <- setdiff(TEST_GENES, present); if(length(miss)) cat("  (missing:", paste(miss,collapse=", "),")\n")
for (g in present) pos[[g]] <- as.integer(mat[g, ])

## RESTRICT to tumor spots with >=1 nucleus  <- the cell-type step
spots <- as.data.frame(pos[tumor == TRUE & Y > 0,
                           c("x_um","y_um","Y","UMI", present), with=FALSE])
cat("tumor spots used (Y>0):", nrow(spots), "\n")


## ---- 4. Aperture-7 (7-spot) pooling  ->  super-cells `agg` -----------------
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
PITCH_UM <- 100
ap <- pool_aperture7_xy(spots$x_um, spots$y_um, PITCH_UM)
spots$sid<-as.integer(factor(ap$sid)); spots$ex<-ap$ex; spots$ey<-ap$ey
## fold orphan pools (<3 spots) into the nearest large pool; area offset (below)
## neutralizes the resulting pool-size differences, so this does not bias density.
sz<-table(spots$sid); large<-as.integer(names(sz[sz>=3])); sm<-!spots$sid%in%large
if(any(sm)&&length(large)>0){ lc<-spots[!sm,]|>group_by(sid)|>summarise(ex=mean(ex),ey=mean(ey))
  j<-apply(fields::rdist(as.matrix(spots[sm,c("ex","ey")]),as.matrix(lc[,c("ex","ey")])),1,which.min)
  spots$sid[sm]<-lc$sid[j] }
spots$sidf<-factor(spots$sid); lev<-levels(spots$sidf)
agg<-data.frame(sid=lev,
  x =as.numeric(tapply(spots$x_um,spots$sidf,mean)[lev]),
  y =as.numeric(tapply(spots$y_um,spots$sidf,mean)[lev]),
  Y =as.integer(tapply(spots$Y,  spots$sidf,sum)[lev]),
  UMI=as.numeric(tapply(spots$UMI,spots$sidf,sum)[lev]),
  n_spots=as.integer(table(spots$sidf)[lev]))
for(g in present) agg[[g]]<-as.integer(tapply(spots[[g]],spots$sidf,sum)[lev])
agg$dens <- agg$Y/agg$n_spots
x0<-mean(range(agg$x));y0<-mean(range(agg$y));s<-max(diff(range(agg$x))/2,diff(range(agg$y))/2)
agg$sx<-(agg$x-x0)/s; agg$sy<-(agg$y-y0)/s
cat(sprintf("super-cells: %d | pool sizes: %s\n", nrow(agg),
            paste(as.integer(table(table(spots$sid))), collapse="/")))


## ---- visualize the 7-spot super-cell tessellation (buffered-union blobs) ---
##  Faithful depiction: each pool's footprint = union of buffered discs around
##  its member spots. No convex-hull bridging across gaps, no overlapping blobs,
##  and 1-2-spot pools still get a shape. (This replaces the misleading hull plot.)
pal <- rep(c("#F4B3B9","#B7D3B0","#C9BBD9","#9EC3E0","#F6D6A8"), 40)
buf <- 0.62 * PITCH_UM

sp <- spots
sp$y_um <- -sp$y_um                       # flip y here instead of scale_y_reverse()

pool_poly <- sp %>%
  st_as_sf(coords = c("x_um","y_um")) %>%
  group_by(sid) %>% summarise(do_union = TRUE, .groups = "drop") %>%
  st_buffer(dist = buf) %>% st_buffer(dist = -buf*0.35)
pool_poly$fill <- pal[as.integer(factor(pool_poly$sid)) %% length(pal) + 1]

fig_pool <- ggplot() +
  geom_sf(data = pool_poly, aes(fill = fill), colour = "grey40", linewidth = .2, alpha = .6) +
  geom_point(data = sp, aes(x_um, y_um), size = .5, colour = "grey20") +
  scale_fill_identity() +
  coord_sf() + theme_void() +             # coord_sf, NOT coord_equal + scale_y_reverse
  ggtitle(sprintf("Aperture-7 pooling: %d super-cells from %d tumor spots",
                  nrow(agg), nrow(spots)))
print(fig_pool)
ggsave(file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", sprintf("Fig_pooling_%s.png", SAMPLE)), fig_pool,
       width = 150, height = 140, units = "mm", dpi = 320)


## ---- 5. DenMark fit -> Astar[2,1] (mean + 95% CI) per TEST gene ------------
#cmdstanr::rebuild_cmdstan(cores = 4)

mod   <- cmdstan_model(stan_file)
m1 <- m2 <- basis; mstar <- m1*m2
S <- unname(as.matrix(expand.grid(S2=1:m1,S1=1:m2)[,2:1])); storage.mode(S) <- "integer"

fit_gene <- function(agg, gene){
  coords <- as.matrix(agg[,c("sx","sy")]); n <- nrow(coords)
  L  <- boundfactor*c(max(abs(coords[,1])),max(abs(coords[,2])))
  la <- log(agg$n_spots*SPOT_AREA); la <- la - mean(la)
  log_lib <- log(agg$UMI + 0.05); log_lib <- log_lib - mean(log_lib)
  input <- list(d=2, L=L, mstar=mstar, n=n, coords=coords,
                log_grid_area=la, log_lib=log_lib,
                y1=as.integer(agg$Y), y2=as.integer(agg[[gene]]),
                mrange=max(fields::rdist(coords))/(2*2.75), indices=S, is_centerted_PHI=0L)
  init <- function() list(Astar=matrix(c(1,0.3,0,sqrt(1-0.09)),2,2), sigma=c(1,1),
                          beta0=0, beta1=-1, ell=c(1,1), betab1=rep(0,mstar), betab2=rep(0,mstar))
  fit <- mod$sample(data=input, chains=4, parallel_chains=4,
                    iter_warmup=1500, iter_sampling=1500,
                    adapt_delta=0.97, init=init, refresh=600, show_messages=TRUE)
  A21 <- as.numeric(as_draws_matrix(fit$draws("Astar[2,1]")))
  sm  <- fit$summary("Astar[2,1]")
  data.frame(gene=gene,
             A21_mean=mean(A21), A21_q2.5=unname(quantile(A21,.025)),
             A21_q50=median(A21), A21_q97.5=unname(quantile(A21,.975)),
             rhat=sm$rhat, ess_bulk=sm$ess_bulk)
}
res <- do.call(rbind, lapply(present, function(g){ cat("fitting", g, "\n"); fit_gene(agg, g) }))
res <- cbind(sample=SAMPLE, res)
cat("\n===== tumor-region Astar[2,1] (log-UMI offset),", SAMPLE, "=====\n")
print(res[,c("gene","A21_mean","A21_q2.5","A21_q50","A21_q97.5","rhat","ess_bulk")],
      row.names=FALSE, digits=3)
write.csv(res, file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", sprintf("realVisium_%s_Astar21.csv", SAMPLE)), row.names=FALSE)


## ---- 6. Figures -----------------------------------------------------------
figA <- ggplot(pos, aes(x_um, y_um, colour=tumor)) + geom_point(size=.6) +
  scale_colour_manual(values=c(`FALSE`="grey85",`TRUE`=NPG[1]), name="tumor spot",
                      na.value="grey85") +
  scale_y_reverse() + coord_equal() + theme_map() + ggtitle("Tumor region (Leiden + markers)")
figB <- ggplot(spots, aes(x_um, y_um, colour=Y)) + geom_point(size=.8) +
  scale_colour_viridis_c(option="magma", name="cells/spot", trans="log1p") +
  scale_y_reverse() + coord_equal() + theme_map() + ggtitle("Cell density in tumor spots (Y)")
figC <- ggplot(agg, aes(x, y, colour=dens)) + geom_point(size=2.4) +
  scale_colour_viridis_c(option="magma", name="cells/spot", trans="log1p") +
  scale_y_reverse() + coord_equal() + theme_map() + ggtitle("Pooled density (7-spot super-cell)")

resF <- res; resF$gene <- factor(resF$gene, levels=rev(present))
figD <- ggplot(resF, aes(A21_mean, gene)) +
  geom_vline(xintercept=0, colour="grey60", linetype=2, linewidth=.4) +
  geom_errorbarh(aes(xmin=A21_q2.5, xmax=A21_q97.5), height=.22, linewidth=.5, colour=NPG[1]) +
  geom_point(size=2, colour=NPG[1]) +
  labs(x=expression(A[21]*" (tumor density-expression, log(UMI) offset): mean + 95% CI"),
       y=NULL, title=sprintf("%s: Astar[2,1]", SAMPLE)) + theme_nature()

fig <- (figA | figB) / (figC | figD) +
  plot_annotation(title=sprintf("DenMark on tumor compartment — %s", SAMPLE),
                  tag_levels="a") &
  theme(plot.tag=element_text(size=10, face="bold"))
ggsave(file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", sprintf("Fig_realVisium_%s.pdf", SAMPLE)), fig, width=190, height=150, units="mm")
ggsave(file.path(OUTPUT_DIR, "Revision-R1/spot-Visium", sprintf("Fig_realVisium_%s.png", SAMPLE)), fig, width=190, height=150, units="mm", dpi=320)

cat("\nDone.", SAMPLE, "- tumor-region Astar[2,1] CSV + figure written.\n")
###############################################################################



###############################################################################
##  Re-plot the DenMark posterior for rho (Astar[2,1]) for MKI67 and MYC
##  from the saved results CSV, without refitting.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(ggsci) })

## ---- paths ----
OUTS_DIR <- file.path(OUTPUT_DIR, "Revision-R1/spot-Visium")
SAMPLE   <- "CRC_tumor"
RES_FILE <- file.path(OUTS_DIR, sprintf("realVisium_%s_Astar21.csv", SAMPLE))
## (the fitting script wrote the CSV to the R working directory; if it is not in
##  OUTS_DIR, set RES_FILE to that location instead)

GENES <- c("MKI67", "MYC")                       # top-to-bottom order in the plot
NPG   <- pal_npg("nrc")(10)
COLS  <- c(MKI67 = NPG[3], MYC = NPG[1])          # NPG green (#00A087), NPG red (#E64B35)

## ---- load ----
res <- read.csv(RES_FILE)
res <- res[res$gene %in% GENES, ]
stopifnot(nrow(res) == length(GENES))
res$gene  <- factor(res$gene, levels = rev(GENES))
res$label <- sprintf("%.2f", res$A21_mean)
res$strip <- "DenMark fitted posterior for \u03c1"
print(res[, c("gene", "A21_mean", "A21_q2.5", "A21_q97.5", "rhat", "ess_bulk")],
      row.names = FALSE, digits = 3)

## ---- plot ----
fig_rho <- ggplot(res, aes(A21_mean, gene, colour = gene)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = A21_q2.5, xmax = A21_q97.5), height = 0.18, linewidth = 0.6) +
  geom_point(size = 2.2) +
  geom_text(aes(label = label), nudge_y = 0.22, size = 2.6, colour = "black") +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  facet_wrap(~ strip) +
  labs(x = expression("Nucleus density" * "\u2013" * "expression correlation (" * rho * ")"),
       y = NULL) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line        = element_line(linewidth = 0.3),
        axis.ticks       = element_line(linewidth = 0.3),
        axis.text        = element_text(colour = "black"),
        axis.text.y      = element_text(face = "bold", size = 8),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text       = element_text(face = "bold", size = 8.5),
        plot.margin      = margin(4, 6, 4, 4))
print(fig_rho)

ggsave(file.path(OUTS_DIR, sprintf("Fig_rho_%s_MKI67_MYC.pdf", SAMPLE)), fig_rho,
       width = 70, height = 75, units = "mm")
ggsave(file.path(OUTS_DIR, sprintf("Fig_rho_%s_MKI67_MYC.png", SAMPLE)), fig_rho,
       width = 70, height = 75, units = "mm", dpi = 400)



