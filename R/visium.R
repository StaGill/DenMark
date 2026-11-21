## ------------------------------------------------------------------------- ##
## Spot-level Visium: seven-spot (aperture-7) super-grids                    ##
## ------------------------------------------------------------------------- ##

# Euclidean distances between the rows of `a` and the rows of `b`, summed
# coordinate by coordinate as fields::rdist() does.
.denmark_rdist <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  d2 <- 0
  for (k in seq_len(ncol(a))) d2 <- d2 + outer(a[, k], b[, k], "-")^2
  sqrt(d2)
}

# Mingchi Xu's aperture-7 pooling from the spot-level Visium analysis of the
# paper (CodeInPaper/Visium-ColorectalCancer), kept verbatim except that
# fields::rdist() is replaced by .denmark_rdist(). Spots are snapped to the
# hexagonal lattice of pitch `pitch`, and each is assigned to the nearest
# point of the coarser lattice spanned by v1 and v2, whose cells hold seven
# spots: a centre and its six neighbours.
.denmark_pool_aperture7_xy <- function(x_um, y_um, pitch){
  ax<-x_um/pitch; ay<-y_um/(pitch*sqrt(3)/2); r<-round(ay); q<-round(ax-r/2)
  ex<-q+r/2; ey<-(sqrt(3)/2)*r
  v1<-c(2.5,sqrt(3)/2); v2<-c(0.5,1.5*sqrt(3)); B<-cbind(v1,v2); P<-cbind(ex,ey)
  anchor<-P[which.min(rowSums(sweep(P,2,colMeans(P))^2)),]
  AB<-t(solve(B,t(sweep(P,2,anchor))))
  G<-as.matrix(expand.grid(A=(floor(min(AB[,1]))-1):(ceiling(max(AB[,1]))+1),
                           B=(floor(min(AB[,2]))-1):(ceiling(max(AB[,2]))+1)))
  list(sid=apply(.denmark_rdist(P, sweep(G%*%t(B),2,anchor,"+")),1,which.min), ex=ex, ey=ey)
}


#' Pool Visium spots into seven-spot super-grids
#'
#' A Visium section holds thousands of spots, and the cost of a DenMark fit
#' grows with the number of units. As in the spot-level analysis of the
#' paper, adjacent spots are pooled into hexagonal super-grids of seven spots
#' (a spot and its six neighbours, an aperture-7 tiling of the spot lattice).
#' Pools with fewer than `min_spots` spots, which arise at the edge of the
#' region, are folded into the nearest pool with at least `min_spots` spots.
#' Counts are summed over the spots of each pool, and each pool gets the area
#' `n_spots * spot_area`, which enters the model as its area offset.
#'
#' The cell (or nucleus) count per spot comes from outside DenMark, for
#' example from nuclei segmented in the H&E image and assigned to the spot
#' that contains them.
#'
#' @param spots Data frame with one row per spot: the coordinates, the number
#'   of focal cells or nuclei per spot, optionally the total count per spot,
#'   and the gene count columns.
#' @param genes Gene count columns to sum (default: every numeric column other
#'   than the coordinates, `cells` and `umi`).
#' @param x,y Coordinate columns, in the same unit as `pitch` (micrometres for
#'   the default).
#' @param cells Column with the number of cells or nuclei per spot.
#' @param umi Column with the total count per spot (library size), summed for
#'   `offset = "umi"`; `NULL` to skip.
#' @param pitch Centre-to-centre distance between neighbouring spots (100
#'   micrometres on Visium).
#' @param min_spots Pools with fewer spots are folded into the nearest larger
#'   pool.
#' @param spot_area Area of one spot (default: a 55-micrometre spot). With a
#'   common value it only shifts the intercept of the density model.
#' @return A data frame with one row per super-grid and the columns `unit`,
#'   `x`, `y` (mean spot coordinates), the summed `cells` and `umi` columns
#'   under their input names, `n_spots`, `area` and the summed gene counts.
#'   The attribute `"spot_unit"` gives the unit of each input spot. Pass the
#'   result to [DenMark()] with `units = TRUE`.
#' @export
#' @examples
#' # a small hexagonal spot lattice, 100 um apart
#' xy <- expand.grid(i = 0:13, j = 0:11)
#' spots <- data.frame(x = 100 * (xy$i + 0.5 * (xy$j %% 2)),
#'                     y = 100 * sqrt(3) / 2 * xy$j)
#' set.seed(1)
#' spots$Y <- rpois(nrow(spots), 8)
#' spots$UMI <- rpois(nrow(spots), 5000)
#' spots$geneA <- rpois(nrow(spots), 3)
#' units <- denmark_pool_visium(spots)
#' table(units$n_spots)
#' head(units)
#' \dontrun{
#' fit <- DenMark(units, gene = "geneA", units = TRUE,
#'                offset = "umi", umi = "UMI", basis = 8)
#' }
denmark_pool_visium <- function(spots, genes = NULL, x = "x", y = "y", cells = "Y",
                                umi = "UMI", pitch = 100, min_spots = 3,
                                spot_area = pi * 27.5^2) {
  spots <- as.data.frame(spots)
  need <- c(x, y, cells, umi)
  miss <- setdiff(need, names(spots))
  if (length(miss)) stop("Missing column(s) in `spots`: ", paste(miss, collapse = ", "), ".")
  if (is.null(genes)) {
    num <- names(spots)[vapply(spots, is.numeric, logical(1))]
    genes <- setdiff(num, need)
  }
  miss <- setdiff(genes, names(spots))
  if (length(miss)) stop("Gene column(s) not in `spots`: ", paste(miss, collapse = ", "), ".")
  for (v in c(need, genes)) {
    if (!is.numeric(spots[[v]])) stop("Column '", v, "' is not numeric.")
    if (anyNA(spots[[v]])) stop("Column '", v, "' has missing values.")
  }
  for (v in c(cells, genes))
    if (any(spots[[v]] < 0) || any(abs(spots[[v]] - round(spots[[v]])) > 1e-8))
      stop("Column '", v, "' must hold non-negative integer counts.")
  if (nrow(spots) < 2) stop("At least two spots are needed.")
  if (pitch <= 0) stop("`pitch` must be positive.")

  ap <- .denmark_pool_aperture7_xy(spots[[x]], spots[[y]], pitch)
  sid <- as.integer(factor(ap$sid)); ex <- ap$ex; ey <- ap$ey
  ## fold orphan pools (< min_spots spots) into the nearest large pool
  sz <- table(sid); large <- as.integer(names(sz[sz >= min_spots])); sm <- !sid %in% large
  if (any(sm) && length(large) > 0) {
    lc <- stats::aggregate(data.frame(ex = ex[!sm], ey = ey[!sm]), list(sid = sid[!sm]), mean)
    j <- apply(.denmark_rdist(cbind(ex[sm], ey[sm]), as.matrix(lc[, c("ex", "ey")])), 1, which.min)
    sid[sm] <- lc$sid[j]
  }

  f <- factor(sid); lev <- levels(f)
  units <- data.frame(unit = as.integer(lev),
                      x = as.numeric(tapply(spots[[x]], f, mean)[lev]),
                      y = as.numeric(tapply(spots[[y]], f, mean)[lev]))
  units[[cells]] <- as.integer(tapply(spots[[cells]], f, sum)[lev])
  if (!is.null(umi)) units[[umi]] <- as.numeric(tapply(spots[[umi]], f, sum)[lev])
  units$n_spots <- as.integer(table(f)[lev])
  units$area <- units$n_spots * spot_area
  for (g in genes) units[[g]] <- as.integer(tapply(spots[[g]], f, sum)[lev])
  attr(units, "spot_unit") <- sid
  units
}
