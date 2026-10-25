## ------------------------------------------------------------------------- ##
## 1. Data preparation                                                       ##
## ------------------------------------------------------------------------- ##

#' Quality-control filter for genes
#'
#' Keeps genes whose variance across cells is at or above the `var_quantile`
#' quantile of the gene variances and that are detected (count > 0) in at
#' least a fraction `min_nonzero` of the cells. With the defaults this is the
#' filter used in the DenMark paper: the genes in the lowest 5 percent of
#' variance are dropped, and at least 5 percent of the cells must have a
#' non-zero count.
#'
#' @param expr Cells x genes matrix or data frame of counts (numeric columns
#'   only), with the gene names as column names.
#' @param var_quantile Quantile of the gene variances below which genes are
#'   dropped.
#' @param min_nonzero Minimum fraction of cells with a non-zero count.
#' @return Character vector of the genes kept.
#' @export
#' @examples
#' toy <- denmark_toy_data()
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' denmark_filter_genes(tumor[, c("geneA", "geneB", "geneC")])
denmark_filter_genes <- function(expr, var_quantile = 0.05, min_nonzero = 0.05) {
  mat <- as.matrix(expr)
  storage.mode(mat) <- "double"
  if (is.null(colnames(mat))) stop("`expr` must have gene names as column names.")
  gene_var <- apply(mat, 2, stats::var)
  keep_var <- gene_var >= stats::quantile(gene_var, var_quantile, na.rm = TRUE)
  keep_det <- colMeans(mat > 0) >= min_nonzero
  colnames(mat)[keep_var & keep_det]
}


#' Tissue boundary as an alpha shape of the cell locations
#'
#' The Delaunay triangulation of the cells is computed and every triangle
#' whose circumradius is smaller than `alpha` is kept; the union of the kept
#' triangles is the tissue region. This follows the boundary step of the
#' real-data analyses: compute it from all cells of the section (not only the
#' cell type analyzed), so that the density of one cell type is modeled over
#' the whole tissue rather than over its own convex hull.
#'
#' @param x,y Cell coordinates (all cells of the section).
#' @param alpha Circumradius threshold; the default is `alpha_mult` times the
#'   median Delaunay edge length.
#' @param alpha_mult Multiplier for the default `alpha` (3, as in the
#'   analyses).
#' @param keep `"largest"` keeps the largest connected piece of tissue,
#'   `"all"` keeps every piece.
#' @param fill_holes Fill holes inside the tissue.
#' @return An `sfc` polygon (package sf) with the attribute `"alpha"`.
#' @export
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' attr(tissue, "alpha")
#' sf::st_area(tissue)
denmark_alpha_shape <- function(x, y, alpha = NULL, alpha_mult = 3,
                                keep = c("largest", "all"), fill_holes = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("denmark_alpha_shape() needs the 'sf' package: install.packages('sf')")
  keep <- match.arg(keep)
  xy <- unique(cbind(as.numeric(x), as.numeric(y)))
  if (nrow(xy) < 3) stop("At least 3 distinct cell locations are needed.")
  tri <- sf::st_collection_extract(
    sf::st_triangulate(sf::st_sfc(sf::st_multipoint(xy))), "POLYGON")
  cc <- sf::st_coordinates(tri)                      # 4 rows per triangle
  id <- cc[, "L2"]
  first <- match(unique(id), id)
  A <- cc[first, 1:2]; B <- cc[first + 1, 1:2]; C <- cc[first + 2, 1:2]
  la <- sqrt(rowSums((A - B)^2)); lb <- sqrt(rowSums((B - C)^2))
  lc <- sqrt(rowSums((C - A)^2))
  area2 <- abs((B[, 1] - A[, 1]) * (C[, 2] - A[, 2]) - (C[, 1] - A[, 1]) * (B[, 2] - A[, 2]))
  R <- ifelse(area2 > 0, la * lb * lc / (2 * area2), Inf)   # circumradius
  if (is.null(alpha)) alpha <- alpha_mult * stats::median(c(la, lb, lc))
  kept <- tri[R < alpha]
  if (!length(kept)) stop("No triangle passes the alpha threshold; increase `alpha`.")
  poly <- sf::st_union(kept)
  parts <- sf::st_cast(poly, "POLYGON")
  if (keep == "largest" && length(parts) > 1)
    parts <- parts[which.max(as.numeric(sf::st_area(parts)))]
  if (fill_holes)
    parts <- sf::st_sfc(lapply(parts, function(p) sf::st_polygon(list(p[[1]]))))
  out <- sf::st_union(parts)
  attr(out, "alpha") <- alpha
  out
}

# Area of every grid rectangle that lies inside the tissue polygon, and the
# grid type: "interior" (fully inside), "boundary" (partly) or "exterior".
.denmark_grid_area <- function(grid, step, boundary) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("A tissue boundary needs the 'sf' package: install.packages('sf')")
  dx <- step[[1]]; dy <- step[[2]]; full <- dx * dy
  rects <- sf::st_sfc(lapply(seq_len(nrow(grid)), function(k) {
    x0 <- grid$x[k] - dx / 2; y0 <- grid$y[k] - dy / 2
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + dx, y0), c(x0 + dx, y0 + dy),
                              c(x0, y0 + dy), c(x0, y0))))
  }))
  boundary <- sf::st_set_crs(sf::st_sfc(boundary), NA)
  hit <- sf::st_intersects(rects, boundary, sparse = FALSE)[, 1]
  area <- numeric(nrow(grid))
  if (any(hit)) {
    area[hit] <- vapply(which(hit), function(k)
      sum(as.numeric(sf::st_area(sf::st_intersection(rects[k], boundary)))), 0)
  }
  tol <- 1e-6 * full
  type <- ifelse(area <= tol, "exterior", ifelse(area >= full - tol, "interior", "boundary"))
  area[type == "interior"] <- full
  data.frame(area = area, type = type)
}


#' Discretize cell locations and expression on a regular grid
#'
#' Each cell is assigned to exactly one grid cell (cells on an interior edge
#' go to the upper or right cell; cells on the outer boundary are kept).
#'
#' @param grid_size Number of grid cells per axis, or `c(nx, ny)`.
#' @param dataset Data frame with columns `x`, `y` and the expression
#'   column(s), one row per cell.
#' @param gene Expression column to use (default: the first non-coordinate
#'   column).
#' @param bbox Optional bounding box `c(xmin, xmax, ymin, ymax)`; the default
#'   is the data range, or the extent of `boundary` when one is given. Cells
#'   outside the box are dropped.
#' @param boundary Optional tissue polygon from [denmark_alpha_shape()]. The
#'   grid then spans the tissue, and every grid cell gets its area inside the
#'   tissue (`area`) and a `type`: interior, boundary or exterior.
#' @return A list. The vectors are ordered with x varying fastest, that is
#'   `grid_index = i + (j - 1) * nx`, identical to `as.vector(Ymat)`. The first
#'   three elements keep the order of the original DenMark function, so
#'   `res[[1]]`, `res[[2]]` and `res[[3]]` still work.
#'   \item{loc_mark_mat}{N x 2 matrix with the columns `point` (= Y) and
#'     `mark` (= M).}
#'   \item{Ymat}{nx x ny matrix of cell counts (rows = x, columns = y).}
#'   \item{Mmat}{nx x ny matrix of summed expression.}
#'   \item{grid}{Data frame with the columns `grid_index`, `i`, `j`, `x`, `y`
#'     (centroids), `Y`, `M`, `area` (area inside the tissue; the full grid
#'     cell area without a boundary) and `type`.}
#'   \item{bbox, step, grid_area, grid_size, gene, boundary}{Bounding box,
#'     grid steps `c(dx, dy)`, area of one full grid cell, `c(nx, ny)`, the
#'     gene name and the tissue polygon (`NULL` without one).}
#' @export
#' @examples
#' toy <- denmark_toy_data()
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' g <- buildgridpp(grid_size = 10, dataset = tumor, gene = "geneA")
#' dim(g$Ymat)
#' head(g$grid)
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' gb <- buildgridpp(grid_size = 10, dataset = tumor, gene = "geneA",
#'                   boundary = tissue)
#' table(gb$grid$type)
buildgridpp <- function(grid_size, dataset, gene = NULL, bbox = NULL,
                        boundary = NULL) {
  if (inherits(boundary, "sfc") && length(boundary) > 1)
    stop("`boundary` has ", length(boundary), " features; combine them into one ",
         "(multi)polygon first, e.g. with sf::st_union(boundary).")
  dat <- .denmark_xyz(dataset, gene)
  gene <- attr(dat, "gene")
  dims <- .grid_dims(grid_size); nx <- dims[1]; ny <- dims[2]
  bbox_from_boundary <- is.null(bbox) && !is.null(boundary)
  if (bbox_from_boundary) {
    b <- sf::st_bbox(boundary)
    bbox <- c(b[["xmin"]], b[["xmax"]], b[["ymin"]], b[["ymax"]])
  }
  bb <- .bbox(dat$x, dat$y, bbox)

  inside <- dat$x >= bb["xmin"] & dat$x <= bb["xmax"] &
            dat$y >= bb["ymin"] & dat$y <= bb["ymax"]
  if (bbox_from_boundary && any(!inside))
    message(sprintf("%d cell(s) lie outside the extent of the tissue boundary and are dropped.",
                    sum(!inside)))
  if (any(!inside)) dat <- dat[inside, , drop = FALSE]
  if (nrow(dat) == 0)
    stop(if (bbox_from_boundary)
           "No cells inside the extent of the tissue boundary; are the cells and the boundary in the same units?"
         else "No cells inside the bounding box.")

  dx <- unname(bb["xmax"] - bb["xmin"]) / nx
  dy <- unname(bb["ymax"] - bb["ymin"]) / ny
  if (dx <= 0 || dy <= 0) stop("All cells share the same x or y coordinate.")
  ix <- pmin(pmax(floor((dat$x - bb["xmin"]) / dx) + 1L, 1L), nx)
  iy <- pmin(pmax(floor((dat$y - bb["ymin"]) / dy) + 1L, 1L), ny)
  idx <- ix + (iy - 1L) * nx

  Y <- tabulate(idx, nbins = nx * ny)
  M <- numeric(nx * ny)
  s <- rowsum(dat$gene, idx)
  M[as.integer(rownames(s))] <- s[, 1]
  if (any(abs(M - round(M)) > 1e-8))
    warning("Expression values are not integers; summed counts are rounded for the Poisson model.",
            call. = FALSE)
  M <- round(M)

  gi <- rep(seq_len(nx), times = ny)
  gj <- rep(seq_len(ny), each = nx)
  grid <- data.frame(grid_index = seq_len(nx * ny), i = gi, j = gj,
                     x = unname(bb["xmin"]) + (gi - 0.5) * dx,
                     y = unname(bb["ymin"]) + (gj - 0.5) * dy,
                     Y = Y, M = M, area = dx * dy, type = "interior")
  if (!is.null(boundary)) {
    ga <- .denmark_grid_area(grid, c(dx, dy), boundary)
    grid$area <- ga$area
    grid$type <- ga$type
  }

  loc_mark_mat <- cbind(point = Y, mark = M)
  list(loc_mark_mat = loc_mark_mat,
       Ymat = matrix(Y, nx, ny),
       Mmat = matrix(M, nx, ny),
       grid = grid,
       bbox = bb, step = c(dx = dx, dy = dy), grid_area = dx * dy,
       grid_size = c(nx = nx, ny = ny), gene = gene, boundary = boundary)
}


#' Grid centroids, grid area and number of grids
#'
#' @param grid_size Number of grid cells per axis, or `c(nx, ny)`.
#' @param dataset Data frame with the columns `x` and `y` (if these are
#'   missing, the first two columns are used).
#' @param bbox Optional bounding box `c(xmin, xmax, ymin, ymax)`; the default
#'   is the data range.
#' @return `list(grid_points, grid_area, N, step, bbox)`; the first three
#'   elements keep the order of the original function. `grid_points` are the
#'   true centroids of the grid cells (x varying fastest).
#' @export
#' @examples
#' toy <- denmark_toy_data()
#' geom <- grid_pts_centroid(grid_size = 10, dataset = toy)
#' c(grid_area = geom$grid_area, n_grids = geom$N)
grid_pts_centroid <- function(grid_size, dataset, bbox = NULL) {
  dataset <- as.data.frame(dataset)
  if (!all(c("x", "y") %in% names(dataset))) names(dataset)[1:2] <- c("x", "y")
  dims <- .grid_dims(grid_size); nx <- dims[1]; ny <- dims[2]
  bb <- .bbox(dataset$x, dataset$y, bbox)
  dx <- unname(bb["xmax"] - bb["xmin"]) / nx
  dy <- unname(bb["ymax"] - bb["ymin"]) / ny
  grid_points <- expand.grid(x = unname(bb["xmin"]) + (seq_len(nx) - 0.5) * dx,
                             y = unname(bb["ymin"]) + (seq_len(ny) - 0.5) * dy)
  list(grid_points = grid_points, grid_area = dx * dy, N = nrow(grid_points),
       step = c(dx = dx, dy = dy), bbox = bb)
}


#' HSGP basis-function index matrix
#'
#' All pairs of basis-function indices, `1..m1` x `1..m2`, with the first
#' index varying fastest, as used by the Stan models.
#'
#' @param basis Number of basis functions per axis, or `c(m1, m2)`.
#' @return Integer matrix with `m1 * m2` rows and two columns, without
#'   dimnames.
#' @export
#' @examples
#' denmark_hsgp_indices(3)
#' dim(denmark_hsgp_indices(c(10, 8)))
denmark_hsgp_indices <- function(basis) {
  basis <- as.integer(basis)
  if (length(basis) == 1) basis <- rep(basis, 2)
  S <- as.matrix(expand.grid(S1 = seq_len(basis[1]), S2 = seq_len(basis[2])))
  dimnames(S) <- NULL
  storage.mode(S) <- "integer"
  S
}


#' Build the Stan data list for DenMark
#'
#' The grid is centered at the origin, so that the region of interest becomes
#' `[-L1, L1] x [-L2, L2]`; the grid centroids are the HSGP input coordinates
#' and the HSGP domain is `boundfactor * c(L1, L2)`. With a tissue
#' `boundary`, the grid cells outside the tissue are removed first and the
#' region of interest is that of the retained grid.
#'
#' @param dataset Data frame with columns `x`, `y` and the expression
#'   column(s), one row per cell.
#' @param gene Expression column (default: the first non-coordinate column).
#' @param grid_size Grid cells per axis, or `c(nx, ny)`.
#' @param basis Number of HSGP basis functions per axis, or `c(m1, m2)`.
#' @param boundfactor HSGP boundary factor `c` (>= 1; 1.5 is a good default).
#' @param bbox Optional region of interest `c(xmin, xmax, ymin, ymax)`; the
#'   default is the data range, or the extent of `boundary` when one is given.
#' @param L1,L2 Optional half-ranges. Only used to warn if they disagree with
#'   the data (the old interface required them; they are now computed).
#' @param boundary Optional tissue polygon from [denmark_alpha_shape()]. Grid
#'   cells outside the tissue are removed, and each remaining grid cell enters
#'   with its own area inside the tissue: `log_grid_area` becomes a vector and
#'   the variable-area model `DenMark-vararea.stan` is used.
#' @param min_area_frac Drop boundary grid cells whose area inside the tissue
#'   is below this fraction of a full grid cell (default 0: keep all).
#' @return A list with the elements
#'   \item{stan_data}{Data list for the Stan models: `n`, `y1` (cell counts),
#'     `y2` (summed expression), `d`, `mstar`, `coords` (centered grid
#'     centroids), `log_grid_area` (one value, or one per grid cell with a
#'     boundary), `mrange`, `indices`, `L` and `is_centerted_PHI`.}
#'   \item{grid}{The grid from [buildgridpp()] (without the exterior grid
#'     cells when a boundary is given) with the centered coordinates `cx` and
#'     `cy`.}
#'   \item{info}{Gene, grid geometry, HSGP settings, number of cells, the
#'     count matrices `Ymat` and `Mmat`, and for the tissue boundary
#'     `vararea`, `boundary`, `n_grid_total`, `n_interior` and `n_boundary`.}
#' @export
#' @examples
#' toy <- denmark_toy_data()
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' prep <- denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10)
#' str(prep$stan_data[c("n", "mstar", "L", "log_grid_area", "mrange")])
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' prep_b <- denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10,
#'                           boundary = tissue)
#' prep_b$info[c("n_grid_total", "n_interior", "n_boundary")]
#' summary(exp(prep_b$stan_data$log_grid_area))
denmark_prepare <- function(dataset, gene = NULL, grid_size = 10, basis = 25,
                            boundfactor = 1.5, bbox = NULL, L1 = NULL, L2 = NULL,
                            boundary = NULL, min_area_frac = 0) {
  if (boundfactor < 1) stop("`boundfactor` must be >= 1.")
  g <- buildgridpp(grid_size, dataset, gene = gene, bbox = bbox, boundary = boundary)
  bb <- g$bbox
  vararea <- !is.null(boundary)
  grid <- g$grid
  if (vararea) {
    keep <- grid$type != "exterior" & grid$area >= min_area_frac * g$grid_area
    lost <- sum(grid$Y[!keep])
    if (lost > 0)
      message(sprintf("%d cell(s) fall in grid cells outside the tissue and are not modelled.", lost))
    grid <- grid[keep, , drop = FALSE]
    rownames(grid) <- NULL
  }
  # centre the retained grid at the origin: ROI = [-L1, L1] x [-L2, L2]
  centre <- c(x = mean(range(grid$x)), y = mean(range(grid$y)))
  half <- c(L1 = max(abs(grid$x - centre["x"])) + g$step[["dx"]] / 2,
            L2 = max(abs(grid$y - centre["y"])) + g$step[["dy"]] / 2)
  for (k in c("L1", "L2")) {
    user <- get(k)
    if (!is.null(user) && abs(user - half[[k]]) > 0.01 * half[[k]])
      warning(sprintf("%s = %.4g supplied but the data give %.4g; the data value is used.",
                      k, user, half[[k]]), call. = FALSE)
  }

  coords <- cbind(grid$x - centre["x"], grid$y - centre["y"])
  dimnames(coords) <- NULL
  S <- denmark_hsgp_indices(basis)
  dmax <- sqrt(diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2)

  stan_data <- list(
    n = nrow(coords),
    y1 = as.integer(grid$Y),                 # outcome 1: cell counts
    y2 = as.integer(grid$M),                 # outcome 2: summed expression
    d = 2L,
    mstar = nrow(S),                         # total number of basis functions
    coords = coords,                         # centred grid centroids
    log_grid_area = if (vararea) log(grid$area) else log(g$grid_area),  # area offset
    mrange = dmax / (2 * 2.75),              # scale of the inv-gamma prior on ell
    indices = S,
    L = as.numeric(boundfactor * half),      # HSGP domain half-widths
    is_centerted_PHI = 0L)

  grid$cx <- coords[, 1]; grid$cy <- coords[, 2]
  list(stan_data = stan_data,
       grid = grid,
       info = list(gene = g$gene, grid_size = g$grid_size, step = g$step,
                   grid_area = g$grid_area, bbox = bb, centre = centre,
                   L1 = half[["L1"]], L2 = half[["L2"]], basis = basis,
                   boundfactor = boundfactor, max_dist = dmax,
                   n_cells = sum(grid$Y), Ymat = g$Ymat, Mmat = g$Mmat,
                   vararea = vararea, boundary = boundary,
                   n_grid_total = nrow(g$grid),
                   n_interior = sum(grid$type == "interior"),
                   n_boundary = sum(grid$type == "boundary")))
}
