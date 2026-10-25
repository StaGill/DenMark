## ------------------------------------------------------------------------- ##
## 5. Visualisation                                                          ##
## ------------------------------------------------------------------------- ##

.theme_denmark <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 1),
                   legend.key.height = ggplot2::unit(0.9, "lines"))
}

# boundary polygon -> data.frame of rings for geom_path
.boundary_df <- function(boundary) {
  if (is.null(boundary)) return(NULL)
  cc <- as.data.frame(sf::st_coordinates(boundary))
  grp_cols <- intersect(c("L1", "L2", "L3"), names(cc))
  cc$ring <- interaction(cc[grp_cols], drop = TRUE)
  cc
}

.add_boundary <- function(p, boundary, colour = "black") {
  bd <- .boundary_df(boundary)
  if (is.null(bd)) return(p)
  p + ggplot2::geom_path(data = bd, ggplot2::aes(x = .data$X, y = .data$Y,
                                                 group = .data$ring),
                         inherit.aes = FALSE, colour = colour, linewidth = 0.5)
}

.tile_plot <- function(df, value, title, limits = NULL, legend = NULL,
                       xlab = "x", ylab = "y", diverging = FALSE, boundary = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                        fill = .data[[value]])) +
    ggplot2::geom_tile(colour = "grey30", linewidth = 0.15) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = xlab, y = ylab, fill = legend) +
    .theme_denmark()
  p <- .add_boundary(p, boundary, colour = "gray60")
  if (diverging) {
    lim <- if (is.null(limits)) max(abs(df[[value]])) * c(-1, 1) else limits
    p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                      limits = lim)
  } else {
    p + ggplot2::scale_fill_viridis_c(option = "magma", limits = limits)
  }
}

.combine <- function(plots, ncol) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = ncol)
  } else {
    message("Install 'patchwork' to combine the panels; returning a list of plots.")
    plots
  }
}

#' Gridded cell counts (Y) and summed expression (M)
#'
#' Tile maps of the cell count and the summed expression of each grid cell.
#' With a tissue boundary, only the grid cells inside the tissue are drawn,
#' together with the tissue outline.
#'
#' @param obj A `denmark` fit, or the result of [buildgridpp()].
#' @param xlab,ylab Axis labels.
#' @return The two panels combined with patchwork, or a list of two ggplot
#'   objects when patchwork is not installed.
#' @export
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' plot_denmark_grid(buildgridpp(grid_size = 10, dataset = tumor, gene = "geneA",
#'                               boundary = tissue))
plot_denmark_grid <- function(obj, xlab = "x", ylab = "y") {
  g <- obj$grid
  gene <- obj$gene
  boundary <- if (inherits(obj, "denmark")) obj$info$boundary else obj$boundary
  if (!is.null(g$type)) g <- g[g$type != "exterior", , drop = FALSE]
  .combine(list(
    .tile_plot(g, "Y", "Gridded cell counts", legend = "cells", xlab = xlab,
               ylab = ylab, boundary = boundary),
    .tile_plot(g, "M", paste0("Gridded expression: ", gene), legend = "counts",
               xlab = xlab, ylab = ylab, boundary = boundary)), ncol = 2)
}

#' Tissue boundary and grid classification
#'
#' Left: all cells (gray), the cells of the analyzed type (colored) and the
#' alpha-shape tissue boundary. Right (when `grid` is given): the grid cells
#' colored by type (interior, boundary, exterior) with the boundary outline.
#'
#' @param cells Data frame with the columns `x` and `y` (all cells).
#' @param boundary Polygon from [denmark_alpha_shape()].
#' @param grid Optional output of `buildgridpp(..., boundary = boundary)`.
#' @param highlight Optional logical vector, one element per row of `cells`:
#'   `TRUE` for the cells of the analyzed type.
#' @param xlab,ylab Axis labels.
#' @param point_size Size of the cell points.
#' @return A ggplot object; with `grid`, the two panels combined with
#'   patchwork, or a list of two ggplot objects when patchwork is not
#'   installed.
#' @export
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' g <- buildgridpp(grid_size = 10, dataset = tumor, gene = "geneA",
#'                  boundary = tissue)
#' plot_denmark_boundary(toy, tissue, grid = g,
#'                       highlight = toy$celltype == "Tumor")
plot_denmark_boundary <- function(cells, boundary, grid = NULL, highlight = NULL,
                                  xlab = "x", ylab = "y", point_size = 0.15) {
  cells <- as.data.frame(cells)
  cells$group <- if (is.null(highlight)) "cells" else
    ifelse(highlight, "analyzed cell type", "other cells")
  cells <- cells[order(cells$group == "analyzed cell type"), ]
  p1 <- ggplot2::ggplot(cells, ggplot2::aes(.data$x, .data$y, colour = .data$group)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(values = c(`analyzed cell type` = "#B2182B",
                                            `other cells` = "grey75", cells = "grey40")) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = sprintf("Alpha-shape tissue boundary (alpha = %.3g)",
                                  attr(boundary, "alpha") %||% NA),
                  x = xlab, y = ylab, colour = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) +
    .theme_denmark() + ggplot2::theme(legend.position = "bottom")
  p1 <- .add_boundary(p1, boundary)
  if (is.null(grid)) return(p1)
  g <- grid$grid
  dx <- grid$step[["dx"]]; dy <- grid$step[["dy"]]
  g$type <- factor(g$type, c("interior", "boundary", "exterior"))
  p2 <- ggplot2::ggplot(g) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - dx / 2, xmax = .data$x + dx / 2,
                                    ymin = .data$y - dy / 2, ymax = .data$y + dy / 2,
                                    fill = .data$type),
                       colour = "grey40", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = c(interior = "#9ECAE1", boundary = "#FDAE6B",
                                          exterior = "white"), drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = sprintf("Grid cells: %d interior, %d boundary, %d exterior",
                                  sum(g$type == "interior"), sum(g$type == "boundary"),
                                  sum(g$type == "exterior")),
                  x = xlab, y = ylab, fill = NULL) +
    .theme_denmark() + ggplot2::theme(legend.position = "bottom")
  p2 <- .add_boundary(p2, boundary)
  .combine(list(p1, p2), ncol = 2)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Latent-field decomposition maps (posterior means)
#'
#' Row 1 (cell density): `beta0 + a11 w1` and `a11 w1`, on a shared color
#' scale. Row 2 (expression): `beta1 + a21 w1 + a22 w2`, `a21 w1` and
#' `a22 w2`, on a shared color scale.
#'
#' @param obj A `denmark` object fitted with `keep_fit = TRUE`.
#' @param fields Optional output of [denmark_fields()]; computed from `obj`
#'   when `NULL`.
#' @param xlab,ylab Axis labels.
#' @return A patchwork figure, or a list of five ggplot objects when
#'   patchwork is not installed.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' plot_denmark_fields(fit)
plot_denmark_fields <- function(obj, fields = NULL, xlab = "x", ylab = "y") {
  if (is.null(fields)) fields <- denmark_fields(obj)
  lim1 <- range(c(fields$eta1, fields$a11w1))
  lim2 <- range(c(fields$eta2, fields$a21w1, fields$a22w2))
  ttl <- list(eta1  = quote(beta[0] + a[11] * w[1]),
              a11w1 = quote(a[11] * w[1]),
              eta2  = quote(beta[1] + a[21] * w[1] + a[22] * w[2]),
              a21w1 = quote(a[21] * w[1]),
              a22w2 = quote(a[22] * w[2]))
  leg <- c(eta1 = "density", a11w1 = "density", eta2 = "expression",
           a21w1 = "expression", a22w2 = "expression")
  p <- lapply(names(ttl), function(v)
    .tile_plot(fields, v, as.expression(ttl[[v]]),
               limits = if (v %in% c("eta1", "a11w1")) lim1 else lim2,
               legend = leg[[v]], xlab = xlab, ylab = ylab,
               boundary = obj$info$boundary))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    row1 <- (p[[1]] | p[[2]] | patchwork::plot_spacer()) +
      patchwork::plot_layout(guides = "collect")
    row2 <- (p[[3]] | p[[4]] | p[[5]]) + patchwork::plot_layout(guides = "collect")
    row1 / row2
  } else p
}

#' Spatial correlation decay curves with posterior bands
#'
#' Posterior mean and band of the spatial correlation of the density and
#' expression fields, with the threshold that defines the practical range as
#' a dashed line.
#'
#' @param obj A `denmark` object from [DenMark()].
#' @param decay Optional output of [denmark_decay()]; computed from `obj`
#'   when `NULL`.
#' @param xlab Label of the distance axis.
#' @return A ggplot object.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' fit <- DenMark(tumor, gene = "geneA", boundary = tissue,
#'                grid_size = 10, basis = 10,
#'                chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                seed = 2024, quiet = TRUE)
#' plot_denmark_decay(fit)
plot_denmark_decay <- function(obj, decay = NULL, xlab = "Distance") {
  if (is.null(decay)) decay <- denmark_decay(obj)
  cv <- decay$curves
  cv$field <- factor(cv$field, c("density", "expression"),
                     c("r1: cell density field", "r2: expression field"))
  ggplot2::ggplot(cv, ggplot2::aes(x = .data$distance)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), linewidth = 1) +
    ggplot2::geom_hline(yintercept = decay$threshold, linetype = 2, colour = "grey40") +
    ggplot2::facet_wrap(~ field) +
    ggplot2::labs(x = xlab, y = "Posterior spatial correlation") +
    .theme_denmark()
}

#' Forest plot of rho across genes from denmark_scan()
#'
#' Posterior mean and credible interval of rho for each gene, colored by the
#' DCG call.
#'
#' @param scan Data frame returned by [denmark_scan()].
#' @param top Optional number of genes to show, those with the largest
#'   absolute posterior mean of rho.
#' @return A ggplot object.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' scan <- denmark_scan(tumor, genes = c("geneA", "geneC"), boundary = tissue,
#'                      grid_size = 10,
#'                      basis = 10, chains = 2, iter_warmup = 300,
#'                      iter_sampling = 300, seed = 2024, verbose = FALSE)
#' plot_denmark_scan(scan)
plot_denmark_scan <- function(scan, top = NULL) {
  s <- scan[!is.na(scan$rho_mean), , drop = FALSE]
  if (!is.null(top) && nrow(s) > top)
    s <- s[order(-abs(s$rho_mean)), ][seq_len(top), ]
  s$gene <- factor(s$gene, levels = s$gene[order(s$rho_mean)])
  s$class <- factor(s$class, c("PCG", "NCG", "not DCG"))
  ggplot2::ggplot(s, ggplot2::aes(x = .data$rho_mean, y = .data$gene,
                                  colour = .data$class)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$rho_lower, xmax = .data$rho_upper),
                            height = 0.25) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(PCG = "#B2182B", NCG = "#2166AC",
                                            `not DCG` = "grey50"), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::labs(x = expression("Density-expression correlation " * rho *
                                 " (posterior mean, credible interval)"),
                  y = NULL, colour = NULL) +
    .theme_denmark()
}
