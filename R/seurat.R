## ------------------------------------------------------------------------- ##
## 7. Seurat objects: input and output                                       ##
## ------------------------------------------------------------------------- ##

#' Cell table from a Seurat object
#'
#' Extracts what DenMark needs from a Seurat object with imaging-based spatial
#' data, such as Xenium, MERSCOPE or CosMx data loaded with `LoadXenium()`,
#' `LoadVizgen()` or `LoadNanostring()`: the cell centroids, a cell-type label
#' and the raw counts, one row per cell. Compute the tissue boundary from all
#' rows with [denmark_alpha_shape()], then analyze the cells of one type.
#'
#' @param object A Seurat object with a field of view holding cell centroids.
#' @param celltype Name of the column of the cell metadata that holds the cell
#'   types, for example the result of a clustering or an annotation.
#' @param genes Genes to include (default: all features of the assay). Each
#'   gene becomes a dense column, so select genes for whole-transcriptome
#'   panels.
#' @param assay Assay holding the raw counts (default: the default assay).
#' @param fov Field of view holding the cell centroids (default: the first).
#' @return Data frame with the columns `x`, `y` (in the units of the field of
#'   view, micrometers for Xenium), `celltype` and one count column per gene,
#'   with the cell names as row names.
#' @seealso [denmark_annotate()] to write the results of [denmark_scan()] back
#'   into the object.
#' @export
#' @examplesIf requireNamespace("SeuratObject", quietly = TRUE)
#' # a small Seurat object built from the toy data
#' toy <- denmark_toy_data()
#' counts <- Matrix::Matrix(t(as.matrix(toy[, c("geneA", "geneB", "geneC")])),
#'                          sparse = TRUE)
#' colnames(counts) <- paste0("cell", seq_len(nrow(toy)))
#' obj <- SeuratObject::CreateSeuratObject(counts = counts)
#' obj$celltype <- toy$celltype
#' obj[["fov"]] <- SeuratObject::CreateFOV(
#'   SeuratObject::CreateCentroids(data.frame(x = toy$x, y = toy$y, cell = colnames(counts))),
#'   type = "centroids", assay = "RNA")
#'
#' cells <- denmark_from_seurat(obj, celltype = "celltype")
#' head(cells)
denmark_from_seurat <- function(object, celltype, genes = NULL, assay = NULL, fov = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("denmark_from_seurat() needs the 'SeuratObject' package: install.packages('SeuratObject')")
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  fovs <- SeuratObject::Images(object)
  if (!length(fovs))
    stop("The object has no field of view with cell centroids (spot-based data such as ",
         "Visium are not supported yet).")
  fov <- fov %||% fovs[1]
  xy <- SeuratObject::GetTissueCoordinates(object[[fov]], which = "centroids")
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  genes <- genes %||% rownames(counts)
  missing <- setdiff(genes, rownames(counts))
  if (length(missing))
    stop("Not features of assay '", assay, "': ", paste(utils::head(missing, 5), collapse = ", "))
  meta <- object[[]]
  if (!celltype %in% names(meta))
    stop("'", celltype, "' is not a column of the cell metadata.")
  cells <- intersect(xy$cell, colnames(counts))
  if (!length(cells)) stop("No cell of field of view '", fov, "' is in assay '", assay, "'.")
  xy <- xy[match(cells, xy$cell), ]
  out <- data.frame(x = xy$x, y = xy$y, celltype = as.character(meta[cells, celltype]),
                    row.names = cells, check.names = FALSE)
  expr <- as.matrix(counts[genes, cells, drop = FALSE])
  cbind(out, as.data.frame(t(expr), check.names = FALSE))
}

#' Write DenMark scan results into a Seurat object
#'
#' Adds, for every feature of the assay, the posterior mean of rho, its
#' credible interval and the call from [denmark_scan()] to the feature
#' metadata (`object[[assay]][[]]`). Features that were not scanned get `NA`.
#'
#' @param object A Seurat object.
#' @param scan Data frame returned by [denmark_scan()].
#' @param assay Assay whose feature metadata receive the results (default: the
#'   default assay).
#' @param prefix Prefix of the new columns: `rho`, `rho_lower`, `rho_upper`
#'   and `class`.
#' @return The Seurat object with the new feature metadata.
#' @seealso [denmark_from_seurat()]
#' @export
#' @examplesIf requireNamespace("SeuratObject", quietly = TRUE)
#' counts <- Matrix::Matrix(matrix(rpois(30, 2), 3, 10,
#'                                 dimnames = list(c("geneA", "geneB", "geneC"), NULL)),
#'                          sparse = TRUE)
#' colnames(counts) <- paste0("cell", 1:10)
#' obj <- SeuratObject::CreateSeuratObject(counts = counts)
#' scan <- data.frame(gene = c("geneA", "geneC"), rho_mean = c(0.6, -0.6),
#'                    rho_lower = c(0.3, -0.8), rho_upper = c(0.8, -0.3),
#'                    class = c("PCG", "NCG"))
#' obj <- denmark_annotate(obj, scan)
#' obj[["RNA"]][[]]
denmark_annotate <- function(object, scan, assay = NULL, prefix = "denmark_") {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("denmark_annotate() needs the 'SeuratObject' package: install.packages('SeuratObject')")
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  features <- rownames(object[[assay]])
  idx <- match(features, scan$gene)
  cols <- c(rho = "rho_mean", rho_lower = "rho_lower", rho_upper = "rho_upper", class = "class")
  md <- as.data.frame(lapply(cols, function(cl) scan[[cl]][idx]), row.names = features)
  names(md) <- paste0(prefix, names(cols))
  object[[assay]] <- SeuratObject::AddMetaData(object[[assay]], md)
  object
}
