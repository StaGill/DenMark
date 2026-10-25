## ------------------------------------------------------------------------- ##
## 4. DCG discovery: scan many genes                                         ##
## ------------------------------------------------------------------------- ##

#' Fit DenMark gene by gene and call density-correlated genes (DCG)
#'
#' A gene is a positively density-correlated gene (PCG) when the
#' `cred_level` credible interval of rho lies above 0, a negatively
#' density-correlated gene (NCG) when it lies below 0, and "not DCG"
#' otherwise.
#'
#' @param data Data frame with columns `x`, `y` and one column per gene
#'   (cells in rows).
#' @param genes Genes to scan (default: all non-coordinate columns).
#' @param save_dir Optional folder; one `.rds` file per gene is written, and
#'   genes already present are loaded instead of refitted (resumable scans).
#' @param cred_level Credible level used to call DCGs.
#' @param quiet Suppress the sampler output, passed to [DenMark()].
#' @param verbose Report one message per gene with rho and its credible
#'   interval.
#' @param ... Further arguments passed to [DenMark()] (`grid_size`, `basis`,
#'   ...).
#' @return Data frame with one row per gene, sorted by decreasing posterior
#'   mean of rho: `gene`, `rho_mean`, `rho_median`, `rho_lower`, `rho_upper`,
#'   `prob_positive`, `total_counts`, `max_rhat`, `min_ess_bulk`,
#'   `divergences`, `time_min`, `error` (`NA` unless the fit failed), `class`
#'   (PCG, NCG or not DCG) and `converged` (maximum R-hat below 1.01 and no
#'   divergent transitions). The credible level is stored in
#'   `attr(, "cred_level")`.
#' @export
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)) && requireNamespace("sf", quietly = TRUE)
#' toy <- denmark_toy_data()
#' tissue <- denmark_alpha_shape(toy$x, toy$y)
#' tumor <- subset(toy, celltype == "Tumor", select = -celltype)
#' scan <- denmark_scan(tumor, genes = c("geneA", "geneC"), boundary = tissue,
#'                      grid_size = 10,
#'                      basis = 10, chains = 2, iter_warmup = 300,
#'                      iter_sampling = 300, seed = 2024, verbose = FALSE)
#' scan[, c("gene", "rho_mean", "rho_lower", "rho_upper", "class")]
denmark_scan <- function(data, genes = NULL, save_dir = NULL,
                         cred_level = 0.95, quiet = TRUE, verbose = TRUE, ...) {
  .denmark_removed_args(...)
  data <- as.data.frame(data)
  if (is.null(genes)) genes <- setdiff(names(data), c("x", "y"))
  genes <- intersect(genes, names(data))
  if (!length(genes)) stop("None of the requested genes are columns of `data`.")
  if (!is.null(save_dir)) dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  a <- (1 - cred_level) / 2

  one_gene <- function(g) {
    f_rds <- if (!is.null(save_dir)) file.path(save_dir, paste0("denmark_", make.names(g), ".rds"))
    if (!is.null(f_rds) && file.exists(f_rds)) return(readRDS(f_rds))
    t0 <- Sys.time()
    row <- tryCatch({
      f2 <- DenMark(data[, c("x", "y", g)], gene = g,
                    keep_fit = FALSE, quiet = quiet, ...)
      rho <- denmark_draws(f2)$rho
      r <- data.frame(gene = g,
                      rho_mean = mean(rho), rho_median = stats::median(rho),
                      rho_lower = stats::quantile(rho, a, names = FALSE),
                      rho_upper = stats::quantile(rho, 1 - a, names = FALSE),
                      prob_positive = mean(rho > 0),
                      total_counts = sum(f2$stan_data$y2),
                      max_rhat = f2$diagnostics$max_rhat,
                      min_ess_bulk = f2$diagnostics$min_ess_bulk,
                      divergences = f2$diagnostics$num_divergent,
                      time_min = f2$time, error = NA_character_)
      r
    }, error = function(e) data.frame(gene = g, error = conditionMessage(e)))
    if (!is.null(f_rds) && is.na(row$error[1])) saveRDS(row, f_rds)
    if (verbose) {
      if (is.na(row$error[1])) {
        message(sprintf("[%s] rho = %+.3f [%+.3f, %+.3f]  R-hat %.3f  (%.1f min)",
                        g, row$rho_mean, row$rho_lower, row$rho_upper,
                        row$max_rhat, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      } else message(sprintf("[%s] failed: %s", g, row$error))
    }
    row
  }

  rows <- lapply(genes, one_gene)
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA; r[cols]
  }))
  res$class <- ifelse(is.na(res$rho_lower), NA_character_,
               ifelse(res$rho_lower > 0, "PCG",
               ifelse(res$rho_upper < 0, "NCG", "not DCG")))
  res$converged <- !is.na(res$max_rhat) & res$max_rhat < 1.01 &
                   !is.na(res$divergences) & res$divergences == 0
  res <- res[order(-res$rho_mean), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "cred_level") <- cred_level
  res
}
