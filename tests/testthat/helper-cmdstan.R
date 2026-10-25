cmdstan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
}

skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
    skip("CmdStan is not installed")
}

skip_if_no_toy_data <- function() {
  if (!nzchar(system.file("extdata", "denmark_toy.csv.gz", package = "DenMark")))
    skip("Toy data not found; run data-raw/make_toy_data.R")
}
