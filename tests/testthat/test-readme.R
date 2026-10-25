# Run the R examples of README.md as written, on the bundled toy data (the
# R analog of PGAlign's tests/test_docs.py). README.md is not part of the
# built package, so this test runs from a source checkout only.

readme_r_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\s*```\\s*[rR]\\s*$", lines)
  ends <- grep("^\\s*```\\s*$", lines)
  lapply(starts, function(s) {
    e <- ends[ends > s][1]
    if (is.na(e)) stop("README.md: code block at line ", s, " is not closed")
    lines[s + seq_len(e - s - 1L)]
  })
}

# Blocks that read the reader's own files or install software are skipped.
# Without CmdStan, the run stops at the first block that fits a model.
run_readme_blocks <- function(blocks, cmdstan_ok) {
  installs <- "\\b(install\\.packages|install_github|install_cmdstan|check_cmdstan_toolchain|set_cmdstan_path)\\("
  fits <- "\\b(DenMark|DenMarkmain|denmark_scan|denmark_model)\\("
  env <- new.env(parent = globalenv())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ran <- 0L
  for (i in seq_along(blocks)) {
    code <- paste(blocks[[i]], collapse = "\n")
    if (grepl('read.csv("', code, fixed = TRUE) || grepl(installs, code, perl = TRUE))
      next
    if (!cmdstan_ok && grepl(fits, code, perl = TRUE))
      return(list(ran = ran, stopped = TRUE))
    tryCatch({
      for (expr in parse(text = code, keep.source = FALSE)) {
        res <- withVisible(eval(expr, env))
        if (res$visible) utils::capture.output(print(res$value))
      }
    }, error = function(e) {
      stop(sprintf("README.md R block %d failed: %s\n%s", i, conditionMessage(e), code),
           call. = FALSE)
    })
    ran <- ran + 1L
  }
  list(ran = ran, stopped = FALSE)
}

test_that("README R examples run as written", {
  readme <- test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not found (not part of the built package)")
  blocks <- readme_r_blocks(readme)
  skip_if(length(blocks) == 0, "README.md has no R code blocks")
  res <- run_readme_blocks(blocks, cmdstan_ok = cmdstan_available())
  if (res$stopped)
    skip(sprintf("CmdStan not installed: ran %d README block(s), stopped at the first model fit",
                 res$ran))
  expect_gt(res$ran, 0)
})
