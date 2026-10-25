# Internal environment: caches compiled models so that a gene scan compiles
# only once.
.denmark_env <- new.env(parent = emptyenv())
.denmark_env$models <- list()
