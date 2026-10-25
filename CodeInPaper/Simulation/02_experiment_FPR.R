## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # not used: the data are simulated
OUTPUT_DIR <- "output"   # fpr_results.rds and fpr_histogram.png are written here
REPO_ROOT  <- "."        # root of the DenMark checkout (holds inst/stan/ and CodeInPaper/)
dir.create(file.path(OUTPUT_DIR, "Revision-R1/FPR+TPR"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

## ============================================================
## 02_experiment_FPR.R  --  FALSE POSITIVE RATE (Type I error)
##
## Design: 1000 independent marked point patterns, each gene has
##         TRUE rho = 0 (no density correlation).
## Metric: FPR = (# genes whose 95% CI for rho excludes 0) / 1000
##         A well-calibrated method gives FPR ~ 0.05.
## ============================================================
source(file.path(REPO_ROOT, "CodeInPaper", "Simulation", "00_functions.R"))
suppressPackageStartupMessages({ library(parallel); library(ggplot2) })

mod    <- compile_model()
coords <- make_grid(n_per_side = 30, half_width = 0.5)   # 144 cells

## ---------------- SPEED KNOBS (read this!) ----------------
## 1000 MCMC fits is heavy. ALWAYS run a pilot first:
##   set N_SIM <- 20, time it, then scale up.
## To go faster, lower: grid size (make_grid), M (basis, below),
## iter_warmup/iter_sampling, or N_SIM. To be more accurate, raise them.
N_SIM   <- 100
M_BASIS <- 6                              # larger = more faithful field, slower
#N_CORES <- max(1, detectCores() - 1)      # Linux/macOS. Windows: see note below.
## ----------------------------------------------------------

run_i <- function(i) {
  out <- tryCatch(
    run_replicate(
      coords, rho = 0, sim_seed = i, fit_seed = 10000 + i, mod = mod,
      stan_args = list(M = M_BASIS, mrange = 0.3),
      fit_args  = list(chains = 2, iter_warmup = 400, iter_sampling = 400)
    ),
    error = function(e) list(detected = NA, error = conditionMessage(e))
  )
  out$i <- i
  out
}

t0 <- Sys.time()
res_list <- mclapply(seq_len(N_SIM), run_i)
## Windows has no fork(): replace the line above with a future backend, e.g.
##   library(future.apply); plan(multisession, workers = N_CORES)
##   res_list <- future_lapply(seq_len(N_SIM), run_i)   # recompile mod inside run_i
cat("elapsed:", format(Sys.time() - t0), "\n")

res <- do.call(rbind, lapply(res_list, function(r)
  data.frame(i        = r$i,
             rho_mean = r$rho_mean %||% NA,
             lower    = r$rho_lower %||% NA,
             upper    = r$rho_upper %||% NA,
             detected = r$detected,
             rhat     = r$rhat %||% NA,
             n_diverg = r$n_diverg %||% NA)))

## keep only converged fits before scoring
ok  <- with(res, !is.na(detected) & (is.na(rhat) | rhat < 1.05))
fpr <- mean(res$detected[ok])

cat(sprintf("\nValid fits ............ %d / %d\n", sum(ok), N_SIM))
cat(sprintf("Detected (false +) .... %d\n", sum(res$detected[ok])))
cat(sprintf("FALSE POSITIVE RATE ... %d / %d = %.3f   (nominal 0.05)\n",
            sum(res$detected[ok]), sum(ok), fpr))

saveRDS(res, file.path(OUTPUT_DIR, "Revision-R1/FPR+TPR/fpr_results.rds"))

## ---------------- visualise ----------------
p <- ggplot(res[ok, ], aes(rho_mean, fill = detected)) +
  geom_histogram(bins = 50, colour = "white") +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d73027"),
                    name = "0 excluded\nfrom 95% CI") +
  labs(title = sprintf("Null genes (true rho = 0):  FPR = %.3f", fpr),
       x = "posterior mean of rho", y = "count") +
  theme_bw()
print(p)
ggsave(file.path(OUTPUT_DIR, "fpr_histogram.png"), p, width = 7, height = 4, dpi = 150)
