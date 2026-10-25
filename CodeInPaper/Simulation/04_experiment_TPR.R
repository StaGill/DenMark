## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # not used: the data are simulated
OUTPUT_DIR <- "output"   # tpr_results.rds and tpr_histogram.png are written here
REPO_ROOT  <- "."        # root of the DenMark checkout (holds inst/stan/ and CodeInPaper/)
dir.create(file.path(OUTPUT_DIR, "Revision-R1/FPR+FDR"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

## ============================================================
## 04_experiment_TPR.R  --  TRUE POSITIVE RATE (power / sensitivity)
##
## Design: 1000 (100) independent marked point patterns, each gene is a TRUE
##         DCG with |rho| = R_EFFECT  (500 (50) positive, 500 (50) negative).
## Metric: TPR = (# detected) / (# valid fits), where a gene is detected
##         when 0 is OUTSIDE its 95% CI (the point-null rule).
##         Every miss is a Type II error (false negative).
##
## ============================================================
source(file.path(REPO_ROOT, "CodeInPaper", "Simulation", "00_functions.R"))
suppressPackageStartupMessages({ library(parallel); library(ggplot2) })

mod    <- compile_model()
coords <- make_grid(n_per_side = 12, half_width = 0.5)   # <- grid is designed to be 12*12, same as FPR

## ---------------- design  -----------------------------
N_SIM    <- 100          # total genes, split 50/50 into +/- rho
R_EFFECT <- 0.5           # true |rho| of every DCG 
M_BASIS  <- 5             # number of basis per direction 
## ------------------------------------------------------

## half positive, half negative rho (check the test is sign-symmetric)
true_rho <- c(rep( R_EFFECT, N_SIM %/% 2),
              rep(-R_EFFECT, N_SIM - N_SIM %/% 2))

run_i <- function(i) {
  cat(sprintf("[%s] fit %d / %d\n", format(Sys.time(), "%H:%M:%S"), i, N_SIM))
  out <- tryCatch(
    run_replicate(
      coords, rho = true_rho[i],
      sim_seed = 200000 + i, fit_seed = 300000 + i, mod = mod,
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
cat("elapsed:", format(Sys.time() - t0), "\n")

## CI bounds are saved 
res <- do.call(rbind, lapply(res_list, function(r)
  data.frame(i         = r$i,
             true_rho  = r$true_rho  %||% NA,
             rho_mean  = r$rho_mean  %||% NA,
             rho_lower = r$rho_lower %||% NA,
             rho_upper = r$rho_upper %||% NA,
             detected  = r$detected,               # point rule: 0 outside 95% CI
             rhat      = r$rhat %||% NA,
             n_diverg  = r$n_diverg %||% NA)))

## keep only converged fits before scoring
ok <- with(res, !is.na(detected) & (is.na(rhat) | rhat < 1.05))
r  <- res[ok, ]

## ---- point-null rule: 0 outside the 95% CI ----
tpr_point <- mean(r$detected)

## ---- sign-correct detections ----
r$sign_correct <- r$detected & (sign(r$rho_mean) == sign(r$true_rho))
tpr_signed <- mean(r$sign_correct)

cat(sprintf("\nValid fits ............. %d / %d\n", nrow(r), N_SIM))
cat(sprintf("TRUE POSITIVE RATE ..... %.3f   (0 outside 95%% CI)\n", tpr_point))
cat(sprintf("  sign-correct ......... %.3f\n", tpr_signed))

## power by sign -- should be ~symmetric
cat("\nTPR by true sign:\n")
print(round(tapply(r$detected, sign(r$true_rho), mean), 3))

## detection breakdown
cat("\nDetection breakdown:\n")
print(table(truth = ifelse(r$true_rho > 0, "rho>0", "rho<0"),
            call  = ifelse(r$detected, "detected", "missed")))

saveRDS(res, file.path(OUTPUT_DIR, "Revision-R1/FPR+FDR/tpr_results.rds"))

## ---------------- visualize ----------------
r$class <- ifelse(r$true_rho > 0,
                  sprintf("true rho = +%.2f", R_EFFECT),
                  sprintf("true rho = -%.2f", R_EFFECT))

p <- ggplot(r, aes(rho_mean, fill = detected)) +
  geom_histogram(bins = 50, colour = "white") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(aes(xintercept = true_rho), colour = "blue", linewidth = .3) +
  facet_wrap(~ class, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#1a9850"),
                    name = "0 excluded\nfrom 95% CI") +
  labs(title = sprintf("True DCGs (|rho| = %.2f):  TPR = %.3f", R_EFFECT, tpr_point),
       subtitle = "blue = true rho | dashed = 0",
       x = "posterior mean of rho", y = "count") +
  theme_bw()
print(p)
ggsave(file.path(OUTPUT_DIR, "tpr_histogram.png"), p, width = 7, height = 5, dpi = 150)


# check the results
res <- readRDS(file.path(OUTPUT_DIR, "Revision-R1/FPR+FDR/tpr_results.rds"))
ok  <- with(res, !is.na(detected) & (is.na(rhat) | rhat < 1.05))
r   <- res[ok, ]
cat("valid fits:", nrow(r), " detected:", sum(r$detected),
    " TPR:", round(mean(r$detected), 3), "\n")

# valid fits: 98  detected: 78  TPR: 0.796 
