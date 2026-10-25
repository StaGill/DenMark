## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds Revision-R1/figures/Figure4Hsheet.xlsx (see "Missing files" in README.md)
OUTPUT_DIR <- "output"   # Fig4H_forest_LKJ.pdf and .png are written here
REPO_ROOT  <- "."        # root of the DenMark checkout
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2); library(ggsci)
})
set.seed(2024)

## ---- 1. Read the sheet and pull the requested genes -----------------------
XLSX  <- file.path(DATA_DIR, "Revision-R1/figures/Figure4Hsheet.xlsx")                     # edit path if needed
GENES <- c("LTB","CCR7","TCF7","CD27")

raw <- read_excel(XLSX, sheet = "Figure4H", col_names = FALSE)
names(raw)[1:4] <- c("gene","mean","lo","hi")
dat <- raw %>%
  mutate(across(c(mean, lo, hi), as.numeric)) %>%
  filter(!is.na(mean), gene %in% GENES) %>%        # drops the "DCG-Immune" title row
  mutate(gene = as.character(gene))
stopifnot(setequal(dat$gene, GENES))              # tripwire: all four genes found
print(dat)

## ---- 2. LKJ(1.2) prior on the correlation, from rethinking ----------------
ETA <- 1.2; NDRAW <- 1e5
if (requireNamespace("rethinking", quietly = TRUE)) {
  R   <- rethinking::rlkjcorr(NDRAW, K = 2, eta = ETA)   # NDRAW x 2 x 2
  rho <- R[, 1, 2]
} else {
  ## exact equivalent for K = 2: (rho + 1)/2 ~ Beta(eta, eta)
  rho <- 2 * rbeta(NDRAW, ETA, ETA) - 1
}
prior <- data.frame(gene = sprintf("LKJ(%.1f) prior", ETA),
                    mean = mean(rho),
                    lo   = unname(quantile(rho, .025)),
                    hi   = unname(quantile(rho, .975)))
print(prior)                                      # ~ 0, [-0.92, 0.92]

## ---- 3. Classify: positive / non-informative / negative -------------------
plot_df <- bind_rows(dat, prior) %>%
  mutate(group = case_when(lo > 0 ~ "Positive",
                           hi < 0 ~ "Negative",
                           TRUE   ~ "Non-informative"),
         group = factor(group, levels = c("Positive","Non-informative","Negative")),
         gene  = factor(gene, levels = rev(c(GENES, prior$gene))))  # prior at bottom

GRP_COL <- c(Positive = "#00A087", `Non-informative` = "grey55", Negative = "#E64B35")

## ---- 4. Forest plot (larger fonts, no "Negative" legend key) --------------
BASE <- 13                                            # base font size for everything

theme_nature <- function(b = BASE) theme_classic(base_size = b, base_family = "sans") +
  theme(axis.line    = element_line(linewidth = .4),
        axis.ticks   = element_line(linewidth = .4),
        axis.text    = element_text(colour = "black", size = b),
        axis.title   = element_text(size = b + 1),
        plot.title   = element_text(size = b + 2, face = "bold"),
        legend.text  = element_text(size = b),
        legend.title = element_text(size = b))

fig <- ggplot(plot_df, aes(x = mean, y = gene, colour = group)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = .5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .25, linewidth = .8) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%.2f", mean)), vjust = -1.1, size = BASE / 2.8,
            show.legend = FALSE) +
  scale_colour_manual(values = GRP_COL, drop = TRUE, name = NULL) +   # only groups present
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
  labs(x = expression("Posterior of " * rho * " (mean and 95% CI)"), y = NULL,
       title = expression("DenMark fitted posterior for " * rho * " vs LKJ(1.2) prior")) +
  theme_nature() +
  theme(legend.position = "top")

print(fig)
ggsave(file.path(OUTPUT_DIR, "Fig4H_forest_LKJ.pdf"), fig, width = 160, height = 120, units = "mm")
ggsave(file.path(OUTPUT_DIR, "Fig4H_forest_LKJ.png"), fig, width = 160, height = 120, units = "mm", dpi = 320)
