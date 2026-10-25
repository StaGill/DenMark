## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # holds full-division.xlsx (see README.md)
OUTPUT_DIR <- "output"   # reads the summaries in BC-whole/ (summarize_BCwhole_Astar21.R) and writes there
REPO_ROOT  <- "."        # root of the DenMark checkout
## -----------------------------------------------------------------------------

###############################################################################
# Merged "tumor" division vs fine DCIS / Invasive-Tumor divisions (whole slice)
# ---------------------------------------------------------------------------
# QUESTION: does pooling several epithelial cell types into one "tumor" density
# field reproduce the DCGs you get from the fine divisions (DCIS, Invasive
# Tumor) analysed separately -- or does merging distort the coupling?
#
# DCG = density-coupled gene; Astar[2,1] sign = coupling direction;
# credible = 95% CI excludes 0.
#
# INPUTS (all present)
#   Tumor_combined      : BC-whole-tumor_Astar21_summary.csv   (the merged run)
#   DCIS_full           : full-division.xlsx, DCIS block         (gene col 2)
#   InvasiveTumor_full  : full-division.xlsx, Invasive-Tumor blk (gene col 12)
#
# OUTPUTS: composition bars, concordance scatter (merged vs each fine type),
#          concordance table, long table.
###############################################################################

DATA_DIR <- normalizePath(DATA_DIR, mustWork = FALSE)
setwd(file.path(OUTPUT_DIR, "BC-whole"))

## ---- 0. packages -----------------------------------------------------------
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl); library(ggplot2)
out_dir <- "."; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. analyses (edit paths if needed) ------------------------------------
analyses <- list(
  Tumor_combined     = list(kind = "csv",        path = "BC-whole-tumor_Astar21_summary.csv"),
  DCIS_full          = list(kind = "xlsx_block", path = file.path(DATA_DIR, "full-division.xlsx"), block = 2),
  InvasiveTumor_full = list(kind = "xlsx_block", path = file.path(DATA_DIR, "full-division.xlsx"), block = 12)
)
ref_analysis   <- "Tumor_combined"
fine_analyses  <- c("DCIS_full", "InvasiveTumor_full")
plot_order     <- c("Tumor_combined", "DCIS_full", "InvasiveTumor_full")

## ---- 2. loaders ------------------------------------------------------------
is_credible <- function(lo, hi) (lo > 0 & hi > 0) | (lo < 0 & hi < 0)

load_astar <- function(spec) {
  if (spec$kind == "csv") {
    d <- read.csv(spec$path, check.names = FALSE)
    out <- data.frame(gene = trimws(as.character(d$gene_name)),
                      mean = d$A21_mean, lower = d$A21_q2.5, upper = d$A21_q97.5,
                      stringsAsFactors = FALSE)
  } else {  # xlsx_block: row1 titles, row2 headers, row3+ data
    raw <- as.data.frame(read_excel(spec$path, sheet = 1, col_names = FALSE, col_types = "text"))
    raw <- raw[-(1:2), , drop = FALSE]; c0 <- spec$block
    out <- data.frame(gene  = trimws(as.character(raw[[c0]])),
                      mean  = suppressWarnings(as.numeric(raw[[c0 + 1]])),
                      lower = suppressWarnings(as.numeric(raw[[c0 + 2]])),
                      upper = suppressWarnings(as.numeric(raw[[c0 + 3]])),
                      stringsAsFactors = FALSE)
  }
  out <- out[!is.na(out$gene) & out$gene != "" & !is.na(out$mean), ]
  out$credible <- is_credible(out$lower, out$upper)
  out
}
loaded <- lapply(analyses, load_astar)

## ---- 3. composition (positive / negative / non-informative) ----------------
cols3 <- c("positive DCG"="#c0392b","negative DCG"="#2c7fb8","non-informative"="#bdc3c7")
comp <- do.call(rbind, lapply(plot_order, function(nm){
  d <- loaded[[nm]]
  data.frame(analysis = nm,
             category = c("positive DCG","negative DCG","non-informative"),
             n = c(sum(d$credible & d$mean>0), sum(d$credible & d$mean<0), sum(!d$credible)))
}))
comp$category <- factor(comp$category, levels = names(cols3))
comp$analysis <- factor(comp$analysis, levels = plot_order)
tot <- do.call(rbind, lapply(plot_order, function(nm)
  data.frame(analysis = nm, n_tested = nrow(loaded[[nm]]), n_dcg = sum(loaded[[nm]]$credible))))
tot$analysis <- factor(tot$analysis, levels = plot_order)

p_comp <- ggplot(comp, aes(analysis, n, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(n>0,n,"")), position = position_stack(vjust=0.5),
            size = 3, color = "white") +
  geom_text(data = tot, inherit.aes = FALSE,
            aes(analysis, n_tested, label = paste0(n_dcg," DCGs")),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cols3) +
  scale_x_discrete(labels = c(
    Tumor_combined     = "Tumor (combined types)",
    DCIS_full          = "DCIS",
    InvasiveTumor_full = "Invasive tumor")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Merged tumor vs fine DCIS / Invasive-Tumor divisions",
       x = NULL, y = "number of genes", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out_dir, "fig_composition_merged_vs_fine.png"), p_comp, width = 8, height = 6, dpi = 200)


## ---- 4. concordance: merged vs each fine type ------------------------------
R <- loaded[[ref_analysis]]
pts <- do.call(rbind, lapply(fine_analyses, function(o){
  m <- merge(R[,c("gene","mean","credible")], loaded[[o]][,c("gene","mean","credible")],
             by = "gene", suffixes = c("_ref","_oth"))
  m$fine <- o
  m$concord <- with(m, ifelse(credible_ref & credible_oth & sign(mean_ref)==sign(mean_oth), "both, same sign",
                       ifelse(credible_ref & credible_oth, "both, opposite",
                       ifelse(credible_ref & !credible_oth, paste0(ref_analysis," only"),
                       ifelse(!credible_ref & credible_oth, "fine only", "neither")))))
  m
}))
stats <- do.call(rbind, lapply(split(pts, pts$fine), function(m){
  both <- m$credible_ref & m$credible_oth
  data.frame(reference = ref_analysis, fine = m$fine[1], n_shared = nrow(m),
             spearman_r = round(cor(m$mean_ref, m$mean_oth, method="spearman"),3),
             n_both_credible = sum(both),
             pct_sign_agree = ifelse(sum(both)==0, NA,
               round(100*mean(sign(m$mean_ref[both])==sign(m$mean_oth[both])),1)),
             n_opposite = sum(both & sign(m$mean_ref)!=sign(m$mean_oth)),
             stringsAsFactors = FALSE)
}))
print(stats); write.csv(stats, file.path(out_dir, "DCG_merged_vs_fine_concordance.csv"), row.names = FALSE)

lab <- setNames(sprintf("%s\n(r=%.2f, sign agree=%s%%)", stats$fine, stats$spearman_r, stats$pct_sign_agree), stats$fine)
pts$fine_lab <- lab[pts$fine]
p_conc <- ggplot(pts, aes(mean_ref, mean_oth, color = concord)) +
  geom_hline(yintercept = 0, color="grey85") + geom_vline(xintercept = 0, color="grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color="grey60") +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ fine_lab) +
  labs(title = "Does merging preserve each gene's coupling?",
       x = "Astar[2,1] — Tumor_combined", y = "Astar[2,1] — fine division", color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "fig_concordance_merged_vs_fine.png"), p_conc, width = 10, height = 5.2, dpi = 200)

## ---- 5. gene-level detail (opposite / merged-only genes worth inspecting) ---
wide <- Reduce(function(a,b) merge(a,b,by="gene",all=TRUE), lapply(plot_order, function(nm){
  d <- loaded[[nm]][,c("gene","mean","credible")]; names(d)[2:3] <- paste0(c("mean_","cred_"), nm); d
}))
write.csv(wide, file.path(out_dir, "DCG_merged_vs_fine_long.csv"), row.names = FALSE)

## ---- 6. verdict ------------------------------------------------------------
verdict <- c(
"=================================================================",
" MERGED TUMOR vs FINE DIVISIONS -- SUMMARY",
"=================================================================",
sprintf("Tumor_combined : %d DCGs (+%d/-%d), %d non-informative of %d",
  sum(R$credible), sum(R$credible & R$mean>0), sum(R$credible & R$mean<0), sum(!R$credible), nrow(R)),
apply(stats, 1, function(z) sprintf("vs %-19s r=%s  sign-agree=%s%%  opposite=%s (of %s both-credible)",
  z[["fine"]], z[["spearman_r"]], z[["pct_sign_agree"]], z[["n_opposite"]], z[["n_both_credible"]])),
"",
"READING",
"  - Merged tumor tracks whichever fine type dominates the pooled density",
"    field: higher correlation with that subtype, lower with the other.",
"  - Merging adds negative DCGs the fine types don't show (genes high in one",
"    compartment become 'negatively coupled' to the pooled density peak).",
"  - => pooling cell types is NOT signal-preserving; prefer the fine division.",
"=================================================================")
writeLines(verdict, file.path(out_dir, "DCG_merged_vs_fine_verdict.txt"))
cat("\n", paste(verdict, collapse = "\n"), "\n")

print(p_comp); print(p_conc)
cat("\nDone. Figures, concordance table, long table, verdict written to:", out_dir, "\n")
