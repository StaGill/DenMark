## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # not used: the input is DCG-Immune.csv from Viz-for-cluster-genes-long-range.R
OUTPUT_DIR <- "output"   # holds DCG-Immune.csv; the GSEA dot plot is written next to it
REPO_ROOT  <- "."        # root of the DenMark checkout
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

#----------------------------------------
# GSEA — immune (CD8+ T cells), matched to the brain/tumour pipeline
#----------------------------------------
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(ggplot2); library(dplyr); library(stringr)

## ---------------------------------------------------------------------------
## 1. Ranked gene list from DenMark Mean (immune, whole tissue)
## ---------------------------------------------------------------------------
d <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Immune.csv'),
              stringsAsFactors = FALSE)
geneList <- setNames(d$Mean, d$Gene)
geneList <- sort(geneList[!duplicated(names(geneList))], decreasing = TRUE)

## human symbols -> ENSEMBL (clean, alignment-safe)
gi  <- bitr(names(geneList), "SYMBOL", "ENSEMBL", org.Hs.eg.db)
gi  <- gi[!duplicated(gi$SYMBOL), ]
s2e <- setNames(gi$ENSEMBL, gi$SYMBOL)
geneList <- geneList[names(geneList) %in% names(s2e)]
names(geneList) <- s2e[names(geneList)]
geneList <- sort(geneList[!duplicated(names(geneList))], decreasing = TRUE)
cat(sprintf("genes into GSEA: %d\n", length(geneList)))

## ---------------------------------------------------------------------------
## 2. GSEA on GO — HUMAN, same settings as the brain run
## ---------------------------------------------------------------------------
set.seed(2026)
gse <- gseGO(geneList      = geneList,
             OrgDb         = org.Hs.eg.db,
             keyType       = "ENSEMBL",
             ont           = "ALL",
             minGSSize     = 5,
             maxGSSize     = 300,
             pvalueCutoff  = 0.05,              # same as brain
             pAdjustMethod = "BH",
             eps           = 0,
             verbose       = FALSE)

cat(sprintf("enriched sets (p.adjust < 0.05): %d\n", nrow(as.data.frame(gse))))

## ---------------------------------------------------------------------------
## 3. THE FIGURE — same faceted dotplot as the brain (only if terms survive)
## ---------------------------------------------------------------------------
if (nrow(as.data.frame(gse)) > 0) {
  p_dot <- dotplot(gse,
                   showCategory = 5,
                   split        = ".sign",
                   x            = "GeneRatio",
                   color        = "p.adjust",
                   size         = "Count") +
    facet_grid(. ~ .sign) +
    scale_colour_gradient(low = "#2166AC", high = "#B2182B", name = "p.adjust") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 28)) +
    labs(x = "GeneRatio", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey60"),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 20, colour = "black"),
      axis.text.x      = element_text(size = 20, colour = "black"),
      axis.title.x     = element_text(size = 20),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size = 18),
      legend.text      = element_text(size = 15)
    )
  
  ggsave(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/12_gsea_dotplot_Immune.png'),
         p_dot, width = 8.5, height = 5, dpi = 300)
  print(p_dot)
} else {
  cat("\nNo GO term passes p.adjust < 0.05 for the immune panel.\n",
      "Check length(geneList): if it is small (~150-170), GSEA is underpowered here too.\n")
}



gse_sym <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
res <- gse_sym@result
for (t in c("lymphocyte activation","immune response","T cell differentiation")) {
  i <- which(tolower(res$Description) == t)[1]
  if (!is.na(i)) cat("\n===", res$Description[i],
                     sprintf("(NES=%.2f, p.adj=%.3g) ===\n", res$NES[i], res$p.adjust[i]),
                     gsub("/", ", ", res$core_enrichment[i]), "\n")
}

