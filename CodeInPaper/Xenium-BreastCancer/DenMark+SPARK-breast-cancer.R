## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # Xenium tables and QC gene lists, under their original sub-folders (see README.md)
OUTPUT_DIR <- "output"   # SPARK-X tables and the rank figures; share it with the other scripts of this folder
REPO_ROOT  <- "."        # root of the DenMark checkout
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

######################################################## 
# Real data analysis 
# In the large-range analysis 
########################################################

# we restrict only on several regions ----
## regions: 


# 1. First, we install SPARK package:--------- 
#install.packages('devtools')
#devtools::install_github('xzhoulab/SPARK')
#install.packages("promises")
#library(promises)
# finally: 
library(SPARK)
library('ggplot2')
library('MASS')
library(Matrix)
library(tidyverse)
library(magrittr)
library(dplyr)

library(doParallel)
registerDoParallel(cores = 5)


# the spatial information & the gene expression information: 
cell_gene <- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv'))
cell_spatial <- read.csv(file.path(DATA_DIR, '10X-breast-cancer/after-preprocess/cell_spatial.csv'))
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000

# check the genes in the dataset:
colnames(cell_gene)

# A. on dcis region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
gene_spatial_DCIS<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]

raw_count_mat <- as.matrix(gene_spatial_DCIS)[,-1]

gene_name_dcis<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/DCIS-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_dcis)]


# run spark-x: 
# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                                 center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 

sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-DCIS.csv') )



# B. on invasive region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_spatial_invasive<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]

raw_count_mat <- as.matrix(gene_spatial_invasive)[,-1]

gene_name_invasive<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Tumor-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_invasive)]


# run spark-x: 
# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                                 center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 

sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Tumor.csv') )


# C. on immune cells region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),]
gene_spatial_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]

raw_count_mat <- as.matrix(gene_spatial_immune)[,-1]

gene_name_immune<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Immune-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_immune)]


# run spark-x: 
# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                                 center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 
sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Immune.csv') )




# The intersection between spark and denmark identified genes -----
df_dcis_spark    <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-DCIS.csv'))
df_invasive_spark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Tumor.csv'))
df_immune_spark  <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Immune.csv'))


# order the p-value in an increasing rank 
df_dcis_order <- df_dcis_spark %>% 
  arrange(adj_pval)

df_invasive_order <- df_invasive_spark %>% 
  arrange(adj_pval)

df_immune_order <- df_immune_spark %>% 
  arrange(adj_pval)


# compare that with the denmark results: 
df_dcis_denmark   <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-DCIS.csv'))
df_invasive_denmark   <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Tumor.csv'))
df_immune_denmark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Immune.csv'))

# the ranking of the density correlated genes 
df_dcg_dcis<- data.frame(df_dcis_denmark$Gene,
                          df_dcis_denmark$Mean)

df_dcg_invasive<- data.frame(df_invasive_denmark$Gene,
                          df_invasive_denmark$Mean)

df_dcg_immune<- data.frame(df_immune_denmark$Gene,
                             df_immune_denmark$Mean)


# ordering: 
df_dcis_order_dcg <- df_dcg_dcis %>% 
  arrange(desc(abs(df_dcis_denmark.Mean)))

df_invasive_order_dcg <- df_dcg_invasive %>% 
  arrange(desc(abs(df_invasive_denmark.Mean)))

df_immune_order_dcg <- df_dcg_immune %>% 
  arrange(desc(abs(df_immune_denmark.Mean)))


df_dcis_order$SVGs005[1:10]
df_invasive_order$SVGs005[1:10]
df_immune_order$SVGs005[1:10]


df_dcis_order_dcg$df_dcis_denmark.Gene[1:10]
df_invasive_order_dcg$df_invasive_denmark.Gene[1:10]
df_immune_order_dcg$df_immune_denmark.Gene[1:10]


# the intersection between two datasets at each ROI: 
intersect(df_dcis_order$SVGs005[1:50], 
          df_dcis_order_dcg$df_dcis_denmark.Gene[1:50])

intersect(df_invasive_order$SVGs005[1:50], 
          df_invasive_order_dcg$df_invasive_denmark.Gene[1:50])

intersect(df_immune_order$SVGs005[1:50], 
          df_immune_order_dcg$df_immune_denmark.Gene[1:50])

# the percentage of the genes:
length(intersect(df_dcis_order$SVGs005[1:50], 
                 df_dcis_order_dcg$df_dcis_denmark.Gene[1:50]))/50

length(intersect(df_invasive_order$SVGs005[1:50], 
                 df_invasive_order_dcg$df_invasive_denmark.Gene[1:50]))/50

length(intersect(df_immune_order$SVGs005[1:50], 
                 df_immune_order_dcg$df_immune_denmark.Gene[1:50]))/50



##############################################################
#  the revised figure for comparison with the Figure 2I-----
################################################################

## ===========================================================================
## Rank–rank comparison: DenMark (DCG) vs SPARK-X (SVG), three cell types
## ===========================================================================
library(dplyr); library(ggplot2)
# library(ggrepel)   # only if you want gene labels (see note)

BASE <- file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes')

spark <- list(
  DCIS   = read.csv(file.path(BASE, 'SPARKSVG-DCIS.csv')),
  Tumor  = read.csv(file.path(BASE, 'SPARKSVG-Tumor.csv')),
  Immune = read.csv(file.path(BASE, 'SPARKSVG-Immune.csv')))
denmark <- list(
  DCIS   = read.csv(file.path(BASE, 'DCG-DCIS.csv')),
  Tumor  = read.csv(file.path(BASE, 'DCG-Tumor.csv')),
  Immune = read.csv(file.path(BASE, 'DCG-Immune.csv')))

LABELS <- c(DCIS = "DCIS", Tumor = "Invasive tumour", Immune = "CD8+ T cells")

## ---- join on genes ranked by BOTH methods, rank within the common set ----
build_ranks <- function(dm, sp, ct) {
  dm2 <- distinct(data.frame(gene = dm$Gene,     mean = dm$Mean))
  sp2 <- distinct(data.frame(gene = sp$SVGs005,  pval = sp$adj_pval))
  j <- merge(dm2, sp2, by = "gene")
  j <- j[is.finite(j$mean) & is.finite(j$pval), ]
  j$dcg_rank <- rank(-abs(j$mean), ties.method = "average")  # 1 = strongest density corr
  j$svg_rank <- rank(j$pval,       ties.method = "average")  # 1 = most spatially variable
  j$celltype <- ct
  j
}
dat <- bind_rows(Map(build_ranks, denmark, spark, names(denmark)))
dat$celltype <- factor(dat$celltype, levels = names(LABELS), labels = LABELS)

## ---- correlation per cell type ----
stats <- dat %>% group_by(celltype) %>% summarise(
  n   = n(),
  rho = cor(dcg_rank, svg_rank, method = "spearman"),
  p   = cor.test(dcg_rank, svg_rank, method = "spearman", exact = FALSE)$p.value,
  tau = cor(dcg_rank, svg_rank, method = "kendall"),
  .groups = "drop") %>%
  mutate(p_lab = ifelse(p < 1e-3, "P < 0.001", sprintf("P = %.3f", p)))
print(stats)     # exact numbers for the text

leg_lab <- setNames(
  sprintf("%s  (\u03c1 = %.2f, %s, n = %d)", stats$celltype, stats$rho, stats$p_lab, stats$n),
  as.character(stats$celltype))

cols <- c("DCIS" = "#E64B35", "Invasive tumour" = "#00A087", "CD8+ T cells" = "#3C5488")

## ---- one panel, coloured by cell type ----
p_rank <- ggplot(dat, aes(dcg_rank, svg_rank, colour = celltype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.55) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +  # optional trend
  scale_colour_manual(values = cols, labels = leg_lab, name = NULL) +
  scale_x_reverse() + scale_y_reverse() +           # rank 1 (best) at top-right
  labs(x = "DCG rank (DenMark, 1 = strongest density-correlation)",
       y = "SVG rank (SPARK-X, 1 = most spatially variable)") +
  theme_bw(base_size = 20) +
  theme(
    legend.position      = c(0.02, 0.98),          # top-left, inside panel
    legend.justification = c(0, 1),                # anchor legend's top-left corner
    legend.direction     = "vertical",             # stack the three cell types
    legend.background     = element_rect(fill = alpha("white", 0.75), colour = "grey80"),
    legend.key            = element_blank(),
    legend.text           = element_text(size = 18),
    panel.grid.minor      = element_blank()
  )

ggsave(file.path(BASE, "rank_dcg_vs_svg_by_celltype.png"),
       p_rank, width = 9, height = 8, dpi = 300)
p_rank



#########################################################
# Redraw into three panels 
#########################################################

## ===========================================================================
## Rank–rank comparison: DenMark (DCG) vs SPARK-X (SVG), one panel per population
##   - all genes DenMark analysed; non-SVGs tied last, drawn in a band below the SVG axis
##   - top-20 DCGs by |rho| highlighted, top few labelled, no regression lines
## ===========================================================================
library(dplyr); library(ggplot2); library(ggrepel); library(patchwork)

BASE <- file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes')

spark <- list(
  DCIS   = read.csv(file.path(BASE, 'SPARKSVG-DCIS.csv')),
  Tumor  = read.csv(file.path(BASE, 'SPARKSVG-Tumor.csv')),
  Immune = read.csv(file.path(BASE, 'SPARKSVG-Immune.csv')))
denmark <- list(
  DCIS   = read.csv(file.path(BASE, 'DCG-DCIS.csv')),
  Tumor  = read.csv(file.path(BASE, 'DCG-Tumor.csv')),
  Immune = read.csv(file.path(BASE, 'DCG-Immune.csv')))

LABELS  <- c(DCIS = "DCIS", Tumor = "Invasive tumour", Immune = "CD8+ T cells")
cols    <- c("DCIS" = "#E64B35", "Invasive tumour" = "#00A087", "CD8+ T cells" = "#3C5488")
TOP_N   <- 20     # highlighted DCGs per panel
N_LABEL <- 5      # of those, how many get a gene label
SVG_P   <- 0.05   # SPARK-X adjusted-P threshold for calling an SVG

## ---- ranks over ALL DenMark genes; non-SVGs tied at the last rank ----
build_ranks <- function(dm, sp, ct) {
  dm2 <- dm %>% transmute(gene = Gene, rho = Mean) %>%
    filter(!is.na(gene), is.finite(rho)) %>% distinct(gene, .keep_all = TRUE)
  sp2 <- sp %>% transmute(gene = SVGs005, pval = adj_pval) %>%
    filter(!is.na(gene), gene != "", is.finite(pval)) %>%
    group_by(gene) %>% summarise(pval = min(pval), .groups = "drop")
  left_join(dm2, sp2, by = "gene") %>%
    mutate(is_svg   = !is.na(pval) & pval < SVG_P,
           dcg_rank = rank(-abs(rho), ties.method = "average"),             # 1 = strongest |rho|
           svg_rank = rank(ifelse(is_svg, pval, Inf), ties.method = "average"),  # SVGs 1..n_svg; rest tied last
           top      = rank(-abs(rho), ties.method = "first") <= TOP_N,
           celltype = ct)
}
dat <- bind_rows(Map(build_ranks, denmark, spark, names(denmark)))
dat$celltype <- factor(dat$celltype, levels = names(LABELS), labels = LABELS)

## ---- numbers for the text ----
stats <- dat %>% group_by(celltype) %>% summarise(
  n_genes   = n(),
  n_svg     = sum(is_svg),
  n_not_svg = sum(!is_svg),
  top_svg   = sum(top & is_svg),                                   # of the top-20 DCGs, how many are SVGs
  r_s       = cor(dcg_rank, svg_rank, method = "spearman"),
  p         = cor.test(dcg_rank, svg_rank, method = "spearman", exact = FALSE)$p.value,
  tau_b     = cor(dcg_rank, svg_rank, method = "kendall"),
  r_s_svg_only = if (sum(is_svg) >= 3)
    cor(dcg_rank[is_svg], svg_rank[is_svg], method = "spearman") else NA_real_,
  .groups = "drop") %>%
  mutate(p_lab = ifelse(p < 1e-3, "P < 0.001", sprintf("P = %.3f", p)))
print(stats, width = Inf)

## ---- one panel per population ----
make_panel <- function(ct, show_ylab = TRUE) {
  d <- filter(dat, celltype == ct); s <- filter(stats, celltype == ct)
  col   <- cols[[as.character(ct)]]
  n_svg <- s$n_svg
  gap    <- max(3, 0.08 * n_svg)          # space between the last SVG and the band
  band_h <- max(2, 0.05 * n_svg)          # half-height of the band
  band_y <- n_svg + gap + band_h          # band centre (reversed axis -> below the SVGs)
  
  set.seed(1)
  d <- d %>% mutate(y = ifelse(is_svg, svg_rank, band_y + runif(n(), -band_h, band_h)))
  lab <- d %>% filter(top) %>% slice_min(dcg_rank, n = N_LABEL, with_ties = FALSE)
  
  yb <- pretty(c(1, max(n_svg, 1))); yb <- unique(c(1, yb[yb >= 1 & yb <= n_svg]))
  if (n_svg == 0) yb <- numeric(0)
  
  ggplot(d, aes(dcg_rank, y)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = n_svg + gap, ymax = band_y + band_h,
             fill = "grey94") +
    geom_point(data = filter(d, !top), colour = "grey70", size = 1.2, alpha = 0.6) +
    geom_point(data = filter(d,  top), colour = col, size = 2.4) +
    geom_text_repel(data = lab, aes(label = gene), colour = col, size = 3.6,
                    fontface = "italic", box.padding = 0.4, min.segment.length = 0,
                    max.overlaps = Inf, seed = 1) +
    scale_x_reverse(expand = expansion(mult = 0.03)) +
    scale_y_reverse(breaks = c(yb, band_y),
                    labels = c(yb, sprintf("Not SVG\n(n = %d)", s$n_not_svg)),
                    expand = expansion(mult = 0.03)) +
    labs(title = as.character(ct),
         subtitle = sprintf("%d genes, %d SVGs; top-%d DCGs: %d SVGs\nSpearman r = %.2f, %s",
                            s$n_genes, n_svg, TOP_N, s$top_svg, s$r_s, s$p_lab),
         x = "DCG rank (DenMark, 1 = strongest density correlated)",
         y = if (show_ylab) "SVG rank (SPARK-X, 1 = most spatially variable)" else NULL) +
    theme_classic(base_size = 14) +
    theme(plot.title    = element_text(face = "bold", colour = col),
          plot.subtitle = element_text(size = 11, colour = "grey25"))
}

cts <- levels(dat$celltype)
p_rank <- wrap_plots(Map(make_panel, cts, seq_along(cts) == 1), nrow = 1) +
  plot_annotation(caption = sprintf(
    "Coloured: top %d DCGs by |\u03c1| per population; labelled: top %d. Grey band: genes not called SVGs by SPARK-X (adj. P \u2265 %.2f), tied at the last SVG rank.",
    TOP_N, N_LABEL, SVG_P))

ggsave(file.path(BASE, "rank_dcg_vs_svg_by_celltype_3panel.png"),
       p_rank, width = 16, height = 5.8, dpi = 300)
p_rank



## ---- counts for the text: where do the top-20 DCGs fall? ----
cnt <- dat %>% group_by(celltype) %>% summarise(
  n_genes       = n(),
  n_svg         = sum(is_svg),
  mid           = n_svg / 2,
  top_upper     = sum(top & is_svg & svg_rank <= mid),
  top_lower_svg = sum(top & is_svg & svg_rank >  mid),   # lower half, SVGs only
  top_not_svg   = sum(top & !is_svg),
  top_lower_all = top_lower_svg + top_not_svg,           # lower half incl. non-SVGs (tied last)
  .groups = "drop")
print(cnt, width = Inf)

# genes behind each number, to check against Data S2
dat %>% filter(top) %>% arrange(celltype, dcg_rank) %>%
  left_join(dplyr::select(cnt, celltype, mid), by = "celltype") %>%
  mutate(position = case_when(!is_svg ~ "not SVG",
                              svg_rank > mid ~ "lower half",
                              TRUE ~ "upper half")) %>%
  dplyr::select(celltype, gene, rho, dcg_rank, svg_rank, position) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

## ---- panel subtitles carry the population-specific count ----
text_lab <- function(ct) {
  k <- filter(cnt, celltype == ct)
  switch(as.character(ct),
         "DCIS"            = sprintf("%d of %d top DCGs in lower half of SVG ranking", k$top_lower_all, TOP_N),
         "Invasive tumour" = sprintf("%d of %d top DCGs in lower half of SVG ranking", k$top_lower_all, TOP_N),
         "CD8+ T cells"    = sprintf("%d of %d top DCGs not SVGs", k$top_not_svg, TOP_N))
}

make_panel2 <- function(ct, show_ylab = TRUE) {
  s <- filter(stats, celltype == ct); k <- filter(cnt, celltype == ct)
  make_panel(ct, show_ylab) +
    geom_hline(yintercept = k$mid, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    labs(subtitle = sprintf("%d genes, %d SVGs\n%s",
                            s$n_genes, s$n_svg, text_lab(ct)))
}

p_rank <- wrap_plots(Map(make_panel2, cts, seq_along(cts) == 1), nrow = 1) +
  plot_annotation(caption = sprintf(
    "Coloured: top %d DCGs by |\u03c1| per population; labelled: top %d. Dashed line: midpoint of the SVG ranking. Grey band: genes not called SVGs by SPARK-X (adj. P \u2265 %.2f), tied at the last rank and counted in the lower half.",
    TOP_N, N_LABEL, SVG_P))

ggsave(file.path(BASE, "rank_dcg_vs_svg_by_celltype_3panel.png"),
       p_rank, width = 16, height = 6.2, dpi = 300)
p_rank





#########################################################
# Real data analysis: Identifying SVGs with SPARK
# In the short-analysis part: 
#########################################################


# 1. First, we install SPARK package:--------- 
#install.packages('devtools')
#devtools::install_github('xzhoulab/SPARK')
#install.packages("promises")
#library(promises)
# finally: 
library(SPARK)
library('ggplot2')
library('MASS')
library(Matrix)
library(tidyverse)
library(magrittr)
library(dplyr)


library(doParallel)
registerDoParallel(cores = 5)



# the spatial information & the gene expression information: 
cell_gene <- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv'))
cell_spatial <- read.csv(file.path(DATA_DIR, '10X-breast-cancer/after-preprocess/cell_spatial.csv'))
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000

# check the genes in the dataset:
colnames(cell_gene)



# B. We then generate the gene sets for three regions for the SPARK-X -----
## 1. DCIS-I
## 2. DCIS-II
## 3. invasive tumor 


## B1. the region of interest: DCIS-I -----
center_slice_cell_location<- cell_spatial[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5 & cell_spatial$y>= 1.8 & cell_spatial$y <= 2.8),]
center_slice_gene_location<- cell_gene[which(cell_spatial$x>=6.5 & cell_spatial$x<= 7.5 & cell_spatial$y>= 1.8 & cell_spatial$y <= 2.8),]

raw_count_mat <- as.matrix(center_slice_gene_location)[,-1]
# QC: keep the genes with criteria: 
## (1) variance across cells ranked in the top 95% 
## (2) non-zero values on at least 5% of locations 


### (1)
#keep_genes_1 <- colSums(raw_count_mat > 0) > 0.05 * nrow(raw_count_mat)
#raw_count_mat_1 <- raw_count_mat[,keep_genes_1]

### (2)
#raw_count_mat_2 <-  raw_count_mat_1[,apply(raw_count_mat_1, 2,var) >= quantile(apply(raw_count_mat_1, 2,var), 0.05)]

gene_name_dcis1<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-I-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_dcis1)]


# locations: 23624 cell locations, with x and y coordinates 
str(center_slice_cell_location)
# gene expression, 649 genes with 23624 cell locations 
str(raw_count_mat_2)

# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                      center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 

sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS1.csv') )







## B2. the region of interest: DCIS-II -----
center_slice_cell_location<- cell_spatial[which(cell_spatial$x>=5 & cell_spatial$x<= 6 & cell_spatial$y>= 0.75 & cell_spatial$y <= 1.75),]
center_slice_gene_location<- cell_gene[which(cell_spatial$x>=5 & cell_spatial$x<= 6 & cell_spatial$y>= 0.75 & cell_spatial$y <= 1.75),]



raw_count_mat <- as.matrix(center_slice_gene_location)[,-1]
# QC: keep the genes with criteria: 
## (1) variance across cells ranked in the top 95% 
## (2) non-zero values on at least 5% of locations 


### (1)
#keep_genes_1 <- colSums(raw_count_mat > 0) > 0.05 * nrow(raw_count_mat)
#raw_count_mat_1 <- raw_count_mat[,keep_genes_1]
### (2)
#raw_count_mat_2 <-  raw_count_mat_1[,apply(raw_count_mat_1, 2,var) >= quantile(apply(raw_count_mat_1, 2,var), 0.05)]


gene_name_dcis2<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-II-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_dcis2)]



# locations: 23624 cell locations, with x and y coordinates 
str(center_slice_cell_location)
# gene expression, 649 genes with 23624 cell locations 
str(raw_count_mat_2)

# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                                 center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 

sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS2.csv') )




## B3. the region of interest: invasive tumor -----
center_slice_cell_location<- cell_spatial[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y <= 3),]
center_slice_gene_location<- cell_gene[which(cell_spatial$x>=0.5 & cell_spatial$x<= 1.5 & cell_spatial$y>= 2 & cell_spatial$y <= 3),]

raw_count_mat <- as.matrix(center_slice_gene_location)[,-1]
# QC: keep the genes with criteria: 
## (1) variance across cells ranked in the top 95% 
## (2) non-zero values on at least 5% of locations 


### (1)
#keep_genes_1 <- colSums(raw_count_mat > 0) > 0.05 * nrow(raw_count_mat)
#raw_count_mat_1 <- raw_count_mat[,keep_genes_1]

### (2)
#raw_count_mat_2 <-  raw_count_mat_1[,apply(raw_count_mat_1, 2,var) >= quantile(apply(raw_count_mat_1, 2,var), 0.05)]


gene_name_invasive<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/Invasive-after-QC.csv'))$x
raw_count_mat_2<- raw_count_mat[,which(colnames(raw_count_mat) %in% gene_name_invasive)]



# locations: 23624 cell locations, with x and y coordinates 
str(center_slice_cell_location)
# gene expression, 649 genes with 23624 cell locations 
str(raw_count_mat_2)

# 2. load the dataset, and identify the SVGs with SPARK in the dataset 
slice_spatial<- data.frame(cbind(center_slice_cell_location$x,
                                 center_slice_cell_location$y))
colnames(slice_spatial)<- c('x','y')

# use spark-X: 

sparkX <- sparkx(t(raw_count_mat_2),
                 slice_spatial,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


adj_pval<-sparkX$res_mtest$adjustedPval[which(sparkX$res_mtest$adjustedPval<0.05)]


SVGs005<-rownames(t(raw_count_mat_2))[which(sparkX$res_mtest$adjustedPval<0.05)]



# save SVGs and their corresponding p-values 
write.csv(data.frame(SVGs005,adj_pval),
          file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-invasive.csv') )




# the ranking of genes from SPARK and from DennMark ----
df_dcis1_spark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS1.csv'))
df_dcis2_spark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS2.csv'))
df_invasive_spark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-invasive.csv'))


# order the p-value in an increasing rank 
df_dcis1_order <- df_dcis1_spark %>% 
  arrange(adj_pval)

df_dcis2_order <- df_dcis2_spark %>% 
  arrange(adj_pval)

df_invasive_order <- df_invasive_spark %>% 
  arrange(adj_pval)


# compare that with the denmark results: 
df_dcis1_denmark   <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS1.csv'))
df_dcis2_denmark   <- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS2.csv'))
df_invasive_denmark<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-invasive.csv'))

# the ranking of the density correlated genes 
df_dcg_dcis1<- data.frame(df_dcis1_denmark$Gene,
                          df_dcis1_denmark$Mean)

df_dcg_dcis2<- data.frame(df_dcis2_denmark$Gene,
                          df_dcis2_denmark$Mean)

df_dcg_invasive<- data.frame(df_invasive_denmark$Gene,
                             df_invasive_denmark$Mean)


# ordering: 
df_dcis1_order_dcg <- df_dcg_dcis1 %>% 
  arrange(desc(abs(df_dcis1_denmark.Mean)))

df_dcis2_order_dcg <- df_dcg_dcis2 %>% 
  arrange(desc(abs(df_dcis2_denmark.Mean)))

df_invasive_order_dcg <- df_dcg_invasive %>% 
  arrange(desc(abs(df_invasive_denmark.Mean)))


df_dcis1_order$SVGs005[1:10]
df_dcis2_order$SVGs005[1:10]
df_invasive_order$SVGs005[1:10]


df_dcis1_order_dcg$df_dcis1_denmark.Gene[1:10]
df_dcis2_order_dcg$df_dcis2_denmark.Gene[1:10]
df_invasive_order_dcg$df_invasive_denmark.Gene[1:10]


# the intersection between two datasets at each ROI: 
intersect(df_dcis1_order$SVGs005[1:50], 
             df_dcis1_order_dcg$df_dcis1_denmark.Gene[1:50])

intersect(df_dcis2_order$SVGs005[1:50], 
          df_dcis2_order_dcg$df_dcis2_denmark.Gene[1:50])

intersect(df_invasive_order$SVGs005[1:50], 
          df_invasive_order_dcg$df_invasive_denmark.Gene[1:50])



