## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # Xenium tables and per-gene fits, under their original sub-folders (see README.md)
OUTPUT_DIR <- "output"   # whole-section DCG tables written here; share it with the other scripts of this folder
REPO_ROOT  <- "."        # root of the DenMark checkout
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

#############################
# viz for BC long range analysis
#############################

# load packages: 
library('ggplot2')
library('rstan')
library('MASS')
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
library(loo)
# for viz:
library(cowplot)
library(viridis)
library(rethinking)
library(reshape2)
library(readxl)

# quantiles:
q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)


#gene_tab<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv'))

#gene_name_ls<-colnames(gene_tab)[-1]


# all gene names: 


# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/large-range/results-dcis")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)



# gene name after filtering: 
gene_name_ls<- read.csv(file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/DCIS-after-QC.csv'))$x

# Extract gene names from filenames
extract_gene <- function(filename) {
  base <- basename(filename)                 # e.g. "ZNF562_M2_313.csv"
  sub("_M2_.*$", "", base)                   # → "ZNF562"
}

file_genes <- sapply(files, extract_gene)

# Keep only files whose gene is in gene_name_ls
keep_idx <- file_genes %in% gene_name_ls
files_filtered <- files[keep_idx]
file_genes_filtered <- file_genes[keep_idx]

# -----------------------------------------------------
# Reorder filtered files to match order of gene_name_ls
# -----------------------------------------------------
order_idx <- match(gene_name_ls, file_genes_filtered)

# Remove genes that have no file
valid <- !is.na(order_idx)
order_idx <- order_idx[valid]
ordered_genes <- gene_name_ls[valid]

files_ordered <- files_filtered[order_idx]

# -----------------------------------------------------
# Initialize result vectors
# -----------------------------------------------------
mean_vec  <- numeric(length(files_ordered))
low95_vec <- numeric(length(files_ordered))
high95_vec <- numeric(length(files_ordered))

# -----------------------------------------------------
# Loop through ordered files
# -----------------------------------------------------
for (i in seq_along(files_ordered)) {
  
  df <- read.csv(files_ordered[i], stringsAsFactors = FALSE)
  
  row <- df[df$variable == "Astar[2,1]", ]
  
  if (nrow(row) == 0) {
    warning("Row Astar[2,1] not found for gene: ", ordered_genes[i])
    next
  }
  
  mean_vec[i]  <- as.numeric(row$mean)
  low95_vec[i] <- as.numeric(row$`X2.5.`)
  high95_vec[i] <- as.numeric(row$`X97.5.`)
}

# -----------------------------------------------------
# Output final data frame
# -----------------------------------------------------
result_df_dcis <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)

colnames(result_df_dcis)<- c('Gene','Mean','Low95','High95')

write.csv(result_df_dcis, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-DCIS.csv'))


#######################################
# for invasive tumor
#######################################


# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/large-range/results-invasive-tumor")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)

# gene name after filtering: 
gene_name_ls<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Tumor-after-QC.csv'))$x

# Extract gene names from filenames
extract_gene <- function(filename) {
  base <- basename(filename)                 # e.g. "ZNF562_M2_313.csv"
  sub("_M2_.*$", "", base)                   # → "ZNF562"
}

file_genes <- sapply(files, extract_gene)

# Keep only files whose gene is in gene_name_ls
keep_idx <- file_genes %in% gene_name_ls
files_filtered <- files[keep_idx]
file_genes_filtered <- file_genes[keep_idx]

# -----------------------------------------------------
# Reorder filtered files to match order of gene_name_ls
# -----------------------------------------------------
order_idx <- match(gene_name_ls, file_genes_filtered)

# Remove genes that have no file
valid <- !is.na(order_idx)
order_idx <- order_idx[valid]
ordered_genes <- gene_name_ls[valid]

files_ordered <- files_filtered[order_idx]

# -----------------------------------------------------
# Initialize result vectors
# -----------------------------------------------------
mean_vec  <- numeric(length(files_ordered))
low95_vec <- numeric(length(files_ordered))
high95_vec <- numeric(length(files_ordered))

# -----------------------------------------------------
# Loop through ordered files
# -----------------------------------------------------
for (i in seq_along(files_ordered)) {
  
  df <- read.csv(files_ordered[i], stringsAsFactors = FALSE)
  
  row <- df[df$variable == "Astar[2,1]", ]
  
  if (nrow(row) == 0) {
    warning("Row Astar[2,1] not found for gene: ", ordered_genes[i])
    next
  }
  
  mean_vec[i]  <- as.numeric(row$mean)
  low95_vec[i] <- as.numeric(row$`X2.5.`)
  high95_vec[i] <- as.numeric(row$`X97.5.`)
}

# -----------------------------------------------------
# Output final data frame
# -----------------------------------------------------
result_df_tumor <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)

colnames(result_df_tumor)<- c('Gene','Mean','Low95','High95')


write.csv(result_df_tumor, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Tumor.csv'))



#######################################
# for Immune cells
#######################################


# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/large-range/results-immune")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)

# gene name after filtering: 
gene_name_ls<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Immune-after-QC.csv'))$x

# Extract gene names from filenames
extract_gene <- function(filename) {
  base <- basename(filename)                 # e.g. "ZNF562_M2_313.csv"
  sub("_M2_.*$", "", base)                   # → "ZNF562"
}

file_genes <- sapply(files, extract_gene)

# Keep only files whose gene is in gene_name_ls
keep_idx <- file_genes %in% gene_name_ls
files_filtered <- files[keep_idx]
file_genes_filtered <- file_genes[keep_idx]

# -----------------------------------------------------
# Reorder filtered files to match order of gene_name_ls
# -----------------------------------------------------
order_idx <- match(gene_name_ls, file_genes_filtered)

# Remove genes that have no file
valid <- !is.na(order_idx)
order_idx <- order_idx[valid]
ordered_genes <- gene_name_ls[valid]

files_ordered <- files_filtered[order_idx]

# -----------------------------------------------------
# Initialize result vectors
# -----------------------------------------------------
mean_vec  <- numeric(length(files_ordered))
low95_vec <- numeric(length(files_ordered))
high95_vec <- numeric(length(files_ordered))

# -----------------------------------------------------
# Loop through ordered files
# -----------------------------------------------------
for (i in seq_along(files_ordered)) {
  
  df <- read.csv(files_ordered[i], stringsAsFactors = FALSE)
  
  row <- df[df$variable == "Astar[2,1]", ]
  
  if (nrow(row) == 0) {
    warning("Row Astar[2,1] not found for gene: ", ordered_genes[i])
    next
  }
  
  mean_vec[i]  <- as.numeric(row$mean)
  low95_vec[i] <- as.numeric(row$`X2.5.`)
  high95_vec[i] <- as.numeric(row$`X97.5.`)
}

# -----------------------------------------------------
# Output final data frame
# -----------------------------------------------------
result_df_immune <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)

colnames(result_df_immune)<- c('Gene','Mean','Low95','High95')


write.csv(result_df_immune, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Immune.csv'))


# Find the intersection part of the genes----

tumor_neg<-result_df_tumor %>%
  filter(Mean < 0)

tumor_pos<-result_df_tumor %>%
  filter(Mean > 0)


dcis_neg<-result_df_dcis %>%
  filter(Mean < 0)

dcis_pos<-result_df_dcis %>%
  filter(Mean > 0)

immune_neg<-result_df_immune %>%
  filter(Mean < 0)

immune_pos<-result_df_immune %>%
  filter(Mean > 0)

Reduce(intersect, list(dcis_neg$Gene, 
                       tumor_neg$Gene,
                       immune_pos$Gene))
Reduce(intersect, list(dcis_pos$Gene, 
                       tumor_pos$Gene,
                       immune_neg$Gene))

Reduce(intersect, list(dcis_neg$Gene, 
                       tumor_pos$Gene,
                       immune_neg$Gene))


#write.csv(dcis_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/dcis-pos.csv'))
#write.csv(dcis_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/dcis-neg.csv'))


#write.csv(tumor_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/tumor-pos.csv'))
#write.csv(tumor_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/tumor-neg.csv'))


#write.csv(immune_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/immune-pos.csv'))
#write.csv(immune_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/immune-neg.csv'))


##########################################################
# Figure for intersection between Immune cells & DCIS
##########################################################
library(ggplot2)
library(dplyr)
library(ggrepel)

# Merge by Gene to align values
merged_df <- merge(
  result_df_dcis[, c("Gene", "Mean")],
  result_df_immune[, c("Gene", "Mean")],
  by = "Gene",
  suffixes = c("_dcis", "_immune")
)

# Define quadrants and select top 5 genes per quadrant
annot_df <- merged_df %>%
  mutate(quadrant = case_when(
    Mean_dcis > 0 & Mean_immune > 0 ~ "Q1 (+,+)",
    Mean_dcis < 0 & Mean_immune > 0 ~ "Q2 (-,+)",
    Mean_dcis < 0 & Mean_immune < 0 ~ "Q3 (-,-)",
    Mean_dcis > 0 & Mean_immune < 0 ~ "Q4 (+,-)"
  )) %>%
  group_by(quadrant) %>%
  slice_max(order_by = abs(Mean_dcis) + abs(Mean_immune), n = 5, with_ties = FALSE) %>%
  ungroup()




# Merge by Gene to align values
merged_df_2 <- merge(
  result_df_tumor[, c("Gene", "Mean")],
  result_df_immune[, c("Gene", "Mean")],
  by = "Gene",
  suffixes = c("_invasive", "_immune")
)

# Define quadrants and select top 5 genes per quadrant
annot_df_2 <- merged_df_2 %>%
  mutate(quadrant = case_when(
    Mean_invasive > 0 & Mean_immune > 0 ~ "Q1 (+,+)",
    Mean_invasive < 0 & Mean_immune > 0 ~ "Q2 (-,+)",
    Mean_invasive < 0 & Mean_immune < 0 ~ "Q3 (-,-)",
    Mean_invasive > 0 & Mean_immune < 0 ~ "Q4 (+,-)"
  )) %>%
  group_by(quadrant) %>%
  slice_max(order_by = abs(Mean_invasive) + abs(Mean_immune), n = 5, with_ties = FALSE) %>%
  ungroup()





##########################################################
# Scatter plot with annotated and highlighted genes
##########################################################
p_dcis_immune<- ggplot(merged_df, aes(x = Mean_dcis, y = Mean_immune)) +
  # All points (background)
  geom_point(size = 2, alpha = 0.7, color = "gray50") +
  
  # Highlighted extreme points
  geom_point(
    data = annot_df,
    aes(x = Mean_dcis, y = Mean_immune),
    color = "#ed1b2f",
    size = 3
  ) +
  
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  
  # Gene labels for the highlighted points
  geom_text_repel(
    data = annot_df,
    aes(label = Gene),
    size = 5,
    color = "#ed1b2f",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "gray60"
  ) +
  
  # Titles and theme
  labs(
    title = expression(paste("Posterior mean of ", rho, " Between DCIS and Immune")),
    x = "Correlation (DCIS)",
    y = "Correlation (Immune)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20)
  )


p_invasive_immune<- ggplot(merged_df_2, aes(x = Mean_invasive, y = Mean_immune)) +
  # All points (background)
  geom_point(size = 2, alpha = 0.7, color = "gray50") +
  
  # Highlighted extreme points
  geom_point(
    data = annot_df_2,
    aes(x = Mean_invasive, y = Mean_immune),
    color = "#ed1b2f",
    size = 2.8
  ) +
  
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  
  # Gene labels for the highlighted points
  geom_text_repel(
    data = annot_df_2,
    aes(label = Gene),
    size = 5,
    color = "#ed1b2f",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "gray60"
  ) +
  
  # Titles and theme
  labs(
    title = expression(paste("Posterior mean of ", rho, " Between Invasive Tumor and Immune")),
    x = "Correlation (Invasive)",
    y = "Correlation (Immune)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20)
  )


p_dcis_immune 
p_invasive_immune




# Adding uncertainty quantification:##########################################################
# Merge credible intervals for DCIS vs Immune plot
##########################################################
merged_df <- merge(
  result_df_dcis[, c("Gene", "Mean", "Low95", "High95")],
  result_df_immune[, c("Gene", "Mean", "Low95", "High95")],
  by = "Gene",
  suffixes = c("_dcis", "_immune")
)

# Define quadrants and select top 5 genes per quadrant
annot_df <- merged_df %>%
  mutate(quadrant = case_when(
    Mean_dcis > 0 & Mean_immune > 0 ~ "Q1 (+,+)",
    Mean_dcis < 0 & Mean_immune > 0 ~ "Q2 (-,+)",
    Mean_dcis < 0 & Mean_immune < 0 ~ "Q3 (-,-)",
    Mean_dcis > 0 & Mean_immune < 0 ~ "Q4 (+,-)"
  )) %>%
  group_by(quadrant) %>%
  slice_max(order_by = abs(Mean_dcis) + abs(Mean_immune), n = 5, with_ties = FALSE) %>%
  ungroup()

##########################################################
# Merge credible intervals for Invasive vs Immune plot
##########################################################
merged_df_2 <- merge(
  result_df_tumor[, c("Gene", "Mean", "Low95", "High95")],
  result_df_immune[, c("Gene", "Mean", "Low95", "High95")],
  by = "Gene",
  suffixes = c("_invasive", "_immune")
)

# Define quadrants and select top 5 genes per quadrant
annot_df_2 <- merged_df_2 %>%
  mutate(quadrant = case_when(
    Mean_invasive > 0 & Mean_immune > 0 ~ "Q1 (+,+)",
    Mean_invasive < 0 & Mean_immune > 0 ~ "Q2 (-,+)",
    Mean_invasive < 0 & Mean_immune < 0 ~ "Q3 (-,-)",
    Mean_invasive > 0 & Mean_immune < 0 ~ "Q4 (+,-)"
  )) %>%
  group_by(quadrant) %>%
  slice_max(order_by = abs(Mean_invasive) + abs(Mean_immune), n = 5, with_ties = FALSE) %>%
  ungroup()

##########################################################
# Set scaling factor for credible intervals (smaller)
##########################################################
ci_scale <- 0.05  # Smaller scale for more compact error bars

##########################################################
# Scatter plot with scaled credible intervals - DCIS vs Immune
##########################################################
p_dcis_immune <- ggplot(merged_df, aes(x = Mean_dcis, y = Mean_immune)) +
  # Scaled credible interval lines (horizontal for immune) - BLACK
  #geom_segment(
  #  aes(x = Mean_dcis, 
  #      xend = Mean_dcis, 
  #      y = Mean_immune - (Mean_immune - Low95_immune) * ci_scale, 
  #      yend = Mean_immune + (High95_immune - Mean_immune) * ci_scale),
  #  color = "black", linewidth = 0.5, alpha = 0.8
  #) +
  ## Scaled credible interval lines (vertical for dcis) - BLACK
  #geom_segment(
  #  aes(x = Mean_dcis - (Mean_dcis - Low95_dcis) * ci_scale, 
  #      xend = Mean_dcis + (High95_dcis - Mean_dcis) * ci_scale, 
  #      y = Mean_immune, 
  #      yend = Mean_immune),
  #  color = "black", linewidth = 0.5, alpha = 0.8
  #) +
  # All points (background)
  geom_point(size = 4, alpha = 0.7, color = "gray50") +
  
  # Highlighted extreme points
  geom_point(
    data = annot_df,
    aes(x = Mean_dcis, y = Mean_immune),
    color = "#ed1b2f",
    size = 4
  ) +
  
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth=1 ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=1) +
  
  # Gene labels for the highlighted points
  geom_text_repel(
    data = annot_df,
    aes(label = Gene),
    size = 8,
    color = "#ed1b2f",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "gray60"
  ) +
  
  # Titles and theme
  labs(
    #title = expression(paste("Posterior correlation, DCIS & Immune")),
    title='',
    x = "Mean Correlation (DCIS)",
    y = "Mean Correlation (Immune)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5, "pt"),
    strip.text.x = element_text(size = 30),
    strip.text.y = element_text(size = 30)
  )

##########################################################
# Scatter plot with scaled credible intervals - Invasive vs Immune
##########################################################
p_invasive_immune <- ggplot(merged_df_2, aes(x = Mean_invasive, y = Mean_immune)) +
  # Scaled credible interval lines (horizontal for immune) - BLACK
  #geom_segment(
  #  aes(x = Mean_invasive, 
  #      xend = Mean_invasive, 
  #      y = Mean_immune - (Mean_immune - Low95_immune) * ci_scale, 
  #      yend = Mean_immune + (High95_immune - Mean_immune) * ci_scale),
  #  color = "black", linewidth = 0.5, alpha = 0.8
  #) +
  # Scaled credible interval lines (vertical for invasive) - BLACK
  #geom_segment(
  #  aes(x = Mean_invasive - (Mean_invasive - Low95_invasive) * ci_scale, 
  #      xend = Mean_invasive + (High95_invasive - Mean_invasive) * ci_scale, 
  #      y = Mean_immune, 
  #      yend = Mean_immune),
  #  color = "black", linewidth = 0.5, alpha = 0.8
  #) +
  # All points (background)
  geom_point(size = 4, alpha = 0.7, color = "gray50") +
  
  # Highlighted extreme points
  geom_point(
    data = annot_df_2,
    aes(x = Mean_invasive, y = Mean_immune),
    color = "#ed1b2f",
    size = 4
  ) +
  
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth=1 ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=1) +
  
  # Gene labels for the highlighted points
  geom_text_repel(
    data = annot_df_2,
    aes(label = Gene),
    size = 8,
    color = "#ed1b2f",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "gray60"
  ) +
  
  # Titles and theme
  labs(
    #title = expression(paste("Posterior correlation, Invasive Tumor & Immune")),
    title='',
    x = "Mean Correlation (Invasive)",
    y = "Mean Correlation (Immune)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 30),
    strip.text.y = element_text(size = 30)
  )


p_invasive_immune
p_dcis_immune




# Gene of interest ranking-----
dcis_cdh1<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-dcis/CDH1_M2_68.csv'))
dcis_postn<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-dcis/POSTN_M2_231.csv'))
dcis_trac<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-dcis/TRAC_M2_299.csv'))

invasive_cdh1<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-invasive-tumor/CDH1_M2_68.csv'))
invasive_postn<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-invasive-tumor/POSTN_M2_231.csv'))
invasive_trac<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-invasive-tumor/TRAC_M2_299.csv'))

immune_cdh1<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-immune/CDH1_M2_68.csv'))
immune_postn<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-immune/POSTN_M2_231.csv'))
immune_trac<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/results-immune/TRAC_M2_299.csv'))





library(ggplot2)
library(dplyr)
library(rethinking)
set.seed(111)
simulated_priors<-rlkjcorr(5000, 2,eta=1.2)[,1,2]
simulated_Q<- c(quantile(simulated_priors, 0.025),
                mean(simulated_priors),
                quantile(simulated_priors, 0.975))
plot_data <- data.frame(
  tissue = rep(c("DCIS", "Invasive", "Immune"), each = 3),
  marker = rep(c("CDH1", "POSTN", "TRAC"), times = 3),
  mean = c(
    dcis_cdh1$mean[dcis_cdh1$variable == "Astar[2,1]"],
    dcis_postn$mean[dcis_postn$variable == "Astar[2,1]"],
    dcis_trac$mean[dcis_trac$variable == "Astar[2,1]"],
    
    invasive_cdh1$mean[invasive_cdh1$variable == "Astar[2,1]"],
    invasive_postn$mean[invasive_postn$variable == "Astar[2,1]"],
    invasive_trac$mean[invasive_trac$variable == "Astar[2,1]"],
    
    immune_cdh1$mean[immune_cdh1$variable == "Astar[2,1]"],
    immune_postn$mean[immune_postn$variable == "Astar[2,1]"],
    immune_trac$mean[immune_trac$variable == "Astar[2,1]"]
  ),
  lower = c(
    dcis_cdh1$X2.5.[dcis_cdh1$variable == "Astar[2,1]"],
    dcis_postn$X2.5.[dcis_postn$variable == "Astar[2,1]"],
    dcis_trac$X2.5.[dcis_trac$variable == "Astar[2,1]"],
    
    invasive_cdh1$X2.5.[invasive_cdh1$variable == "Astar[2,1]"],
    invasive_postn$X2.5.[invasive_postn$variable == "Astar[2,1]"],
    invasive_trac$X2.5.[invasive_trac$variable == "Astar[2,1]"],
    
    immune_cdh1$X2.5.[immune_cdh1$variable == "Astar[2,1]"],
    immune_postn$X2.5.[immune_postn$variable == "Astar[2,1]"],
    immune_trac$X2.5.[immune_trac$variable == "Astar[2,1]"]
  ),
  upper = c(
    dcis_cdh1$X97.5.[dcis_cdh1$variable == "Astar[2,1]"],
    dcis_postn$X97.5.[dcis_postn$variable == "Astar[2,1]"],
    dcis_trac$X97.5.[dcis_trac$variable == "Astar[2,1]"],
    
    invasive_cdh1$X97.5.[invasive_cdh1$variable == "Astar[2,1]"],
    invasive_postn$X97.5.[invasive_postn$variable == "Astar[2,1]"],
    invasive_trac$X97.5.[invasive_trac$variable == "Astar[2,1]"],
    
    immune_cdh1$X97.5.[immune_cdh1$variable == "Astar[2,1]"],
    immune_postn$X97.5.[immune_postn$variable == "Astar[2,1]"],
    immune_trac$X97.5.[immune_trac$variable == "Astar[2,1]"]
  )
)

plot_data$tissue <- factor(plot_data$tissue, 
                           levels = c("DCIS", "Invasive", "Immune", "Prior"))


### --- Build PRIOR category separately (1 box only) --- ###

prior_data <- data.frame(
  tissue = factor("Prior", levels = c("DCIS", "Invasive", "Immune", "Prior")),
  marker = "Prior",
  mean = simulated_Q[2],
  lower = simulated_Q[1],
  upper = simulated_Q[3]
)


### --- Plot --- ###
pd <- position_dodge(width = 0.6)
ggplot() +
  
  ## COLORED ERRORBARS
  geom_errorbar(
    data = plot_data,
    aes(x = tissue, y = mean, ymin = lower, ymax = upper, color = marker),
    width = 0.2,
    position = pd,
    linewidth = 2
  ) +
  
  ## BLACK POINTS ATTACHED TO EACH ERRORBAR
  geom_point(
    data = plot_data,
    aes(x = tissue, y = mean, group = marker),
    position = pd,
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +
  
  ## PRIOR
  geom_errorbar(
    data = prior_data,
    aes(x = tissue, y = mean, ymin = lower, ymax = upper),
    width = 0.2,
    color = "grey40",
    linewidth = 2
  ) +
  geom_point(
    data = prior_data,
    aes(x = tissue, y = mean),
    size = 4,
    color = "black"
  ) +
  
  geom_abline(
    intercept = 0, slope = 0,
    color = '#ed1b2f',
    linetype = 'dashed',
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = c("CDH1" = "#8B6914",
               "POSTN" = "#00B050",
               "TRAC" = "#8B3A6E"),
    name = "Gene",
    labels = c("CDH1" = expression(italic("CDH1")),
               "POSTN" = expression(italic("POSTN")),
               "TRAC" = expression(italic("TRAC")))
  ) +
  
  ## CUSTOM LEGEND WITH LONGER SEGMENTS
  guides(
    color = guide_legend(
      keywidth = 3,       # longer segments
      keyheight = 1.5,
      override.aes = list(
        linetype = 1,
        linewidth = 2,
        color = c("CDH1" = "#8B6914",
                  "POSTN" = "#00B050",
                  "TRAC" = "#8B3A6E"),
        shape = 19,
        size = 3
      )
    )
  ) +
  
  labs(x = "",
       y = "Posterior of Correlation",
       title = "") +
  
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30, face = "bold"),  # Bold for emphasis
    legend.text = element_text(size = 25, face = "italic"),  # Italic gene names
    legend.position = "top"
  )

# Plot for the gene expression patterns-----

cell_gene <- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv'))
cell_spatial <- read.csv(file.path(DATA_DIR, '10X-breast-cancer/after-preprocess/cell_spatial.csv'))
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000


center_slice_cell_location_DCIS<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
center_slice_cell_location_tumor<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
center_slice_cell_location_immune<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),]


gene_dcis<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
gene_tumor<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]

# check the genes in the dataset:

dcis_df<- cbind(center_slice_cell_location_DCIS,
                gene_dcis$CDH1, 
                gene_dcis$POSTN,
                gene_dcis$TRAC)
colnames(dcis_df)[5]<- 'CDH1'
colnames(dcis_df)[6]<- 'POSTN'
colnames(dcis_df)[7]<- 'TRAC'


tumor_df<- cbind(center_slice_cell_location_tumor,
                gene_tumor$CDH1, 
                gene_tumor$POSTN,
                gene_tumor$TRAC)

colnames(tumor_df)[5]<- 'CDH1'
colnames(tumor_df)[6]<- 'POSTN'
colnames(tumor_df)[7]<- 'TRAC'


immune_df<- cbind(center_slice_cell_location_immune,
                gene_immune$CDH1, 
                gene_immune$POSTN,
                gene_immune$TRAC)

colnames(immune_df)[5]<- 'CDH1'
colnames(immune_df)[6]<- 'POSTN'
colnames(immune_df)[7]<- 'TRAC'


dcis_df$lCDH1<- log(dcis_df$CDH1+0.05)

# DCIS----
ggplot(dcis_df, aes(x = x, y = y, color = lCDH1)) +
  geom_point(size = 0.3, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#B58716",
    #trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "CDH1") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(dcis_df, aes(x = x, y = y, color = POSTN)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#B58716",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "POSTN") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(dcis_df, aes(x = x, y = y, color = TRAC)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#B58716",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "TRAC") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )





# invasive----
ggplot(tumor_df, aes(x = x, y = y, color = CDH1)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#00B050",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "CDH1") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(tumor_df, aes(x = x, y = y, color = POSTN)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#00B050",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "POSTN") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(tumor_df, aes(x = x, y = y, color = TRAC)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#00B050",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "TRAC") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )

# Immune cells
ggplot(immune_df, aes(x = x, y = y, color = CDH1)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#8B3A6E",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "CDH1") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(immune_df, aes(x = x, y = y, color = POSTN)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#8B3A6E",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "POSTN") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )


ggplot(immune_df, aes(x = x, y = y, color = TRAC)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_gradient(
    low = "snow2",
    high = "#8B3A6E",
    trans = "sqrt",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,           # center title (was 1, which is right-aligned)
      barheight = unit(8, "cm"),   # fixed height
      barwidth = unit(0.7, "cm"),  # fixed width
      label.theme = element_text(size = 28),
      title.theme = element_text(size = 30)
    )
  ) +
  labs(title = "", x = "", y = "", color = "TRAC") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),  # distance between plot and legend
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"),  # margin around legend
    axis.text = element_blank()
  )




# Let's check if the grids will help in the visualization-----

##  DCIS dataframe: ----
center_slice_cell_location<- data.frame(cbind(dcis_df$x, dcis_df$y))
colnames(center_slice_cell_location)<- c('x','y')
df<- cbind( center_slice_cell_location ,dcis_df$TRAC)
colnames(df)<- c('x','y','Gene')



grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


# prepare the grids: 
init_x <- min(center_slice_cell_location$x)
max(center_slice_cell_location$x)
init_y <- min(center_slice_cell_location$y)
max(center_slice_cell_location$y)

# the delta: 
step_x <- (max(loc$x)-min(loc$x))/grid_size
step_y <- (max(loc$y)-min(loc$y))/grid_size
#grid_area<- step_x*step_y

# sum up the grids
loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
for (i in 1:grid_size){
  for (j in 1:grid_size){
    x_bound_1 <- init_x + (i-1) * step_x
    x_bound_2 <- init_x + i * step_x
    y_bound_1 <- init_y + (j-1) * step_y
    y_bound_2 <- init_y + j * step_y
    
    pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                     loc$x <= x_bound_2),
                             which(loc$y >= y_bound_1 & 
                                     loc$y <=  y_bound_2))
    marks_100[i,j]<- sum(eg_gene_express[pts_index])
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}
# the log-transformed items 
layout(matrix(1:2, nrow = 1), widths = c(1, 1), respect = F)
# the gene expression 
par(mfrow=c(1,1), mar = c(5, 4, 5, .5))  # small right margin to keep close to left plot
image.plot(
  log(marks_mat+0.05),
  breaks = make_breaks( log(marks_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = '',
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 3,
  cex.main = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(marks_mat), ny = nrow(marks_mat), col = "black", lwd = 0.5)
# Reset
par(mfrow = c(1, 1))




# points 
par(mfrow=c(1,1), mar = c(4, 4, 4, 0.5))  # small right margin
image.plot(
  log(loc_mat+0.05),
  breaks = make_breaks(log(loc_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  cex.main = 3,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)



##  invasive dataframe:----- 
center_slice_cell_location<- data.frame(cbind(tumor_df$x, tumor_df$y))
colnames(center_slice_cell_location)<- c('x','y')
df<- cbind( center_slice_cell_location ,tumor_df$TRAC)
colnames(df)<- c('x','y','Gene')



grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


# prepare the grids: 
init_x <- min(center_slice_cell_location$x)
max(center_slice_cell_location$x)
init_y <- min(center_slice_cell_location$y)
max(center_slice_cell_location$y)

# the delta: 
step_x <- (max(loc$x)-min(loc$x))/grid_size
step_y <- (max(loc$y)-min(loc$y))/grid_size
#grid_area<- step_x*step_y

# sum up the grids
loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
for (i in 1:grid_size){
  for (j in 1:grid_size){
    x_bound_1 <- init_x + (i-1) * step_x
    x_bound_2 <- init_x + i * step_x
    y_bound_1 <- init_y + (j-1) * step_y
    y_bound_2 <- init_y + j * step_y
    
    pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                     loc$x <= x_bound_2),
                             which(loc$y >= y_bound_1 & 
                                     loc$y <=  y_bound_2))
    marks_100[i,j]<- sum(eg_gene_express[pts_index])
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}
# the log-transformed items 
layout(matrix(1:2, nrow = 1), widths = c(1, 1), respect = F)
# the gene expression 
par(mfrow=c(1,1), mar = c(5, 4, 5, .5))  # small right margin to keep close to left plot
image.plot(
  log(marks_mat+0.05),
  breaks = make_breaks( log(marks_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = '',
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 3,
  cex.main = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(marks_mat), ny = nrow(marks_mat), col = "black", lwd = 0.5)
# Reset
par(mfrow = c(1, 1))


# points 
par(mfrow=c(1,1), mar = c(4, 4, 4, 0.5))  # small right margin
image.plot(
  log(loc_mat+0.05),
  breaks = make_breaks(log(loc_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  cex.main = 3,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)



##  immune dataframe:----- 
center_slice_cell_location<- data.frame(cbind(immune_df$x, immune_df$y))
colnames(center_slice_cell_location)<- c('x','y')
df<- cbind( center_slice_cell_location ,immune_df$TRAC)
colnames(df)<- c('x','y','Gene')



grid_size <- 50 # 100*100
len.grid.x <- ( max(center_slice_cell_location$x)- min(center_slice_cell_location$x) )/grid_size
len.grid.y <- ( max(center_slice_cell_location$y)- min(center_slice_cell_location$y) )/grid_size
grid_area <- (len.grid.x * len.grid.y)
grid_points <- expand.grid(x = seq(  min(center_slice_cell_location$x)+len.grid.x ,
                                     max(center_slice_cell_location$x),
                                     length.out=grid_size)-0.5*len.grid.x, 
                           y = seq(min(center_slice_cell_location$y)+len.grid.y ,
                                   max(center_slice_cell_location$y),
                                   length.out=grid_size)-0.5*len.grid.y)
N <- nrow(grid_points)  # Total number of points


loc<- center_slice_cell_location
eg_gene_express<- df$Gene


# prepare the grids: 
init_x <- min(center_slice_cell_location$x)
max(center_slice_cell_location$x)
init_y <- min(center_slice_cell_location$y)
max(center_slice_cell_location$y)

# the delta: 
step_x <- (max(loc$x)-min(loc$x))/grid_size
step_y <- (max(loc$y)-min(loc$y))/grid_size
#grid_area<- step_x*step_y

# sum up the grids
loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
for (i in 1:grid_size){
  for (j in 1:grid_size){
    x_bound_1 <- init_x + (i-1) * step_x
    x_bound_2 <- init_x + i * step_x
    y_bound_1 <- init_y + (j-1) * step_y
    y_bound_2 <- init_y + j * step_y
    
    pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                     loc$x <= x_bound_2),
                             which(loc$y >= y_bound_1 & 
                                     loc$y <=  y_bound_2))
    marks_100[i,j]<- sum(eg_gene_express[pts_index])
    loc_100[i,j] <- length(loc[pts_index,][,1])
  }
}

loc_mat <- loc_100
marks_mat<- marks_100

# Helper to make breaks
make_breaks <- function(mat, n_colors = 30) {
  seq(min(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE),
      length.out = n_colors + 1)
}
# the log-transformed items 
layout(matrix(1:2, nrow = 1), widths = c(1, 1), respect = F)
# the gene expression 
par(mfrow=c(1,1), mar = c(5, 4, 5, .5))  # small right margin to keep close to left plot
image.plot(
  log(marks_mat+0.05),
  breaks = make_breaks( log(marks_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = '',
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 3,
  cex.main = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(marks_mat), ny = nrow(marks_mat), col = "black", lwd = 0.5)
# Reset
par(mfrow = c(1, 1))



# points 
par(mfrow=c(1,1), mar = c(4, 4, 4, 0.5))  # small right margin
image.plot(
  log(loc_mat+0.05),
  breaks = make_breaks(log(loc_mat+0.05), 30),
  col = viridis(30, option='magma'),
  main = "",
  xlab = "",
  ylab = "",
  useRaster = FALSE,
  legend.width = 0.6,
  legend.mar = 2,
  legend.cex = 2,
  legend.args = list(
    text = "",       
    side = 4,                         
    font = 5,                         
    cex = 10,
    cex.axis=2 
  ),
  axis.args = list(
    cex.axis = 2.5
  ),
  cex.main = 3,
  xaxt = "n",
  yaxt = "n"
)
grid(nx = ncol(loc_mat), ny = nrow(loc_mat), col = "black", lwd = 0.5)


