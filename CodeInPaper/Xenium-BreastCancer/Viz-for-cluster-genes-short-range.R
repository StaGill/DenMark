## ---- CONFIGURE ME ------------------------------------------------------------
DATA_DIR   <- "data"     # per-gene ROI fits and QC gene lists, under their original sub-folders (see README.md)
OUTPUT_DIR <- "output"   # ROI DCG tables and gene lists written here; share it with the other scripts of this folder
REPO_ROOT  <- "."        # root of the DenMark checkout
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "Cell-report-summary-codes-results/dataset-multiple-genes-Xenium"), recursive = TRUE, showWarnings = FALSE)
## -----------------------------------------------------------------------------

#############################
# viz for BC short range analysis
# gene plots 
#############################

##########################
# DCIS-I -----
##########################

library(readxl)
library(cowplot)
#gene_tab<-read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv'))

#gene_name_ls<-colnames(gene_tab)[-1]


# all gene names: 

# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/short-range/DCIS-I")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)

# Extract gene names (everything before "_M2_")
#gene_name_ls <- sub("_M2_.*", "", basename(files))

# gene name after filtering: 
gene_name_ls<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-I-after-QC.csv'))$x

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
result_df_dcis1 <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)

colnames(result_df_dcis1)<- c('Gene','Mean','Low95','High95')

write.csv(result_df_dcis1, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS1.csv'))


##########################
# DCIS-2 -----
##########################


# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/short-range/DCIS-II")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)


# gene name after filtering: 
gene_name_ls<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-II-after-QC.csv'))$x

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
result_df_dcis2 <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)

colnames(result_df_dcis2)<- c('Gene','Mean','Low95','High95')

write.csv(result_df_dcis2, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS2.csv'))


#######################################
# for invasive cells
#######################################


# Set the path of dcis
path <- file.path(DATA_DIR, "Cell-report-summary-codes-results/cluster-BC/short-range/Invasive")

# Get all file names that match your pattern (e.g., [GeneName]_M2_1.csv etc.)
files <- list.files(path, pattern = "_M2_\\d+\\.csv$", full.names = TRUE)


# gene name after filtering: 
gene_name_ls<- read.csv(file.path(DATA_DIR, 'Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/Invasive-after-QC.csv'))$x

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
result_df_invasive <- data.frame(
  gene = ordered_genes,
  mean = mean_vec,
  low95 = low95_vec,
  high95 = high95_vec
)
colnames(result_df_invasive)<- c('Gene','Mean','Low95','High95')


write.csv(result_df_invasive, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-invasive.csv'))

##########################################################
# Figure for valcano plots
##########################################################



library(ggplot2)
library(dplyr)
library(ggrepel)

# Calculate credible interval length
result_df_invasive <- result_df_invasive %>%
  mutate(CI_length = High95 - Low95)

# Calculate median (50th percentile) of CI length
ci_threshold <- quantile(result_df_invasive$CI_length, 0.50)

# Identify top 5 genes with large positive correlation and LOW CI (below threshold)
top_positive <- result_df_invasive %>%
  filter(Mean > 0, CI_length < ci_threshold) %>%
  arrange(desc(Mean), CI_length) %>%
  head(5)

# Identify top 5 genes with large negative correlation and LOW CI (below threshold)
top_negative <- result_df_invasive %>%
  filter(Mean < 0, CI_length < ci_threshold) %>%
  arrange(Mean, CI_length) %>%
  head(5)

# Add color labels
result_df_invasive <- result_df_invasive %>%
  mutate(color_group = case_when(
    Gene %in% top_positive$Gene ~ "Positive",
    Gene %in% top_negative$Gene ~ "Negative",
    TRUE ~ "Other"
  ))

# Create the plot
p_invasive <- ggplot(result_df_invasive, aes(x = Mean, y = CI_length, color = color_group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = ci_threshold, 
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(data = filter(result_df_invasive, color_group != "Other"),
                  aes(label = Gene),
                  size = 10,
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  min.segment.length = 0,
                  force = 2) +
  scale_color_manual(values = c("Positive" = '#3C6F67', 
                                "Negative" = "#732F00", 
                                "Other" = "gray70")) +
  labs(x = "Mean Correlation", 
       y = "Credible Interval Length",
       title = "Invasive Tumor",
       color = "Gene Group") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 30, face = "bold"),
    legend.position = 'none'
  )

print(p_invasive)


########################
# for dcis-1:----
########################

# Calculate credible interval length
result_df_dcis1 <- result_df_dcis1 %>%
  mutate(CI_length = High95 - Low95)

# Calculate median (50th percentile) of CI length
ci_threshold <- quantile(result_df_dcis1$CI_length, 0.50)

# Identify top 5 genes with large positive correlation and LOW CI (below threshold)
top_positive <- result_df_dcis1 %>%
  filter(Mean > 0, CI_length < ci_threshold) %>%
  arrange(desc(Mean), CI_length) %>%
  head(5)

# Identify top 5 genes with large negative correlation and LOW CI (below threshold)
top_negative <- result_df_dcis1 %>%
  filter(Mean < 0, CI_length < ci_threshold) %>%
  arrange(Mean, CI_length) %>%
  head(5)

# Add color labels
result_df_dcis1 <- result_df_dcis1 %>%
  mutate(color_group = case_when(
    Gene %in% top_positive$Gene ~ "Positive",
    Gene %in% top_negative$Gene ~ "Negative",
    TRUE ~ "Other"
  ))

# Create the plot
p_dcis1 <- ggplot(result_df_dcis1, aes(x = Mean, y = CI_length, color = color_group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = ci_threshold, 
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(data = filter(result_df_dcis1, color_group != "Other"),
                  aes(label = Gene),
                  size = 10,
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  min.segment.length = 0,
                  force = 2) +
  scale_color_manual(values = c("Positive" = '#3C6F67', 
                                "Negative" = "#732F00", 
                                "Other" = "gray70")) +
  labs(x = "Mean Correlation", 
       y = "Credible Interval Length",
       title = "DCIS 1",
       color = "Gene Group") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 30, face = "bold"),
    legend.position = 'none'
  )

print(p_dcis1)


########################
# for dcis-2:----
########################

# Calculate credible interval length
result_df_dcis2 <- result_df_dcis2 %>%
  mutate(CI_length = High95 - Low95)

# Calculate median (50th percentile) of CI length
ci_threshold <- quantile(result_df_dcis2$CI_length, 0.50)

# Identify top 5 genes with large positive correlation and LOW CI (below threshold)
top_positive <- result_df_dcis2 %>%
  filter(Mean > 0, CI_length < ci_threshold) %>%
  arrange(desc(Mean), CI_length) %>%
  head(5)

# Identify top 5 genes with large negative correlation and LOW CI (below threshold)
top_negative <- result_df_dcis2 %>%
  filter(Mean < 0, CI_length < ci_threshold) %>%
  arrange(Mean, CI_length) %>%
  head(5)

# Add color labels
result_df_dcis2 <- result_df_dcis2 %>%
  mutate(color_group = case_when(
    Gene %in% top_positive$Gene ~ "Positive",
    Gene %in% top_negative$Gene ~ "Negative",
    TRUE ~ "Other"
  ))

# Create the plot
p_dcis2 <- ggplot(result_df_dcis2, aes(x = Mean, y = CI_length, color = color_group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = ci_threshold, 
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(data = filter(result_df_dcis2, color_group != "Other"),
                  aes(label = Gene),
                  size = 9,
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  min.segment.length = 0,
                  force = 2) +
  scale_color_manual(values = c("Positive" = '#3C6F67', 
                                "Negative" = "#732F00", 
                                "Other" = "gray70")) +
  labs(x = "Mean Correlation", 
       y = "Credible Interval Length",
       title = "DCIS 2",
       color = "Gene Group") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 30, face = "bold"),
    legend.position = 'none'
  )

print(p_dcis2)

invasive_neg<-result_df_invasive %>%
  filter(Mean < 0)

invasive_pos<-result_df_invasive %>%
  filter(Mean > 0)


dcis1_neg<-result_df_dcis1 %>%
  filter(Mean < 0)

dcis1_pos<-result_df_dcis1 %>%
  filter(Mean > 0)

dcis2_neg<-result_df_dcis2 %>%
  filter(Mean < 0)

dcis2_pos<-result_df_dcis2 %>%
  filter(Mean > 0)



Reduce(intersect, list(invasive_neg$Gene, 
                       dcis1_neg$Gene,
                       dcis2_neg$Gene))

Reduce(intersect, list(invasive_pos$Gene, 
                       dcis1_pos$Gene,
                       dcis2_pos$Gene))


length(c(Reduce(intersect, list(dcis1_neg$Gene,
                       dcis2_neg$Gene)),Reduce(intersect, list(dcis1_pos$Gene,
                       dcis2_pos$Gene))))



length(c(Reduce(intersect, list(dcis1_neg$Gene,
                                invasive_neg$Gene)),Reduce(intersect, list(dcis1_pos$Gene,
                                                                        invasive_pos$Gene))))


length(c(Reduce(intersect, list(dcis2_neg$Gene,
                                invasive_neg$Gene)),Reduce(intersect, list(dcis2_pos$Gene,
                                                                           invasive_pos$Gene))))

# save the gene name list: 

dcis_1_gene<- Reduce(union, list(dcis1_pos$Gene,
                       dcis1_neg$Gene))


dcis_2_gene<- Reduce(union, list(dcis2_pos$Gene,
                                 dcis2_neg$Gene))


invasive_gene<- Reduce(union, list(invasive_pos$Gene,
                                 invasive_neg$Gene))


write.csv(dcis1_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis1-pos.csv'))
write.csv(dcis1_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis1-neg.csv'))


write.csv(dcis2_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis2-pos.csv'))
write.csv(dcis2_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/dcis2-neg.csv'))


write.csv(invasive_pos$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/invasive-pos.csv'))
write.csv(invasive_neg$Gene, file=file.path(OUTPUT_DIR, 'Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/invasive-neg.csv'))





