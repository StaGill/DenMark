########################################
# GSEA: gene set enrichment analysis 
# On the breast cancer dataset 
########################################


# 1. load the package 
library(org.Hs.eg.db)
library('DOSE')
library('clusterProfiler')
library(msigdbr)
library(fgsea)
library(enrichplot)

# 2. load the DenMark results from the breast cancer dataset: this file only show for the invasive tumor cells
gene_corr_DenMark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Tumor.csv')


# Make named numeric vector
gene_list <- gene_corr_DenMark$Mean
names(gene_list) <- gene_corr_DenMark$Gene

# Sort
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert SYMBOL -> ENSEMBL
gene_info <- bitr(names(gene_list),
                  fromType = "SYMBOL",
                  toType = "ENSEMBL",
                  OrgDb = org.Hs.eg.db)

# Keep only successfully mapped genes
gene_list <- gene_list[names(gene_list) %in% gene_info$SYMBOL]

# Replace names with Ensembl IDs
names(gene_list) <- gene_info$ENSEMBL[match(names(gene_list), gene_info$SYMBOL)]
# 3. do the enrichment analysis based on the tutorial GSEA: 
# Run GESA:
gse <- gseGO(
  geneList = gene_list,
  ont = "ALL",                 # BP, MF, CC, or ALL
  keyType = "ENSEMBL",         # or "SYMBOL" depending on what you have
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,              # filter very small GO terms
  maxGSSize = 1000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",        # use BH correction instead of "none"
  verbose = TRUE,
  eps = 0                      # prevents numerical issues
)


# check the output and visualize the results: 
dim(gse@result)
head(gse@result)


dotplot(gse, showCategory = 10)
ridgeplot(gse)
gseaplot2(gse, geneSetID = gse@result$ID[1])

ridgeplot(gse) + labs(x = "enrichment distribution")
# from tutorials: 
dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)+
  theme(strip.text = element_text(size = 16, face = "bold"))




df <- gse@result

# Keep only top 10 terms (or as you like)
#df <- df[1:10, ]

df$direction <- ifelse(df$NES > 0, "activation", "suppression")

# Make a dotplot with shape for activation/suppression
ggplot(df, aes(x = reorder(Description, NES), 
               y = NES, 
               size = abs(NES),         # optional: size proportional to |NES|
               shape = direction,
               color = p.adjust)) +
  geom_point() +
  coord_flip() +                # horizontal bars
  scale_shape_manual(values = c(16, 17)) +  # circle vs triangle
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "", y = "Normalized Enrichment Score (NES)", 
       size = "|NES|", color = "Adjusted P") +
  theme_bw(base_size = 16) +    # larger text
  theme(legend.position = "right")


ggplot(df, aes(x = reorder(Description, NES), 
               y = NES, 
               size = abs(NES),       # optional: effect size
               shape = direction,     # binary shape for activation/suppression
               color = p.adjust)) +   # optional: color by adjusted p-value
  geom_point() +
  coord_flip() +                        # horizontal layout
  scale_shape_manual(values = c(16, 17)) +  # circle = activation, triangle = suppression
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "", y = "Normalized Enrichment Score (NES)", 
       size = "|NES|", color = "Adjusted P") +
  theme_bw(base_size = 16) +            # increase text size
  theme(legend.position = "right")





# compute pairwise term similarity
gse_imp <- pairwise_termsim(gse)

# now plot
emapplot(gse_imp, showCategory = 10)


cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)


cnetplot(
  gse,
  showCategory = 5,
  categorySize = "geneNum",
  color.params = list(foldChange = gene_list)
)


library(enrichplot)
library(ggrepel)


cnetplot(
  gse,
  showCategory = 5,
  node_label = "category",
  layout = "kk",
  color.params = list(foldChange = gene_list),
  cex.params = list(category_label = 0.6)
)





p <- cnetplot(
  gse,
  showCategory = 5,                                # top 5 pathways
  node_label = "all",                              # label both categories & genes
  layout = "kk",                                   # Kamada-Kawai layout
  color.params = list(foldChange = gene_list),
  cex.params = list(category_label = 0.7, gene_label = 0.6)  # adjust font size
)

# Optional: tweak text size further and export high-resolution
p + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 9),
  plot.title = element_text(size = 10, face = "bold")
)


# GSEA plot: 
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)




terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)



