#####################################################
# For cell type annotation in two dataset analysis  #
#####################################################


library("sf")
library(rhdf5)
library(rgeos)
library('ggplot2')
library(nimble)
library('rstan')
library('MASS')
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(cmdstanr)
# 1. the cell type annotation for MERFISH dataset -----
#
# Step 1: load the dataset 

### read in the list of spatial transcriptomics dataset across multiple tissue slices
sp_input = readRDS("C:/Users/xumc7/OneDrive/Desktop/IRIS-example/countList_mouseBrain_Vizgen.RDS")
#### extract the count data list and location list
spatial_countMat_list = sp_input$spatial_countMat_list
spatial_location_list = sp_input$spatial_location_list




# the gene expressions: 
# 1. marker genes for excitatory neurons: Slc17a6, Slc17a7 (Spatially patterned excitatory neuron
#subtypes and projections of the claustrum)
# 2. marker genes for inhibitory neurons: Gad1, Slc32a1
# 3. marker genes for astrocyte cells: Aqp4, Gfap, Aldh1l1
# 4. marker genes for microglia cells: Cx3cr1, P2ry12


marker_gene_df<-cbind(spatial_location_list[[1]]$x,
                      spatial_location_list[[1]]$y,
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Slc17a6'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Slc17a7'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Gad1'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Slc32a1'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Aqp4'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Aldh1l1'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='P2ry12'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Gfap'),],
                      spatial_countMat_list[[1]][which(rownames(spatial_countMat_list[[1]])=='Cx3cr1'),]
)


colnames(marker_gene_df)<- c('x','y','Slc17a6',
                             'Slc17a7','Gad1','Slc32a1','Aqp4','Aldh1l1',
                             'P2ry12', 'Gfap','Cx3cr1')
marker_gene_df<- as.data.frame(marker_gene_df)


# test the coverage: 
t_cells<- marker_gene_df[which(marker_gene_df$Slc17a6 >=1 |
                               marker_gene_df$Slc17a7>=1|
                               marker_gene_df$Gad1>=1|
                                 marker_gene_df$Slc32a1>=1| 
                                 marker_gene_df$Aqp4>=1|
                                 marker_gene_df$Gfap>=1|
                                 marker_gene_df$Aldh1l1>=1|
                                 marker_gene_df$Cx3cr1>=1|
                                 marker_gene_df$P2ry12>=1) ,]



library(cowplot)

p00<- ggplot(marker_gene_df, aes(x=x, y=y))+
  geom_point(size=.45)+
  ggtitle('Whole Region')+
  theme_classic()

ptest<-ggplot(t_cells, aes(x=x, y=y))+
  geom_point(size=.45)+
  ggtitle('Test Cells')+
  theme_classic()
plot_grid(p00, ptest,  ncol=2)

# the coverage of the cell types: 
celltype_cover<-1-(length(marker_gene_df$x)-length(t_cells$x))/(length(marker_gene_df$x))
celltype_cover # already cover 95.6% of the cell types 
# it suffices to visualize these 4 cell types with the other one type as 'Other cell type'
# cell type marking: 
# --- Cell type assignment ---
astrocyte_cells <- marker_gene_df[
  ( (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Gfap >= 1) |
      (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Aldh1l1 >= 1) |
      (marker_gene_df$Gfap >= 1 & marker_gene_df$Aldh1l1 >= 1) ),
]

microglia_cells <- marker_gene_df[
  (marker_gene_df$P2ry12 >= 1 & marker_gene_df$Cx3cr1 >= 1),
]

neuron_cells <- marker_gene_df[
  (marker_gene_df$Slc17a6 >= 1 | marker_gene_df$Slc17a7 >= 1 |
     marker_gene_df$Gad1 >= 1 | marker_gene_df$Slc32a1 >= 1),
]

# Label all cells
marker_gene_df$cell_type <- "Other"
marker_gene_df$cell_type[
  (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Gfap >= 1) |
    (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Aldh1l1 >= 1) |
    (marker_gene_df$Gfap >= 1 & marker_gene_df$Aldh1l1 >= 1)
] <- "Astrocyte"

marker_gene_df$cell_type[
  (marker_gene_df$P2ry12 >= 1 & marker_gene_df$Cx3cr1 >= 1)
] <- "Microglia"

marker_gene_df$cell_type[
  (marker_gene_df$Slc17a6 >= 1 | marker_gene_df$Slc17a7 >= 1 |
     marker_gene_df$Gad1 >= 1 | marker_gene_df$Slc32a1 >= 1)
] <- "Neuron"



library(ggplot2)

ggplot(marker_gene_df, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c(
    "Neuron" = "#1f78b4",
    "Astrocyte" = "#33a02c",
    "Microglia" = "#e31a1c",
    "Other" = "gray80"
  )) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(color = "Cell type", x = "x", y = "y", title = "MERFISH cell-type annotation")


p_compare<- ggplot(astrocyte_cells, aes(x=x, y=y))+
  geom_point(size = 0.6, alpha = 0.8, color='red')





# check by intersection:  
astro_idx <- which( (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Gfap >= 1) |
                      (marker_gene_df$Aqp4 >= 1 & marker_gene_df$Aldh1l1 >= 1) |
                      (marker_gene_df$Gfap >= 1 & marker_gene_df$Aldh1l1 >= 1) )

microglia_idx <- which( (marker_gene_df$P2ry12 >= 1) |
                          (marker_gene_df$Cx3cr1 >= 1) )

neuron_idx <- which( (marker_gene_df$Slc17a6 >= 1) | 
                       (marker_gene_df$Slc17a7 >= 1) |
                       (marker_gene_df$Gad1 >= 1) |
                       (marker_gene_df$Slc32a1 >= 1) )


# Pairwise intersections
intersect_astro_micro <- intersect(astro_idx, microglia_idx)
intersect_astro_neuron <- intersect(astro_idx, neuron_idx)
intersect_micro_neuron <- intersect(microglia_idx, neuron_idx)


cat("Astrocyte & Microglia overlap:", length(intersect_astro_micro), "\n")
cat("Astrocyte & Neuron overlap:", length(intersect_astro_neuron), "\n")
cat("Microglia & Neuron overlap:", length(intersect_micro_neuron), "\n")

marker_gene_df$cell_type <- "Other"

# Assign hierarchy
marker_gene_df$cell_type[neuron_idx] <- "Neuron"
marker_gene_df$cell_type[microglia_idx] <- "Microglia"
marker_gene_df$cell_type[astro_idx] <- "Astrocyte"

library(ggplot2)

p_whole<- ggplot(marker_gene_df, aes(x = x/1000, y = y/1000, color = cell_type)) +
  geom_point(size = 1.5, alpha = 0.8) +
  labs(title = '', x = 'X coordinate (mm)', y = 'Y coordinate (mm)', color = 'Cell Type') +   
  scale_color_manual(values = c("Astrocyte" = "#E2A7CC",
                                "Microglia" = "#C8EAF5",
                                "Neuron"   = "#FFF193",
                                "Other"    = "#D4D4D4")) +
  #theme_void() +
  theme_classic()+
  theme(
    text = element_text(size = 40),
    legend.title = element_text(size = 30, margin = margin(b = 30)),
    legend.text  = element_text(size = 30),
    legend.margin = margin(t = 10, b = 10),
    legend.spacing.y = unit(0.5, "cm"),
    strip.text.x = element_text(size = 40),
    strip.text.y = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(size = 14)))



plot_grid(p_whole, p_compare, nrow=1, align='v')

p_whole




# 2. the cell type annotation for 10x Xenium dataset ----

# 10X Xenium provides the cell type annotation of the cell types 













# Appendix for checking of the codes ------

# test some gene expression: -----
ggplot(marker_gene_df, aes(x = x, y = y)) +
  geom_point(aes(color = P2ry12), size = 1, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", trans = "sqrt") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(title = "Expression of P2ry12 across spatial locations",
       color = "P2ry12 (counts)") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

ggplot(marker_gene_df, aes(x = x, y = y)) +
  geom_point(aes(color = Cx3cr1), size = 1, alpha = 0.9) +
  scale_color_viridis_c(option = "magma", trans = "sqrt") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(title = "Expression of Cx3cr1 across spatial locations",
       color = "Cx3cr1 (counts)") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

microglia_cells <- subset(marker_gene_df, cell_type == "Microglia")

ggplot(microglia_cells, aes(x = x, y = y)) +
  geom_point(aes(color = P2ry12 + Cx3cr1), size = 1, alpha = 0.9) +
  scale_color_viridis_c(option = "inferno", trans = "sqrt") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(title = "Combined expression of P2ry12 and Cx3cr1 in Microglia",
       color = "Expression sum") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())














