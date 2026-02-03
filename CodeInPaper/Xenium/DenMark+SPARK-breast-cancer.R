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
cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/cell_spatial.csv')
# 0A. translation:-----
cell_spatial$x<- cell_spatial$x/1000
cell_spatial$y<- cell_spatial$y/1000

# check the genes in the dataset:
colnames(cell_gene)

# A. on dcis region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]
gene_spatial_DCIS<- cell_gene[which(cell_spatial$celltype=='DCIS_1'|cell_spatial$celltype=='DCIS_2'),]

raw_count_mat <- as.matrix(gene_spatial_DCIS)[,-1]

gene_name_dcis<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/DCIS-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-DCIS.csv' )



# B. on invasive region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='Invasive_Tumor'),]
gene_spatial_invasive<- cell_gene[which(cell_spatial$celltype=='Invasive_Tumor'),]

raw_count_mat <- as.matrix(gene_spatial_invasive)[,-1]

gene_name_invasive<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Tumor-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Tumor.csv' )


# C. on immune cells region ----
center_slice_cell_location<- cell_spatial[which(cell_spatial$celltype=='CD8+_T_Cells'),]
gene_spatial_immune<- cell_gene[which(cell_spatial$celltype=='CD8+_T_Cells'),]

raw_count_mat <- as.matrix(gene_spatial_immune)[,-1]

gene_name_immune<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/Immune-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Immune.csv' )




# The intersection between spark and denmark identified genes -----
df_dcis_spark    <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-DCIS.csv')
df_invasive_spark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Tumor.csv')
df_immune_spark  <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/SPARKSVG-Immune.csv')


# order the p-value in an increasing rank 
df_dcis_order <- df_dcis_spark %>% 
  arrange(adj_pval)

df_invasive_order <- df_invasive_spark %>% 
  arrange(adj_pval)

df_immune_order <- df_immune_spark %>% 
  arrange(adj_pval)


# compare that with the denmark results: 
df_dcis_denmark   <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-DCIS.csv')
df_invasive_denmark   <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Tumor.csv')
df_immune_denmark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/large-range-genes/DCG-Immune.csv')

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
cell_gene <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv')
cell_spatial <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/10X-breast-cancer/after-preprocess/cell_spatial.csv')
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

gene_name_dcis1<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-I-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS1.csv' )







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


gene_name_dcis2<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/DCIS-II-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-DCIS2.csv' )




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


gene_name_invasive<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/Invasive-after-QC.csv')$x
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
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-invasive.csv' )




# the ranking of genes from SPARK and from DennMark ----
df_dcis1_spark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-dcis1.csv')
df_dcis2_spark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-dcis2.csv')
df_invasive_spark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/SPARKSVG-invasive.csv')


# order the p-value in an increasing rank 
df_dcis1_order <- df_dcis1_spark %>% 
  arrange(adj_pval)

df_dcis2_order <- df_dcis2_spark %>% 
  arrange(adj_pval)

df_invasive_order <- df_invasive_spark %>% 
  arrange(adj_pval)


# compare that with the denmark results: 
df_dcis1_denmark   <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS1.csv')
df_dcis2_denmark   <- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-DCIS2.csv')
df_invasive_denmark<- read.csv('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/Cell-report-summary-codes-results/DenMark+SPARK/Xenium/short-range-genes/DCG-invasive.csv')

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



