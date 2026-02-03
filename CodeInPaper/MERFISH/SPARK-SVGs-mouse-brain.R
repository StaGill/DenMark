# Real data analysis: Identifying SVGs with SPARK
# Aim: to identify SVGs in the ST datasets, and Venn plots in comparison to the DenMark 

# 1. First, we install SPARK package: 
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


library(doParallel)
registerDoParallel(cores = 5)

load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-location.Rdata')
load('C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/clusters-computecanada/Dataset-Analysis/Astrocyte-DenMark/Dataset/cell-gene-expression.Rdata')

# locations: 23624 cell locations, with x and y coordinates 
str(center_slice_cell_location)
# gene expression, 649 genes with 23624 cell locations 
str(center_slice_gene_location)


# 2. load the dataset, and identify the SVGs with SPARK in the dataset 


# remove all the blank label genes 
raw_count_mat <- as.matrix(center_slice_gene_location)[1:483,]

# only keep the genes expressed in at least 5% of the cells 
#keep_genes <- rowSums(raw_count_mat > 0) > 0.05 * ncol(raw_count_mat)
#raw_count_mat <- raw_count_mat[keep_genes, ]

rownames(raw_count_mat)




# use spark-X: 

sparkX <- sparkx(raw_count_mat,
                 center_slice_cell_location,
                 numCores=1,option="mixture")
head(sparkX$res_mtest)
order(sparkX$res_mtest$adjustedPval,decreasing=T)


which(sparkX$res_mtest$adjustedPval<0.05)

SVGs005<-rownames(raw_count_mat)[which(sparkX$res_mtest$adjustedPval<0.05)]


# 350 SVGs
write.csv(as.data.frame(SVGs005, ncol=1),
          file='C:/Users/xumc7/OneDrive/Desktop/PhD Thesis 1/final-summarize-codes/Other-Methods-Comparison/SPARKSVG.csv' )




## use spark: but this is too slow (>6 hrs)
# Note that this cell type is 23624 of 483 genes 
#spark <- CreateSPARKObject(counts=raw_count_mat, 
#                           location=center_slice_cell_location[,1:2],
#                           percentage = 0.1, 
#                           min_total_counts = 10)
#spark@lib_size <- apply(spark@counts, 2, sum)
#
#rownames(spark@counts)
#
## Estimating Parameter Under Null
#spark <- spark.vc(spark, 
#                  covariates = NULL, 
#                  lib_size = spark@lib_size, 
#                  num_core = 5,
#                  verbose = F)
## Calculating pval
#spark.pval <- spark.test(spark, 
#                    check_positive = T, 
#                    verbose = F)
#head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])

























