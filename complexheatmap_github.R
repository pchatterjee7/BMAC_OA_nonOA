#####GSEA HEATMAP COMPLEXHEATMAP##########
library(dplyr)
library(ComplexHeatmap)
setwd("~/final_fig")
nci <- read.csv('top10_transformed.csv')
rownames(nci) <- nci$CellType
nci <- select(nci, -c('CellType'))
#col_fun = circlize::colorRamp2(c(-8, 0, 8), c("blue", "white", "red"))
col_fun = circlize::colorRamp2(c(-8, -4, 0, 4, 8), c("blue", "mediumpurple2", "white", "lightsalmon","red"))
lgd = Legend(col_fun = col_fun, title = "GSEA_NES")
Heatmap(nci, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

## get the data matrix with pathways and scores
my_mat <- as.matrix(nci) 

Heatmap(my_mat, name = "GSEA_NES", cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2), column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10))




nci <- read.csv('top10_transformed.csv')
rownames(nci) <- nci$CellType
nci <- nci[, -1] # Remove the CellType column

col_fun = circlize::colorRamp2(c(-8, -4, 0, 4, 8), c("blue", "mediumpurple2", "white", "lightsalmon","red"))
lgd = Legend(col_fun = col_fun, title = "GSEA_NES")

my_mat <- as.matrix(nci)

Heatmap(my_mat, name = "GSEA_NES", cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2), column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10))