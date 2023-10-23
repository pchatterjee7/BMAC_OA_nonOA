### monocle3 bmac
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)

bmac <- readRDS("~/MILES/bmac_noerythrocytes_renamed.rds")
OA <- subset(bmac, subset = Condition == "OA")
oa.cds <- as.cell_data_set(OA)
oa.cds <- cluster_cells(cds = oa.cds, reduction_method = "UMAP")
oa.cds <- learn_graph(oa.cds, use_partition = TRUE)

nonOA <- subset(bmac, subset = Condition == "non OA")
nonoa.cds <- as.cell_data_set(nonOA)
nonoa.cds <- cluster_cells(cds = nonoa.cds, reduction_method = "UMAP")
nonoa.cds <- learn_graph(nonoa.cds, use_partition = TRUE)


#To compute pseudotime estimates for each trajectory we need to decide what the start of each trajectory is. In our case, we know that the hematopoietic stem cells are the progenitors of other cell types in the trajectory, so we can set these cells as the root of the trajectory. Monocle 3 includes an interactive function to select cells as the root nodes in the graph. This function will be launched if calling order_cells() without specifying the root_cells parameter. Here we've pre-selected some cells as the root, and saved these to a file for reproducibility. This file can be downloaded here.

# load the pre-selected HSCs
hsc <- readLines("../vignette_data/hsc_cells.txt")
## for miles I chose HSPcs as a root for OA and nonOA subsets
## omit root_cells = hsc
oa.cds <- order_cells(oa.cds, reduction_method = "UMAP")
nonoa.cds <- order_cells(nonoa.cds, reduction_method = "UMAP")


plot_cells(
  cds = oa.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)


plot_cells(
  cds = nonoa.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

bmac <- AddMetaData(
  object = bmac,
  metadata = oa.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "oa"
)

bmac <- AddMetaData(
  object = bmac,
  metadata = nonoa.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "nonoa"
)

FeaturePlot(bmac, c("oa", "nonoa"), pt.size = 0.2) & scale_color_viridis_c()

OAcds_subset <- choose_cells(oa.cds)
subset_pr_test_res <- graph_test(OAcds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))




#### find the root cells programmatically ####

# find all possible partitions
all_partitions <- unique(cds@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]

# set all partitions to 1
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions %in% all_partitions] <- "1"

get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))