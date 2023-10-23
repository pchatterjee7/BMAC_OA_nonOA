# large dataset merge and integrate #####

############ too large to merge################

# Define the names of all the Seurat objects you want to include in the list
seurat_names <- c("sample1", "sample2", "sample3") #add all the sample names in order

# Create a list of Seurat objects based on their names
seurat_list <- mget(seurat_names)

setwd("~/dir")

saveRDS(seurat_list, file = "object_final.list.rds")
object_final.list <- readRDS("~/dir/object_final.list.rds")

# Create a list of Seurat objects (assuming you have already run the 'mget()' command)
# seurat_list <- mget(seurat_names)

# Step 1: Perform QC and normalization on each Seurat object
seurat_list_qc <- lapply(object_final.list, function(x) {
  # QC based on mitochondrial gene expression or other criteria (example here)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  
  # Add log10GenesPerUMI to metadata
  x[["log10GenesPerUMI"]] <- log10(x[["nFeature_RNA"]]) / log10(x[["nCount_RNA"]])
  
  # Filter cells based on QC metrics (example here)
  x <- subset(x, subset = nFeature_RNA >= 250 & nFeature_RNA <= 2500 & percent.mt <= 20 & nCount_RNA >= 500)
  
  # Normalize the data
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  
  return(x)
})


# Assuming that each Seurat object has metadata column "Site"
all_data <- merge(x = seurat_list_qc[[1]], y = seurat_list_qc[-1])




counts <- GetAssayData(object = all_data, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
all_data <- CreateSeuratObject(filtered_counts, meta.data = all_data@meta.data)



metadata <- all_data@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


# Create sample column 
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^incmplt_"))] <- "incmplt"
metadata$sample[which(str_detect(metadata$cells, "^control_"))] <- "ctrl"



qc_plot <-  VlnPlot(all_data, features = c("nUMI",
                                              "percent.mt",
                                              "nGene"), ncol = 3, split.by = "Sample")


# Rename columns
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

all_data@meta.data <- metadata

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=Sample, fill=Sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") + NoLegend()




# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell cycle scores
all_data <- CellCycleScoring(all_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

saveRDS(all_data, file = "filtered_obj.rds")
# DimPlot(bmac_phase,
#                 reduction = "pca",
#                  group.by= "Phase",
#                  split.by = "Phase")



# Step 2: Re-split the samples after merging by the column to correct on
split_by_Site <- SplitObject(all_data, split.by = "Site")

# Step 3: Run the integration with RPCA

split_seurat <- lapply(X = split_by_Site, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = split_seurat)
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

saveRDS(immune.combined, file = "object_integrated_rpca.rds")

