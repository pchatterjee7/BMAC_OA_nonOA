# ---- 1. Data Processing ----
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
filtered_data <- subset(data, subset = nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 50)

# ---- 2. Visualization ----
VlnPlot(filtered_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ---- 3. Integration ----
data_list <- SplitObject(filtered_data, split.by = "Condition")
data_list <- lapply(data_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data_list)
anchors <- FindIntegrationAnchors(object.list = data_list, anchor.features = features, reduction = "rpca", k.anchor = 20)
integrated_data <- IntegrateData(anchorset = anchors)
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)

# ---- 4. Cell Annotation ----
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5) # Adjust resolution based on the granularity you want.
# Assign identities based on known marker genes
# (This is a hypothetical list; replace with actual marker genes for your experiment)
marker_genes <- list(
  Neurons = c("GeneA", "GeneB"),
  Astrocytes = c("GeneC", "GeneD"),
  Microglia = c("GeneE", "GeneF")
)
for(cluster in names(marker_genes)){
  integrated_data <- RenameIdents(integrated_data, which(startsWith(rownames(integrated_data@meta.data), marker_genes[[cluster]])) <- cluster)
}

# ---- 5. Analysis ----
# Differential Expression
# (Example for a specific comparison, can be adapted for others)
msc_response <- FindMarkers(integrated_data, ident.1 = "MSC_responder", ident.2 = "MSC_non_responder")

# Proportions
cell_prop <- prop.table(table(Idents(integrated_data), integrated_data$Condition))

# ---- 6. Saving Results ----
write.csv(msc_response, file = "msc.response.res_nonres.csv")
write.csv(cell_prop, file = "cell_proportions.csv")
saveRDS(integrated_data, file = "integrated_data.rds")
