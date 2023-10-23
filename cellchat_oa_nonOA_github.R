library(CellChat)
library(Seurat)
library(patchwork)

Examples
## Not run: 
# Input is a data matrix
## create a dataframe consisting of the cell labels
# meta = data.frame(labels = cell.labels, row.names = names(cell.labels))
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# input is a Seurat object
## use the default cell identities of Seurat object
# cellchat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
## use other meta information as cell groups
# cellchat <- createCellChat(object = seurat.obj, group.by = "seurat.clusters")

# input is a SingleCellExperiment object
# cellchat <- createCellChat(object = sce.obj, group.by = "sce.clusters")


##################################bmac#################################

bmac <- readRDS("~/MILES/bmac_noerythrocytes_renamed.rds")
########################### NONOA ############################

bmac_nonOA <- readRDS("~/MILES/Subsets/nonOA_group.rds")

## if using as a seurat V3 object
data.input <- GetAssayData(bmac_nonOA, assay = "RNA", slot = "data")
labels <- Idents(bmac_nonOA)
meta <- data.frame(group = labels, row.names = names(labels))
View(meta)

## input is a Seurat object
## use the default cell identities of Seurat object
#cellchat <- createCellChat(object = bmac_nonOA, group.by = "ident", assay = "RNA")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
## If cell mata information is not added when creating CellChat object, 
#USERS can also add it later using addMeta, and set the default cell identities using setIdent.
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "group") # set "labels" as default cell identity
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))


CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", top = 0.1)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", top = 0.1)
saveRDS(cellchat, file = "cellchat.nonOA.rds")



mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, vertex.label.cex = 0.5, vertex.size = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], edge.label.cex = 0.1)
}
#########PATHWAY VISUALIZATION########################
cellchat@netP$pathways
#[1] "MHC-I"    "MIF"      "CLEC"     "CD99"     "MHC-II"   "GALECTIN" "APP"      "ANNEXIN" 
#[9] "LCK"      "CD22"     "CD45"     "RESISTIN" "THBS"     "ITGB2"    "MK"       "ADGRE5"  
#[17] "BAFF"     "SELL"     "SELPLG"   "PECAM1"   "CXCL"     "CD23"     "COLLAGEN" "ICAM"    
#[25] "FLT3"     "FN1"      "IL16"     "BAG"      "ALCAM"    "CD6"      "LAMININ"  "IL2"     
#[33] "SEMA4"    "CSF"      "THY1"     "VCAM"     "MPZ"      "ANGPTL"   "TENASCIN" "FGF" 
pathways.show <- cellchat@netP$pathways

path_sub <- function(ref, total = pathways.show){
  CellChatDB <- CellChatDB.human
  db_sub <- subsetDB(CellChatDB, search = ref)
  path_names <- unique(db_sub[["interaction"]][["pathway_name"]])
  path_intersect <- intersect(path_names, total)
}

sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)
par(mfrow=c(1,1))

netVisual_aggregate(cellchat, signaling = sec_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)
netVisual_aggregate(cellchat, signaling = ecm_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)
netVisual_aggregate(cellchat, signaling = cc_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)


pathways.show <- c("BAFF") 
vertex.receiver = seq(1,4)

netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.5)

#circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Chord diagram
group.cellType <- c(rep("TC", 1, 3, 4, 6), rep("BC", 2, 10, 13), rep("Mono", 5, 8, 11), rep("NK", 7), rep("HSPC", 9), rep("DC", 12, 14), rep("MEGA", 15), rep("MSC", 16))
# grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")



#> [[1]]
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


##BUBBLE PLOT##
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CD45","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)


netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = 8, legend.pos.x = 10, lab.cex = 0.3)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CD99","CXCL"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 3, legend.pos.y = 3, lab.cex = 0.5)
#> Note: The first link end is drawn out of sector 'CXCL '.

plotGeneExpression(cellchat, signaling = "CXCL")
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
plotGeneExpression(cellchat, signaling = "IL2", enriched.only = FALSE)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "IL2"))
gg1 + gg2



ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 10, font.size = 8, font.size.title = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 10, font.size = 8, font.size.title = 10)
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "MIF"))
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, font.size = 8, height = 9)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

###incoming pattern
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

netAnalysis_dot(cellchat, pattern = "incoming")
saveRDS(cellchat, file = "cellchat_nonOA.rds")




########################### OA ############################
bmac_OA <- readRDS("~/10xGenomics/MILES/OA_group.rds")
data.input <- GetAssayData(bmac_OA, assay = "RNA", slot = "data")
labels <- Idents(bmac_OA)
meta <- data.frame(group = labels, row.names = names(labels))
#View(meta)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
#cellchat <- createCellChat(object = bmac_nonOA, group.by = "ident", assay = "RNA")

#cellchat <- addMeta(cellchat, meta = meta)
#cellchat <- setIdent(cellchat, ident.use = "group") # set "labels" as default cell identity
#levels(cellchat@idents)
#groupSize <- as.numeric(table(cellchat@idents))


CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", top = 0.1)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", top = 0.1)
saveRDS(cellchat, file = "cellchat.OA.rds")



mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, vertex.label.cex = 0.5, vertex.size = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], edge.label.cex = 0.1)
}
#########PATHWAY VISUALIZATION########################
cellchat@netP$pathways
#[1] "MHC-I"    "MIF"      "CLEC"     "CD99"     "MHC-II"   "GALECTIN" "APP"      "ANNEXIN" 
#[9] "LCK"      "CD22"     "CD45"     "RESISTIN" "THBS"     "ITGB2"    "MK"       "ADGRE5"  
#[17] "BAFF"     "SELL"     "SELPLG"   "PECAM1"   "CXCL"     "CD23"     "COLLAGEN" "ICAM"    
#[25] "FLT3"     "FN1"      "IL16"     "BAG"      "ALCAM"    "CD6"      "LAMININ"  "IL2"     
#[33] "SEMA4"    "CSF"      "THY1"     "VCAM"     "MPZ"      "ANGPTL"   "TENASCIN" "FGF" 

path_sub <- function(ref, total = pathways.show){
  CellChatDB <- CellChatDB.human
  db_sub <- subsetDB(CellChatDB, search = ref)
  path_names <- unique(db_sub[["interaction"]][["pathway_name"]])
  path_intersect <- intersect(path_names, total)
}

sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)
par(mfrow=c(1,1))

netVisual_aggregate(cellchat, signaling = sec_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)
netVisual_aggregate(cellchat, signaling = ecm_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)
netVisual_aggregate(cellchat, signaling = cc_paths, top = 0.1, layout = "circle", pt.title = 12, title.space = 6)



pathways.show <- c("CXCL") 
vertex.receiver = seq(1,4)

netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.5)

#circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("TC", 1, 3, 4, 6), rep("BC", 2, 10, 13), rep("Mono", 5, 8, 11), rep("NK", 7), rep("HSPC", 9), rep("DC", 12, 14), rep("MEGA", 15), rep("MSC", 16)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,8) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")



#> [[1]]
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,5)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


##BUBBLE PLOT##
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CD45","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)


netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = 8, legend.pos.x = 10, lab.cex = 0.3)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CD99","CXCL"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 3, legend.pos.y = 3, lab.cex = 0.5)
#> Note: The first link end is drawn out of sector 'CXCL '.

plotGeneExpression(cellchat, signaling = "CXCL")
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
plotGeneExpression(cellchat, signaling = "IL2", enriched.only = FALSE)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "IL2"))
gg1 + gg2



ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 10, font.size = 8, font.size.title = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 10, font.size = 8, font.size.title = 10)
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "MIF"))
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, font.size = 8, height = 9)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

###incoming pattern
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

netAnalysis_dot(cellchat, pattern = "incoming")
saveRDS(cellchat, file = "cellchat_OA.rds")







##Merge and compare OA and nonOA##


#First dataset = nonOA, Second dataset = OA
cellchat_nonOA <- readRDS("~/10xGenomics/MILES/cellchat_miles/nonOA_cellchat/cellchat_nonOA.rds")
cellchat_OA <- readRDS("~/10xGenomics/MILES/cellchat_miles/OA_cellchat/cellchat_OA.rds")
object.list <- list(nonOA = cellchat_nonOA, OA = cellchat_OA)
#View(object.list)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

####Compare the number of interactions and interaction strength among different cell populations###
##The differential number of interactions or interaction 
#strength in the cell-cell communication network between two datasets 
#can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) 
#signaling in the second dataset compared to the first one##
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

####We can also show differential number of interactions or interaction strength in a 
##greater details using a heatmap. The top colored bar plot represents 
##the sum of column of values displayed in the heatmap (incoming signaling). 
##The right colored bar plot represents the sum of row of values (outgoing signaling). 
##In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset 
##compared to the first one.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


####The differential network analysis only works for pairwise datasets. If there are more datasets for comparison, we can directly show the number of interactions or interaction strength between any two cell populations in each dataset.

###To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#group.cellType <- c(rep("TC", 1, 3, 4, 6), rep("BC", 2, 10, 13), rep("Mono", 5, 8, 11), rep("NK", 7), rep("HSPC", 9), rep("DC", 12, 14), rep("MEGA", 15), rep("MSC", 16))
#group.cellType <- factor(group.cellType, levels = c("TC", "BC", "Mono", "NK", "HSPC", "DC", "MEGA", "MSC"))
group.cellType <- c('NK', 'Naive T cells', 'CD8 T cells', 'CD4 T cells', 'DNTLGALS3_high', 'MSC', 'HSPC', 'Mature B', 'cDC', 'pDC')
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


##compare celltypes#### Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MSC", font.size = 4)
#Visualizing differential outgoing and incoming signaling changes from nonOA to OA
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK", font.size = 4)

#Visualizing differential outgoing and incoming signaling changes from nonOA to OA
patchwork::wrap_plots(plots = list(gg1,gg2))


#CellChat performs joint manifold learning and classification of the inferred communication networks based on their functional and topological similarity. NB: Such analysis is applicable to more than two datasets.

#Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.

#Structural similarity: A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers. NB: Structural similarity analysis is applicable to multiple datsets with the same cell type composition or the vastly different cell type composition.

#Here we can run the manifold and classification learning analysis based on the functional similarity because the two datasets have the the same cell type composition.

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#Compute signaling network similarity for datasets 1 2 
cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")

#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space




# Hierarchy plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Circle plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

##### Compare the major sources and targets in 2D space
###Comparing the outgoing and incoming interaction strength in 2D space 
###allows ready identification of the cell populations with significant changes 
##in sending or receiving signals between different datasets.


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

##Furthermore, we can identify the specific signaling changes between OA and nonOA. 
## Identify signaling changes associated with one cell group

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mono CD14", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))                                                                                                                                                                   


### Part II: Identify the conserved and context-specific signaling pathways

### CellChat then can identify signaling networks with larger (or less) difference, signaling groups, and the conserved and context-specific 
###signaling pathways based on their cell-cell communication networks among multiple biological conditions.
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


library(ComplexHeatmap)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 8, font.size = 4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 8, font.size = 4)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 8, color.heatmap = "GnBu", font.size = 4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 8, color.heatmap = "GnBu", font.size = 4)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 8, color.heatmap = "OrRd", font.size = 4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 8, color.heatmap = "OrRd", font.size = 4)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


############Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "OA"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net.up <- subsetCommunication(cellchat, net = net, datasets = "OA",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "nonOA",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2
