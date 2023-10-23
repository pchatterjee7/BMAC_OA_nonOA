library(CellChat)
library(Seurat)
msc_ct <- readRDS("~/miles_outcome/msc/msc_temp_cellid.rds")
bmac <- readRDS("~/MILES/bmac_combined_noerythrocytes.rds")
data.inputOA <- GetAssayData(bmac_OA, assay = "RNA", slot = "data")
labels <- Idents(bmac_OA)
meta <- data.frame(group = labels, row.names = names(labels))
#View(meta)
cellchat.OA <- createCellChat(object = data.inputOA, meta = meta, group.by = "group")

data.inputnonOA <- GetAssayData(bmac_nonOA, assay = "RNA", slot = "data")
labels <- Idents(bmac_nonOA)
meta <- data.frame(group = labels, row.names = names(labels))
#View(meta)
cellchat.nonOA <- createCellChat(object = data.inputnonOA, meta = meta, group.by = "group")


cellchat.OA <- addMeta(cellchat.OA, meta = meta, meta.name = "labels")
cellchat.OA <- setIdent(cellchat.OA, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.OA@idents)
groupSize <- as.numeric(table(cellchat.OA@idents))

cellchat.nonOA <- addMeta(cellchat.nonOA, meta = meta, meta.name = "labels")
cellchat.nonOA <- setIdent(cellchat.nonOA, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.nonOA@idents)
groupSize <- as.numeric(table(cellchat.nonOA@idents))


library(patchwork)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat.OA@DB <- CellChatDB.use
cellchat.nonOA@DB <- CellChatDB.use

cellchat.OA <- subsetData(cellchat.OA)
cellchat.nonOA <- subsetData(cellchat.nonOA)

cellchat.OA <- identifyOverExpressedGenes(cellchat.OA)
cellchat.nonOA <- identifyOverExpressedInteractions(cellchat.nonOA)


###optional
cellchat.OA <- projectData(cellchat.OA, PPI.human)
cellchat.nonOA <- projectData(cellchat.nonOA, PPI.human)

computeAveExpr(cellchat.OA, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1)

cellchat.OA <- computeCommunProb(cellchat.OA, population.size = TRUE)
cellchat.nonOA <- computeCommunProb(cellchat.nonOA, population.size = TRUE)

cellchat.OA <- filterCommunication(cellchat.OA, min.cells = 10)
cellchat.nonOA <- filterCommunication(cellchat.nonOA, min.cells = 10)















options(stringsAsFactors = FALSE)
load("ligand_receptor.RData")
load("lung_dat.Robj")
lung<-NormalizeData(lung_dat)
#Preprocess and create a Cell chat object
data.input<-GetAssayData(lung,assay = "RNA",slot = "data")
meta<-lung_dat@meta.data
cellchat <- createCellChat(object = data.input,meta=meta,group.by = "cell.type.ident")
#Define a database
CellChatDB1.use<-data.frame(ligand_receptor[[1]])
cellchat@DB<-CellChatDB1.use
cellchat <- subsetData(cellchat)





######## re-order celltype names ##########

#The function you posted works and reorders cell groups for some visualizations - it works well for chords, circles, and dotplots but does not reorder them for heatmaps. Here is the code I used:

#Set desired cell order
cell.levels <- c("your.cell.order1", "your.cell.order2")

#stash cell labels in scRNA object metadata
object$cellgroup <- Idents(object)

#Add cell labels to cellchat metadata
identity = data.frame(group = object$cellgroup, row.names = names(object$cellgroup))
unique(identity$group)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell_type")
levels(cellchat@idents) #check idents are correct

#Reorder cells
cellchat <- setIdent(cellchat, ident.use = "cell_type", levels = cell.levels)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

#I was then able to do the non-heatmap visualizations functions (all variants of netVisual_aggregate, etc) with the new order.

#Have you re-run computeCommunProbPathway , cellchat <- aggregateNet(cellchat), cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") .

#All the calculations after computeCommunProb should be re-run!


#a somewhat related question, I have two seurat object with the same cell types, but when I ran
#cellchat <- setIdent(cellchat, ident.use = "labels")
#two cellchat objects have different orders of the same cell types, the subsequent circle plot shows different color for the same cell type.
#Do I need to run this modified setIdent() to avoid this situation?
#just do cellchat <- setIdent(cellchat) on both objects or just one?

#You can simply address this issue by defining the factor levels in your labels.
object1@meta$labels = factor(object1@meta$labels, levels = labels.levels) 
object2@meta$labels = factor(object2@meta$labels, levels = labels.levels)





############# Let's re-order the 71 bmac dataset as the 21 sample ##################
levels(cellchat_71@idents)
# [1] "Naive CD4+T"  "NK"           "Mono CD14"    "CD8 T cells"  "Mono Progeni"
# [6] "Progenitor B" "Mature B"     "DNT"          "Mono CD16"    "cDC"         
# [11] "Plasma B"     "pDC"          "HSPC"         "MEGA"         "MSC"         
cellchat_20 <- readRDS("~/miles_outcome/bmac/group71/cellchat_71/cellchat_20_bmac.rds")
levels(cellchat_20@idents)
#[1] "DNTLGALS3_high"       "Mature B"             "Naive T cells"       
# [4] "CD8 T cells"          "Mono CD14"            "CD4 T cells"         
# [7] "NK"                   "Mono/Mono progenitor" "HSPC"                
# [10] "Progenitor B"         "Mono CD16"            "cDC"                 
# [13] "Plasma B"             "pDC"                  "MEGA"                
# [16] "MSC"                 

cell.levels_71 <- c("DNT", "Mature B", "CD8 T cells", "Mono CD14", "Naive CD4+T", 
                    "NK", "Mono Progeni", "HSPC", "Progenitor B", "Mono CD16", 
                    "cDC", "Plasma B", "pDC", "MEGA", "MSC")
bmac_71 <- readRDS("~/miles_outcome/bmac/group71/bmac_final_with_cellid.rds")
bmac_71$cellgroup <- Idents(bmac_71)
identity71 = data.frame(group = bmac_71$cellgroup, row.names = names(bmac_71$cellgroup))
unique(identity71$group)
cellchat71 <- addMeta(cellchat_71, meta = identity71, meta.name = "cell_type")
levels(cellchat71@idents) #check idents are correct
#Reorder cells
cellchat71 <- setIdent(cellchat71, ident.use = "cell_type", levels = cell.levels_71)
levels(cellchat71@idents)
groupSize <- as.numeric(table(cellchat71@idents))
library(patchwork)

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat71@DB <- CellChatDB.use
cellchat71 <- subsetData(cellchat71)
cellchat71 <- identifyOverExpressedGenes(cellchat71)
cellchat71 <- identifyOverExpressedInteractions(cellchat71)
cellchat71 <- computeCommunProb(cellchat71)
cellchat71 <- filterCommunication(cellchat71, min.cells = 10)
cellchat71 <- computeCommunProbPathway(cellchat71)
cellchat71 <- aggregateNet(cellchat71)
groupSize <- as.numeric(table(cellchat71@idents))
par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(cellchat20@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions 21 OA BMAC", top = 0.1)
p2 <- netVisual_circle(cellchat71@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions 71 OA BMAC", top = 0.1)

cellchat_nonOA <- readRDS("~/cellchat_miles/nonOA_cellchat/cellchat_nonOA.rds")
updateCellChat(cellchat_nonOA)
# An object of class CellChat created from a single dataset 
# 28980 genes.
# 66931 cells. 
# CellChat analysis of single cell RNA-seq data! 
p3 <- netVisual_bubble(cellchat_nonOA, sources.use = c(1:15), targets.use = 16, remove.isolate = FALSE, title.name = "Significant interaction from cell groups to MSCs in 21 healthy BMAC")
#Comparing communications on a single object 

p2 <- netVisual_bubble(cellchat, sources.use = c(1:14), targets.use = 15, remove.isolate = FALSE, title.name = "Significant interaction from cell groups to MSCs in 71 OA BMAC")
#Comparing communications on a single object 

p1 <- netVisual_bubble(cellchat_20, sources.use = c(1:15), targets.use = 16, remove.isolate = FALSE, title.name = "Significant interaction from cell groups to MSCs in 21 OA BMAC")
