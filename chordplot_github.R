setwd("../../Volumes/mc3m-special/scRNA-seq Analysis/Combined/NoErythrocytes/")
install.packages("circlize")
library(circlize)

#Disease
Immune_system_disease = c("MYOM2", "STAT3", "REL", "ITGAL", "BCL11B", "PRDM1", "ITGB2", "CXCR4", "ZFP36L1", "ITGAM", "TNFAIP3")
Primary_immunodeficiency_disease = c("STAT3", "REL", "ITGAL", "BCL11B", "PRDM1", "ITGB2", "CXCR4", "ZFP36L1", "ITGAM", "TNFAIP3")
Autoimmune_disease_of_musculoskeletal_system =  c("REL", "PRDM1", "ZFP36L1", "ITGAM", "TNFAIP3")
#Autoimmune_disease = c("STAT3", "REL", "PRDM1", "ZFP36L1", "ITGAM", "TNFAIP3")
Rheumatoid_arthritis = c("REL", "PRDM1", "ZFP36L1", "TNFAIP3", "TGFB1", "ITGAL", "JUN", "ITGB2", "CCL5")
Leukocyte_adhesion_deficiency = c("ITGAL", "ITGB2")
Connective_tissue_disease = c("TGFB1", "MYOM2", "REL", "PRDM1", "ZFP36L1", "ITGAM", "TNFAIP3")
Combined_immunodeficiency = c("ITGAL", "BCL11B", "ITGB2")


# KEGG
Osteoclast_differentiation = c("MAPK1", "NFKBIA", "TGFB1", "JUND", "JUNB", "JUN")
# Rheumatoid_arthritis = c("TGFB1", "ITGAL", "JUN", "ITGB2", "CCL5")
Regulation_of_actin_cytoskeleton = c("MAPK1", "MYH9", "ITGAL", "ITGB2", "CXCR4", "ITGAM")
IL17_signaling_pathway = c("MAPK1", "NFKBIA", "JUND", "JUN", "TNFAIP3")
Cellular_senescence = c("MAPK1", "TGFB1", "ETS1", "ZFP36L1")
NK_cell_mediated_cytotoxicity = c("MAPK1", "ITGAL", "ITGB2", "PTPN6") 
NODlike_receptor_signaling_pathway = c("MAPK1", "NFKBIA", "JUN", "CCL5", "TNFAIP3")


# RCTM
Signaling_by_Interleukins = c("MAPK1", "NFKBIA", "TGFB1", "STAT3", "JUNB", "JUN", "ITGB2", "PTPN6", "ITGAM", "CCL5")
TLR4_Cascade = c("MAPK1", "NFKBIA", "JUN", "ITGB2", "ITGAM")
Cell_surface_interactions_vascular_wall = c("TGFB1", "ITGAL", "ITGB2", "PTPN6", "ITGAM")
Growth_hormone_receptor_signaling = c("MAPK1", "STAT3", "PTPN6")
Signaling_by_NTRK1_TRKA = c("MAPK1", "JUND", "STAT3", "JUNB")











genes = c(Immune_system_disease, Primary_immunodeficiency_disease, Autoimmune_disease_of_musculoskeletal_system, Rheumatoid_arthritis, Leukocyte_adhesion_deficiency, 
          Connective_tissue_disease, Combined_immunodeficiency, Osteoclast_differentiation, Regulation_of_actin_cytoskeleton, 
          IL17_signaling_pathway, Cellular_senescence, NK_cell_mediated_cytotoxicity, NODlike_receptor_signaling_pathway, Signaling_by_Interleukins, 
          TLR4_Cascade, Cell_surface_interactions_vascular_wall, Growth_hormone_receptor_signaling, Signaling_by_NTRK1_TRKA)

# Pathways
path=c(rep("Immune_system_disease", 11), rep("Primary_immunodeficiency_disease", 10), 
       rep("Autoimmune_disease_of_musculoskeletal_system", 5), rep("Rheumatoid_arthritis", 9), rep("Leukocyte_adhesion_deficiency", 2),
       rep("Connective_tissue_disease", 7), rep("Combined_immunodeficiency", 3), rep("Osteoclast_differentiation", 6),
       rep("Regulation_of_actin_cytoskeleton", 6), rep("IL17_signaling_pathway", 5), rep("Cellular_senescence", 4), rep("NK_cell_mediated_cytotoxicity", 4),
       rep("NODlike_receptor_signaling_pathway", 5), rep("Signaling_by_Interleukins", 10), rep("TLR4_Cascade", 5), rep("Cell_surface_interactions_vascular_wall", 5), 
       rep("Growth_hormone_receptor_signaling", 3), rep("Signaling_by_NTRK1_TRKA", 4))


#path=c(rep("MHCII antigen", 5), rep("Phosphorylation CD3 and TCR", 4), rep("ZAP-70 to Immunological synapse", 4), rep("Interferon gamma signaling", 5), rep("TCR signaling", 4), rep("PD-1 signaling", 4), rep("Cytokine Signaling in Immune systemg", 11), rep("Adaptive Immune System", 7), rep("Costimulation by the CD28 family", 4), rep("CellCycle-Mitotic", 6), rep("G0 and EarlyG1", 2), rep("Mitotic G1phase and G1_S transition", 3), rep("IL6 family signaling", 2), rep("IL4 and IL13 signaling", 3))



# dataframe
dat = data.frame(genes, path)

# determine colors





grid.colors = c(MYOM2 = "skyblue3", STAT3 = "skyblue3", REL = "skyblue3", ITGAL = "skyblue3", BCL11B = "skyblue3", PRDM1 = "skyblue3", ITGB2 = "skyblue3", CXCR4 = "skyblue3", ZFP36L1 = "skyblue3", ITGAM = "skyblue3", 
                TNFAIP3 = "skyblue3", STAT3 = "skyblue3", REL = "skyblue3", ITGAL = "skyblue3", BCL11B = "skyblue3", PRDM1 = "skyblue3", ITGB2 = "skyblue3", CXCR4 = "skyblue3", ZFP36L1 = "skyblue3", ITGAM = "skyblue3", TNFAIP3 = "skyblue3",
                REL = "skyblue3", PRDM1 = "skyblue3", ZFP36L1 = "skyblue3", ITGAM = "skyblue3", TNFAIP3 = "skyblue3", REL = "skyblue3", PRDM1 = "skyblue3", ZFP36L1 =  "skyblue3", TNFAIP3 = "skyblue3", TGFB1 = "skyblue3", ITGAL = "skyblue3", JUN = "skyblue3", ITGB2 = "skyblue3", CCL5 = "skyblue3", ITGAL = "skyblue3", ITGB2 = "skyblue3",
                TGFB1 = "skyblue3", MYOM2 = "skyblue3", REL = "skyblue3", PRDM1 = "skyblue3", ZFP36L1 = "skyblue3", ITGAM = "skyblue3", TNFAIP3 = "skyblue3", ITGAL = "skyblue3", BCL11B = "skyblue3", ITGB2 = "skyblue3", 
                
                MAPK1 = "yellowgreen", NFKBIA= "yellowgreen", TGFB1 = "yellowgreen", JUND = "yellowgreen", JUNB = "yellowgreen", JUN = "yellowgreen", 
                MAPK1 = "yellowgreen", MYH9 = "yellowgreen", ITGAL = "yellowgreen", ITGB2 = "yellowgreen", CXCR4 = "yellowgreen", ITGAM = "yellowgreen", MAPK1 = "gray", NFKBIA = "yellowgreen", JUND = "yellowgreen", JUN = "yellowgreen", TNFAIP3 = "yellowgreen", MAPK1 = "yellowgreen", TGFB1 = "yellowgreen", ETS1 = "yellowgreen", ZFP36L1 = "yellowgreen",
                MAPK1 = "yellowgreen", ITGAL = "yellowgreen", ITGB2 = "yellowgreen", PTPN6 = "yellowgreen", MAPK1 = "yellowgreen", NFKBIA = "yellowgreen", JUN = "yellowgreen", CCL5 = "yellowgreen", TNFAIP3 = "yellowgreen", 
                
                
                MAPK1 = "violetred1", NFKBIA = "violetred1", TGFB1 = "violetred1", STAT3 = "violetred1", JUNB = "violetred1", JUN = "violetred1", ITGB2 = "violetred1", PTPN6 = "violetred1", ITGAM = "violetred1", CCL5 = "violetred1", 
                MAPK1 = "violetred1", NFKBIA = "violetred1", JUN = "violetred1", ITGB2 = "violetred1", ITGAM = "violetred1", TGFB1 = "violetred1", ITGAL = "violetred1", ITGB2 = "violetred1", PTPN6 = "violetred1", ITGAM = "violetred1", MAPK1 = "violetred1", STAT3 = "violetred1", PTPN6 = "violetred1", MAPK1 = "violetred1", JUND = "violetred1", STAT3 = "violetred1", JUNB = "violetred1", 
                
                Immune_system_disease = "skyblue3", Primary_immunodeficiency_disease = "skyblue3", 
                Autoimmune_disease_of_musculoskeletal_system = "skyblue3", 
                Rheumatoid_arthritis = "skyblue3", Leukocyte_adhesion_deficiency = "skyblue3", 
                Connective_tissue_disease = "skyblue3", Combined_immunodeficiency = "skyblue3", 
                
                Osteoclast_differentiation = "yellowgreen", Regulation_of_actin_cytoskeleton = "yellowgreen", 
                IL17_signaling_pathway = "yellowgreen", Cellular_senescence = "yellowgreen", NK_cell_mediated_cytotoxicity = "yellowgreen", 
                NODlike_receptor_signaling_pathway = "yellowgreen", 
                
                Signaling_by_Interleukins = "violetred1", TLR4_Cascade = "violetred1",
                Cell_surface_interactions_vascular_wall = "violetred1", Growth_hormone_receptor_signaling = "violetred1", Signaling_by_NTRK1_TRKA = "violetred1")










#circos.par(start.degree = 90, clock.wise = FALSE)
circos.par(start.degree = 135, track.margin = c(.01,0.01))
# Create chord diagram
#pdf("Chord diagram_nk.pdf")
#chordDiagram(as.data.frame(dat), transparency = 0.2, grid.col = grid.colors, 
#annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dat))))))

chordDiagram(as.data.frame(dat), grid.col = grid.colors, transparency = 0.5, 
             annotationTrack = "grid", annotationTrackHeight = 0.08,
             direction.type = c("diffHeight", "arrows"), preAllocateTracks = list(list(track.height = 0.03)))


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, 
              adj = c(0, 0.3), cex = 0.5)
}, bg.border = NA)

#circos.track(track.index = 1, panel.fun = function(x, y) {
#circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#facing = "clockwise", niceFacing = TRUE, 
#adj = c(0, 0.3), cex = 0.8)
#}, bg.border = NA) 

legend(-0.75, -0.55, legend = c("Disease", "KEGG", "Reactome"), bty = "n",
       fill = c("skyblue3", "yellowgreen", "violetred1"), cex = 0.5, adj = c(0.1, 0.5))
dev.off()
