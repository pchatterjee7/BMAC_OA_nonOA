####### progeny pathway scoring #######

bmac <- readRDS("~/bmac_final_with_cellid.rds")
bmac_20 <- readRDS("~/bmac_noerythrocytes_renamed.rds")

levels(bmac_20)
cell_types <- bmac@meta.data$cell_type
new_order <- c("DNT", "Mature B", "CD8 T cells", "Mono CD14", "Naive CD4+T", 
               "NK", "Mono Progeni", "HSPC", "Progenitor B", "Mono CD16", 
               "cDC", "Plasma B", "pDC", "MEGA", "MSC")
cell_types <- factor(cell_types, levels = new_order)
bmac@meta.data$cell_type <- cell_types

levels(bmac)

bmac$cell_type
Idents(bmac) <- "cell_type"

Idents(bmac)

CellsClusters <- data.frame(Cell = names(Idents(bmac)), 
                            CellType = as.character(Idents(bmac)),
                            stringsAsFactors = FALSE)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
bmac <- progeny(bmac, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
bmac <- Seurat::ScaleData(bmac, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(bmac, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

#We plot the different pathway activities for the different cell populations

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)
