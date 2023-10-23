

# Required Libraries
# ==================

# Install CRAN packages
install_required_packages <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
}

# Install Bioconductor packages
install_bioconductor_packages <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) BiocManager::install(new_pkgs, dependencies = TRUE)
}

# CRAN packages
cran_packages <- c("msigdbr", "dplyr", "purrr", "stringr", "magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager", "data.table", "doParallel",
                   "doRNG")
install_required_packages(cran_packages)

# Bioconductor packages
bioconductor_packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree", "Nebulosa")
install_bioconductor_packages(bioconductor_packages)

# Install GitHub packages
if (!requireNamespace("UCell", quietly = TRUE)) devtools::install_github("carmonalab/UCell")
if (!requireNamespace("irGSEA", quietly = TRUE)) devtools::install_github("chuiqin/irGSEA")

# Load Libraries
# ==============

library(Seurat)
library(SeuratData)
library(UCell)
library(irGSEA)

# Analysis
# ========

# ... (Your analysis code)

# For example, to read the seurat.obj data, update Seurat object and so on.
# data("seurat.obj")
# DefaultAssay(seurat.obj) <- "RNA"
# seurat.obj <- UpdateSeuratObject(seurat.obj)

seurat.obj <- irGSEA.score(object = seurat.obj, assay = "RNA", 
                           slot = "data", seeds = 123, ncores = 3,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "C5",  
                           subcategory = "GO:MF", geneid = "symbol",
                           method = c("AUCell", "UCell", "singscore", 
                                      "ssgsea"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = seurat.obj, 
                               group.by = "Condition",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot

#Show co-upregulated or co-downregulated gene sets per cluster in RRA.

#If error (argument “caller_env” is missing, with no default) occurs : please uninstall ggtree and run “remotes::install_github(”YuLab-SMU/ggtree”)“.

irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)
irGSEA.bubble.plot


irGSEA.upset.plot <- irGSEA.upset(object = result.dge, 
                                  method = "RRA")









