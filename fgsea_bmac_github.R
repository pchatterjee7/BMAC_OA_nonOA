####FGSEA

OA <- subset(bmac, subset = Condition == "OA")
nonOA <- subset(bmac, subset = Condition == "non OA")
oa.genes <- wilcoxauc(OA, 'ident')
head(oa.genes)
dplyr::count(oa.genes, group)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
library(presto)

msigdbr_show_species()
m_df<- msigdbr(species = "Homo sapiens", category = "H") ### choose appropriate database from msigdb
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets$GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP
oa.genes %>%
  dplyr::filter(group == "MSC") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)


msc.genes<- oa.genes %>%
  dplyr::filter(group == "MSC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks<- deframe(msc.genes)

head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()




ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()




#########nonOA

nonoa.genes <- wilcoxauc(nonOA, 'ident')
head(nonoa.genes)
dplyr::count(nonoa.genes, group)

msigdbr_show_species()
m_df<- msigdbr(species = "Homo sapiens", category = "H") ### choose appropriate database from msigdb
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets$GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP
nonoa.genes %>%
  dplyr::filter(group == "MSC") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)


msc.genes<- nonoa.genes %>%
  dplyr::filter(group == "MSC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks<- deframe(msc.genes)

head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()




ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()





############### FSGSEA OA #############

#### MSC OA ########
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tibble)
library(presto)

bmac <- readRDS("~/MILES/bmac_noerythrocytes_renamed.rds")
bmac$celltype <- Idents(bmac)


msc <- subset(bmac, idents = "MSC")
msc.genes <- wilcoxauc(msc, group_by = "Condition")
head(msc.genes)
dplyr::count(msc.genes, group)

msigdbr_species()

m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
msc.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_msc.genes<- msc.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_msc.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for MSC") + 
  theme_minimal()

msc.gsea <- as.data.frame(fgseaResTidy)

fwrite(msc.gsea, file ="msc_gsea.csv")

##### DNT ##########

dnt <- subset(bmac, idents = "DNTLGALS3_high")
dnt.genes <- wilcoxauc(dnt, group_by = "Condition")
head(dnt.genes)
dplyr::count(dnt.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
dnt.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_dnt.genes<- dnt.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_dnt.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for DNT") + 
  theme_minimal()

dnt.gsea <- as.data.frame(fgseaResTidy)

fwrite(dnt.gsea, file ="dnt_gsea.csv")


############## NK ####################

nk <- subset(bmac, idents = "NK")
nk.genes <- wilcoxauc(nk, group_by = "Condition")
head(nk.genes)
dplyr::count(nk.genes, group)

msigdbr_species()

m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
nk.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_nk.genes<- nk.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_nk.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p3 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for NK") + 
  theme_minimal()

nk.gsea <- as.data.frame(fgseaResTidy)

fwrite(nk.gsea, file ="nk_gsea.csv")





####### MEGA ############

mega <- subset(bmac, idents = "MEGA")
mega.genes <- wilcoxauc(mega, group_by = "Condition")
head(mega.genes)
dplyr::count(mega.genes, group)

msigdbr_species()

m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
mega.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_mega.genes<- mega.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_mega.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p4 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for MEGA") + 
  theme_minimal()

mega.gsea <- as.data.frame(fgseaResTidy)

fwrite(mega.gsea, file ="mega_gsea.csv")



############## pDC #################


pDC <- subset(bmac, idents = "pDC")
pDC.genes <- wilcoxauc(pDC, group_by = "Condition")
head(pDC.genes)
dplyr::count(pDC.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pDC.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_pDC.genes<- pDC.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_pDC.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p5 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for pDC") + 
  theme_minimal()


pdc.gsea <- as.data.frame(fgseaResTidy)

fwrite(pdc.gsea, file ="pdc_gsea.csv")

################# MONO CD16 #################

Mono_CD16 <- subset(bmac, idents = "Mono CD16")
Mono_CD16.genes <- wilcoxauc(Mono_CD16, group_by = "Condition")
head(Mono_CD16.genes)
dplyr::count(Mono_CD16.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
Mono_CD16.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_Mono_CD16.genes<- Mono_CD16.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_Mono_CD16.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p6 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for Mono_CD16") + 
  theme_minimal()



mono16.gsea <- as.data.frame(fgseaResTidy)

fwrite(mono16.gsea, file ="mono16_gsea.csv")


################# MONO CD14 #################

Mono_CD14 <- subset(bmac, idents = "Mono CD14")
Mono_CD14.genes <- wilcoxauc(Mono_CD14, group_by = "Condition")
head(Mono_CD14.genes)
dplyr::count(Mono_CD14.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
Mono_CD14.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_Mono_CD14.genes<- Mono_CD14.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_Mono_CD14.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p7 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for Mono_CD14") + 
  theme_minimal()


mono14.gsea <- as.data.frame(fgseaResTidy)

fwrite(mono14.gsea, file ="mono14_gsea.csv")



############### PLASMA CELLS ##############

plasmab <- subset(bmac, idents = "Plasma B")
plasmab.genes <- wilcoxauc(plasmab, group_by = "Condition")
head(plasmab.genes)
dplyr::count(plasmab.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
plasmab.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_plasmab.genes<- plasmab.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_plasmab.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p8 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for Plasma B") + 
  theme_minimal()

plasmab.gsea <- as.data.frame(fgseaResTidy)

fwrite(plasmab.gsea, file ="plasmab_gsea.csv")




############## cDC ##############


cDC <- subset(bmac, idents = "cDC")
cDC.genes <- wilcoxauc(cDC, group_by = "Condition")
head(cDC.genes)
dplyr::count(cDC.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
cDC.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_cDC.genes<- cDC.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_cDC.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p9 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for cDC") + 
  theme_minimal()

cdc.gsea <- as.data.frame(fgseaResTidy)

fwrite(cdc.gsea, file ="cdc_gsea.csv")





############### MATURE B ###################

matureb <- subset(bmac, idents = "Mature B")
matureb.genes <- wilcoxauc(matureb, group_by = "Condition")
head(matureb.genes)
dplyr::count(matureb.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
matureb.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_matureb.genes<- matureb.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_matureb.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p10 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for Mature B") + 
  theme_minimal()

matureb.gsea <- as.data.frame(fgseaResTidy)

fwrite(matureb.gsea, file ="matureb_gsea.csv")




p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10





############ CD4 T cells #################


cd4t <- subset(bmac, idents = "CD4 T cells")
cd4t.genes <- wilcoxauc(cd4t, group_by = "Condition")
head(cd4t.genes)
dplyr::count(cd4t.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
cd4t.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_cd4t.genes<- cd4t.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_cd4t.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p11 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for CD4T") + 
  theme_minimal()


cd4t.gsea <- as.data.frame(fgseaResTidy)

fwrite(cd4t.gsea, file ="cd4t_gsea.csv")



############ CD8 T cells #################


cd8t <- subset(bmac, idents = "CD8 T cells")
cd8t.genes <- wilcoxauc(cd8t, group_by = "Condition")
head(cd8t.genes)
dplyr::count(cd8t.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
cd8t.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_cd8t.genes<- cd8t.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_cd8t.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p12 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for CD8T") + 
  theme_minimal()


cd8t.gsea <- as.data.frame(fgseaResTidy)

fwrite(cd8t.gsea, file ="cd8t_gsea.csv")



############ NAIVE T cells #################


naive <- subset(bmac, idents = "Naive T cells")
naive.genes <- wilcoxauc(naive, group_by = "Condition")
head(naive.genes)
dplyr::count(naive.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
naive.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_naive.genes<- naive.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_naive.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p13 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for NAIVE-T") + 
  theme_minimal()


naivet.gsea <- as.data.frame(fgseaResTidy)

fwrite(naivet.gsea, file ="naivet_gsea.csv")




############ PROGENITOR B cells #################


prob <- subset(bmac, idents = "Progenitor B")
prob.genes <- wilcoxauc(prob, group_by = "Condition")
head(prob.genes)
dplyr::count(prob.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
prob.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_prob.genes<- prob.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_prob.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p14 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for PROGENITOR B") + 
  theme_minimal()


progenitorb.gsea <- as.data.frame(fgseaResTidy)

fwrite(progenitorb.gsea, file ="progenitorb_gsea.csv")




############ Mono/Mono progenitor #################


mono_pro <- subset(bmac, idents = "Mono/Mono progenitor")
mono_pro.genes <- wilcoxauc(mono_pro, group_by = "Condition")
head(mono_pro.genes)
dplyr::count(mono_pro.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
mono_pro.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_mono_pro.genes<- mono_pro.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_mono_pro.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p15 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for MONO PROGENITOR") + 
  theme_minimal()

monopro.gsea <- as.data.frame(fgseaResTidy)

fwrite(monopro.gsea, file ="mono_progeni_gsea.csv")




############ HSPC #################


hspc <- subset(bmac, idents = "HSPC")
hspc.genes <- wilcoxauc(hspc, group_by = "Condition")
head(hspc.genes)
dplyr::count(hspc.genes, group)

msigdbr_species()

#m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
hspc.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)



oa_hspc.genes<- hspc.genes %>%
  dplyr::filter(group == "OA") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(oa_hspc.genes)

head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p16 <- ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA for HSPC") + 
  theme_minimal()

hspc.gsea <- as.data.frame(fgseaResTidy)

fwrite(hspc.gsea, file ="hspc_gsea.csv")




p10 + p11 + p12 + p13 + p15 + p16


grid.arrange(p10,p11,p12,p13,p15,p16,nrow = 4, widths = c(1,1.5))



#####GSEA HEATMAP COMPLEXHEATMAP
library(dplyr)
setwd("~/MILES/fgsea/final_fig")
nci <- read.csv('top10_transformed.csv')
rownames(nci) <- nci$CellType
nci <- select(nci, -c('CellType'))
#col_fun = circlize::colorRamp2(c(-8, 0, 8), c("blue", "white", "red"))
col_fun = circlize::colorRamp2(c(-8, -4, 0, 4, 8), c("blue", "mediumpurple2", "white", "lightsalmon","red"))
lgd = Legend(col_fun = col_fun, title = "GSEA_NES")
Heatmap(nci, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))


my_mat <- as.matrix(nci) 

Heatmap(my_mat, name = "GSEA_NES", cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2), column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10))


library(ComplexHeatmap)

nci <- read.csv('top10_transformed.csv')
rownames(nci) <- nci$CellType
nci <- nci[, -1] # Remove the CellType column

col_fun = circlize::colorRamp2(c(-8, -4, 0, 4, 8), c("blue", "mediumpurple2", "white", "lightsalmon","red"))
lgd = Legend(col_fun = col_fun, title = "GSEA_NES")

my_mat <- as.matrix(nci)

Heatmap(my_mat, name = "GSEA_NES", cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2), column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10))






library(ComplexHeatmap)

# Create a matrix of gene expression values (log2 transformed)
expr_mat <- matrix(rnorm(1000), ncol = 20)
rownames(expr_mat) <- paste0("Gene", 1:100)
colnames(expr_mat) <- c(paste0("Male_Sample", 1:5), paste0("Female_Sample", 1:5), paste0("Normal_Sample", 1:5), paste0("Tumor_Sample", 1:5))
expr_mat <- log2(expr_mat)

# Define colors for the heatmap
my_col <- colorRamp2(c(-2, 0, 2), c("#386cb0", "#f0f0f0", "#d95f02"))

# Create factor variables for the groups (male and female) and (normal and tumor)
groups1 <- factor(c(rep("Male", 5), rep("Female", 5)))
groups2 <- factor(c(rep("Normal", 5), rep("Tumor", 5)))

# Create the heatmap
Heatmap(expr_mat, col = my_col, show_row_names = TRUE, name = "Log2 Expression", column_title = "Groups", column_names_gp = gpar(fontsize = 12), row_title = "Genes", row_names_gp = gpar(fontsize = 8), column_km = 2, km_alpha = 0.2, group_column = c(groups1, groups2), group_names = c("Male", "Female", "Normal", "Tumor"), heatmap_legend_param = list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 10)))

