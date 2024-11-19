# **Single-Cell Transcriptome Analysis of Osteoarthritis and Non-OA Bone Marrow**

This repository contains the scripts and data used for the analyses described in the paper:  
**"Single-cell transcriptome and crosstalk analysis reveals immune alterations and key pathways in the bone marrow of knee OA patients"**  
Published in iScience (2024).  
[Access the article here](https://doi.org/10.1016/j.isci.2024.110827)

## **Overview**
This study investigates the systemic immune alterations in bone marrow aspirate concentrates (BMAC) from knee osteoarthritis (OA) patients. Using single-cell RNA sequencing (scRNA-seq) and computational analyses, we highlight changes in immune and stromal cell populations, intercellular communication, and key signaling pathways involved in OA progression.

## **Repository Contents**
1. **Data Preprocessing Scripts**
   - Quality control and normalization of scRNA-seq data.
   - Integration and clustering using the Seurat pipeline.
2. **Differential Expression Analysis**
   - Identification of key genes differentiating OA and non-OA samples.
3. **Pathway and Cell-Cell Communication Analysis**
   - Pathway enrichment using GSEA and ligand-receptor interaction analysis via NicheNet and CellChat.
4. **Visualization Scripts**
   - UMAP plots for cell clustering and marker expression.
   - Stacked bar plots for cell type proportions.
5. **Supplementary Analysis**
   - Variance partitioning, pseudotime trajectory, and validation dataset comparisons.

## **Getting Started**
### **Prerequisites**
- Python (≥3.8)
- R (≥4.0) with the following packages:
  - Seurat
  - CellChat
  - NicheNet
  - GSEA
- Software for single-cell data processing:
  - 10X Genomics Cell Ranger

### **Setup**
Clone this repository and install the required dependencies:
```bash
git clone https://github.com/pchatterjee7/BMAC_OA_nonOA.git
cd BMAC_OA_nonOA

Data Access
Raw and processed scRNA-seq data are deposited in GEO (GSE274018) and SRA (PRJNA1144164). Access requests can be directed to the lead contact as per publication guidelines.

Acknowledgments
This study was funded by The Billie and Bernie Marcus Foundation and the Georgia Research Alliance.


