setwd("./results/Processing")
#####Seurat Object generation and QC metrics
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(limma)
library(clustree)
library(scMayoMap)
library(scDblFinder)
library(SingleR)
library(argparser)
library(biomaRt)
library(scuttle)
library(tidyverse)

prefix <- "9823"
################################################################################################
#Integration after filtering. We used the non filtered for one set of analyses and the filtered for another. Both processes are shown though they are functionally the same after subsetting. 
#####INTEGRATION FILTER APPLIED
#SCT Transform regressing out library size, percentage intergenic, and mitochondrial genes.
combined <- readRDS(paste0("./",prefix,"_pre_QC.rds"))
combinedMAD <- subset(combined,  subset = QC_MAD == 'Pass')
combinedMAD <- SCTransform(combinedMAD, vars.to.regress = c("percent.mt", "pct_intergenic"), variable.features.n = 3000)
combinedMAD <- RunPCA(combinedMAD, npcs = 100, verbose = FALSE)	#determine linear dimensionality
combinedMAD <- RunUMAP(combinedMAD, dims = 1:40) #Projection map
DimPlot(combinedMAD, reduction = "umap", group.by = "orig.ident") + theme(text = element_text(size = 28))
#Integrate layers
combinedMAD <- IntegrateLayers(object = combinedMAD, method = HarmonyIntegration, normalization.method = "SCT", verbose = F, .options = harmony_options(max_iter = 50))
combinedMAD <- RunUMAP(combinedMAD, dims = 1:40, reduction = "harmony") #Projection map
png('OrigIdent_UMAP.png', width = 800, height = 700)
DimPlot(combinedMAD, reduction = "umap", group.by = "orig.ident", label = F) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
ids <- c("S23_105-2-Thalamus_THA", "S25_105-3-Thalamus_THA", "S27_105-10-Thalamus_THA", "S22_105-2-Striatum_STR", "S24_105-3-Striatum_STR", "S26_105-10-Striatum_STR") #Individuals since we don't want differences between individuals
tissue_vec <- c(rep("Thalamus", 3), rep("Striatum", 3)) #Tissue assignment sicne we don't want difference due to tissue difference
names(tissue_vec) <- ids #Name the tissues by sample ID
cell_tissue <- tissue_vec[combinedMAD$orig.ident] # tissue per cell
names(cell_tissue) <- colnames(combinedMAD) # names must be cell IDs
combinedMAD <- AddMetaData(combinedMAD, metadata = cell_tissue, col.name = "tissue") # Assign to Seurat
png('Tissue_UMAP.png', width = 800, height = 700)
DimPlot(combinedMAD, reduction = "umap", group.by = "tissue", label = F) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
saveRDS(combinedMAD, paste0(prefix,"_POST_INT_MAD.rds")) #Save it
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


################################################################################################
#####CLUSTERING FILTER APPLIED
combinedMAD <- FindNeighbors(combinedMAD, reduction = "harmony", dims = 1:40)	#neighbor detection function
combinedMAD <- FindClusters(object = combinedMAD, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) #running at mutiple resolutions for subsequent rerunning
clus_res <- "0.8"
ident_select <- paste0("SCT_snn_res.", clus_res)	#select an indent to work with. 
Idents(object = combinedMAD) <- ident_select
png('Cluster_UMAP.png', width = 800, height = 700)
DimPlot(combinedMAD, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + NoLegend()
dev.off()
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


################################################################################################
#####MARKER DETECTION FILTER APPLIED
combinedMAD <- PrepSCTFindMarkers(combinedMAD)
ALL_markers <- FindAllMarkers(combinedMAD, only.pos = TRUE)
ALL_markers %>%
	group_by(cluster) %>%
	top_n(n = 20, wt = avg_log2FC) -> top20
write.csv(top20, paste0(prefix, "_Top_20_Feats_MAD.csv"))	#output top 20 features
write.csv(ALL_markers, paste0(prefix, "_All_Feats_MAD.csv"))	#output all features
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
