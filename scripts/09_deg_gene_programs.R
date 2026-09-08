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
setwd("./results/comparison")
################################################################################################
#CELL TYPE EXPRESSION Integration of each level of filtering to ensure cell groupings are still meaningful as cell partioning.
#NOTE Integration should be run on each level of filtering. This would be combined, combined3. These are saved and the right point to run the steps below. Update as needed to run each through these steps.
combined <- readRDS("./ALL_INT_STR.rds") #Read multispeices object in
combined <- SCTransform(combined, ncells = 15000) #perform SCT normalization and transformation
combined <- RunPCA(combined, assay = "SCT", npcs=100) #Run dimensionality
combined <- IntegrateLayers(object = combined, method = HarmonyIntegration, normalization.method = "SCT", assay = "SCT", verbose = F, orig.reduction = "pca", new.reduction = 'harmony') 
combined <- RunUMAP(combined, assay = "SCT", dims = 1:50, reduction = "harmony") #Projection visulaization post integration
#Diagnostic plots post integration
DimPlot(combined, reduction = "umap", group.by = "orig.ident", raster=FALSE, alpha = 0.1) + theme(text = element_text(size = 28)) 
DimPlot(combined, reduction = "umap", group.by = "CT_toplevel", raster=FALSE) + theme(text = element_text(size = 28))  + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
#Save the objects after running.
saveRDS(combined, file = "ALL_INT_STR.rds")
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+

################################################################################################
#CELL TYPE EXPRESSION Integration of each level of filtering to ensure cell groupings are still meaningful as cell partioning.
#NOTE Integration should be run on each level of filtering. This would be combined, combined2, combined3. These are saved and the right point to run the steps below. Update as needed to run each through these steps.
combined <- readRDS("./PREINT_HIGHEXP3_INTV2_STR.rds") #Read multispeices object in
combined <- SCTransform(combined, ncells = 15000) #perform SCT normalization and transformation
combined <- RunPCA(combined, assay = "SCT", npcs=100) #Run dimensionality
combined <- IntegrateLayers(object = combined, method = HarmonyIntegration, normalization.method = "SCT", assay = "SCT", verbose = F, orig.reduction = "pca", new.reduction = 'harmony') 
combined <- RunUMAP(combined, assay = "SCT", dims = 1:40, reduction = "harmony") #Projection visulaization post integration
#Diagnostic plots post integration
#DimPlot(combined, reduction = "umap", group.by = "orig.ident", raster=FALSE, alpha = 0.1) + theme(text = element_text(size = 28)) 
#DimPlot(combined, reduction = "umap", group.by = "CT_toplevel", raster=FALSE) + theme(text = element_text(size = 28))  + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
combined$CT_toplevel <- ifelse(combined$orig.ident == 'PIG', combined$CT, combined$CT_toplevel)
combined$CT_toplevel2 <- collapse_to_toplevel(combined$CT_toplevel)
combined$CT_toplevel_species <- paste0(combined$species, "_", combined$CT_toplevel2)
#DimPlot(combined, reduction = "umap", group.by = "CT_toplevel2", raster=FALSE) + theme(text = element_text(size = 28))  + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
png("INT_species_UMAP3_STR.png", width=1200, height=500)
DimPlot(combined, reduction = "umap", group.by = "species", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
png("INT_CT_UMAP3_STR.png", width=1200, height=500)
DimPlot(combined, reduction = "umap", group.by = "CT_toplevel2", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
dev.off()
png("INT_CT_species_UMAP3_STR.png", width=1200, height=500)
DimPlot(combined, reduction = "umap", group.by = "CT_toplevel_species", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
dev.off()
#Save the objects after running.
saveRDS(combined, file = "HIGHEXP3_INT_STR.rds")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+

################################################################################################


combined <- readRDS("./ALL_INT_STR.rds") #Read multispeices object in
combined$CT_species <- paste0(combined$species, "_", combined$CT)
DimPlot(combined, reduction = "umap", group.by = "species", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
DimPlot(combined, reduction = "umap", group.by = "CT_toplevel2", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
DimPlot(combined, reduction = "umap", group.by = "CT_toplevel_species", label = F, raster=FALSE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + guides(color = guide_legend(override.aes = list(size=4), ncol=1))

library(Seurat)
library(dplyr)

marker_lists <- list()

for(sp in unique(combined$species)) {
    message("Species: ", sp)
    obj.sp <- subset(combined, subset = species == sp)
    Idents(obj.sp) <- "Cluster"
    cts <- levels(Idents(obj.sp))
    for(ct in cts) {
        ncells <- sum(Idents(obj.sp) == ct)
        if(ncells < 10) {
            message("Skipping ", ct)
            next
        }
        message("  ", ct)
        markers <- FindMarkers(obj.sp, ident.1 = ct, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.10)
        sig <- markers %>%
            filter(p_val_adj < 0.05, avg_log2FC > 0.25)
        sig$gene <- rownames(sig)
        sig <- sig %>%
            mutate(score = avg_log2FC * -log10(pmax(p_val_adj, 1e-300))) %>%
            arrange(desc(score)) %>%
            slice_head(n = 200)
        marker_lists[[paste(sp, ct, sep="_")]] <- sig
    }
}

names(marker_lists)


jaccard <- function(df1, df2){
    g1 <- df1$gene
    g2 <- df2$gene
    if(length(g1) == 0 || length(g2) == 0)
        return(NA)
    length(intersect(g1, g2)) /
        length(union(g1, g2))
}

logfc_cor <- function(df1, df2){
    shared <- intersect(df1$gene, df2$gene)
    if(length(shared) < 3)
        return(NA)
    x <- df1$score[match(shared, df1$gene)]
    y <- df2$score[match(shared, df2$gene)]
    cor(x, y, method = "spearman") #ranks more stable than actual FC
}

combined_score <- function(jac, cor){
    if(is.na(jac) || is.na(cor))
        return(NA)
    jac * ((cor + 1) / 2)
}

combined_score <- function(jac, cor){
    if(is.na(jac) || is.na(cor))
        return(NA)
    jac * cor
}

sigs <- names(marker_lists)

jac_mat <- matrix(NA, length(sigs), length(sigs), dimnames=list(sigs,sigs))

fc_mat <- matrix(NA, length(sigs), length(sigs), dimnames=list(sigs,sigs))

combo_mat <- matrix(NA, length(sigs), length(sigs), dimnames=list(sigs,sigs))

for(i in seq_along(sigs)){
    for(j in seq_along(sigs)){
        jac <- jaccard(marker_lists[[i]], marker_lists[[j]])
        fc <- logfc_cor(marker_lists[[i]], marker_lists[[j]])
        combo <- combined_score(jac, fc)
        jac_mat[i,j] <- jac
        fc_mat[i,j] <- fc
        combo_mat[i,j] <- combo
    }
}
saveRDS(marker_lists, file="12K_CT_DEGs.rds")
write.csv(jac_mat, "12K_jac_DEG.csv", quote=F)
write.csv(fc_mat, "12K_fc_DEG.csv", quote=F)
write.csv(combo_mat, "12K_combo_DEG.csv", quote=F)


#PigMSN to HMSN
length(intersect(marker_lists[[6]]$gene, marker_lists[[18]]$gene)) / length(union(marker_lists[[6]]$gene, marker_lists[[18]]$gene))
#PigMSN to HMSN
length(intersect(marker_lists[[39]]$gene, marker_lists[[18]]$gene)) / length(union(marker_lists[[39]]$gene, marker_lists[[18]]$gene))