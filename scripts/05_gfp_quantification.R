setwd("./results/GFP")
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
#######
SO <- readRDS("../Processing/9823_POST_INT_MAD.rds") #Filtered and integrated Seruat Object.
Labels <- read.csv("../Processing/Full_Labels.csv")
SO$Cluster <- Labels$Cluster
Labels$Final_celltype2 <- ifelse(Labels$Final_celltype == "Unclassified_Neuron", "Neuron_Unclassified", Labels$Final_celltype)
SO$CT <- Labels$Final_celltype2
SO$tissue <- Labels$Tissue
library(RColorBrewer)
my_colors <- brewer.pal(n = 10, name = "Paired")  # Paired has 12, good contrast
my_colors

# Save UMAP plots with custom colors
png('ANN_Final_celltype_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "CT", label = TRUE, label.size = 4, repel = TRUE, cols = my_colors) +
  theme(text = element_text(size = 16)) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()

png('ANN_Final_celltype_Tiss_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "CT", split.by = 'tissue', label = TRUE, label.size = 4, repel = TRUE, cols = my_colors) +
  theme(text = element_text(size = 16)) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()

Idents(object = SO) <- "Cluster" #Target grouping
DefaultAssay(SO) <- "RNA" #Target counts
SO <- JoinLayers(SO) #merge sample counts together
SO <- NormalizeData(SO) #Normalize counts
SO <- ScaleData(SO, verbose = FALSE) #Scale normalized counts

#####DEG ANALYSIS SET UP
SO[["percent.GFP"]] <- PercentageFeatureSet(SO, features = 'GFP')
SO$GFP_binary <- SO$percent.GFP > 0 # create binary 
SO$GFP_binary <- ifelse(SO$GFP_binary == "TRUE", "GFP+", "GFP-") #Convert to grouping

SO[["percent.TUBB"]] <- PercentageFeatureSet(SO, features = 'TUBB')
SO$TUBB_binary <- SO$percent.TUBB > 0 # create binary 
SO$TUBB_binary <- ifelse(SO$TUBB_binary == "TRUE", "TUBB+", "TUBB-") #Convert to grouping

SO[["percent.NUCKS1"]] <- PercentageFeatureSet(SO, features = 'NUCKS1')
SO$NUCKS1_binary <- SO$percent.NUCKS1 > 0 # create binary 
SO$NUCKS1_binary <- ifelse(SO$NUCKS1_binary == "TRUE", "NUCKS1+", "NUCKS1-") #Convert to grouping

png("GFP_prevelence_UMAP.png", width=800, height=700)
DimPlot(SO, reduction = "umap", group.by = "GFP_binary", label = F, cols = c("pink", "darkgreen")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
png("GFP_intensity_UMAP.png", width=800, height=700)
FeaturePlot(SO, feature = 'GFP', cols = c("lightgrey", "darkgreen")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()

png("TUBB_prevelence_UMAP.png", width=800, height=700)
DimPlot(SO, reduction = "umap", group.by = "TUBB_binary", label = F, cols = c("pink", "blue")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
png("TUBB_intensity_UMAP.png", width=800, height=700)
FeaturePlot(SO, feature = 'TUBB', cols = c("lightgrey", "blue")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()

png("NUCKS1_prevelence_UMAP.png", width=800, height=700)
DimPlot(SO, reduction = "umap", group.by = "NUCKS1_binary", label = F, cols = c("pink", "royalblue4")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
png("NUCKS1_intensity_UMAP.png", width=800, height=700)
FeaturePlot(SO, feature = 'NUCKS1', cols = c("lightgrey", "royalblue4")) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()


SO <- FindVariableFeatures(SO, nfeatures = 3000)
top_genes <- VariableFeatures(SO) #Vector list

#####DEG ANALYSIS
# Compare GFP+ vs GFP- within all cells. Differnece has to be at least greater than |0.25| with genes in at least 10% of the cells sampled. We used the MAST algortithm as it is a GLMM and used expeirmental variables as covariates
for (s in c("GFP", "TUBB", "NUCKS1")) {
	gene_target <- s
	gene_pos <- paste0(gene_target, "+")
	gene_neg <- paste0(gene_target, "-")
	DEG_all <- FindMarkers(SO, group.by = paste0(gene_target, "_binary"), ident.1 = gene_pos, ident.2 = gene_neg, logfc.threshold = 0.25, min.pct = 0.1, features = top_genes, test.use = "MAST", latent.vars = c("orig.ident", "nCount_RNA", "tissue"))
	write.csv(DEG_all, paste0("DEGs_", gene_target, "_MAST.csv"), row.names = TRUE) #Save results

	celltypes <- names(which(table(SO$CT) > 10)) # Cell types with >10 cells
	latent_all <- c("orig.ident", "nCount_RNA", "tissue")

	DEG_list <- lapply(celltypes, function(ct) {
	  message("Processing cluster: ", ct)
	  subset_cells <- subset(SO, CT == ct) # Subset Seruat object to only those cell IDs extracted
	  meta <- subset_cells@meta.data #Extract metadata
	  # Check GFP groups
	  gfp_tab <- table(meta$GFP_binary)
	  if (length(gfp_tab) < 2 || min(gfp_tab) < 3) {
		message("  Skipping ", ct, ": insufficient GFP groups")
		return(NULL)
	  }
	  # Identify valid latent variables SAFELY
	  valid_latent <- latent_all[
		sapply(latent_all, function(v) {
		  v %in% colnames(meta) && length(unique(meta[[v]])) > 1
		})
	  ]
	  if (length(valid_latent) == 0) {
		message("  No valid latent variables for ", ct)
	  } else {
		message("  Using latent vars: ", paste(valid_latent, collapse = ", "))
	  }
	  # Variable features
	  subset_cells <- FindVariableFeatures(subset_cells, nfeatures = 3000)
	  top_genes <- VariableFeatures(subset_cells)
	  # Run DE
	  res <- tryCatch(
		{
		  if (length(valid_latent) > 0) {
			FindMarkers(subset_cells, ident.1 = gene_pos, ident.2 = gene_neg, group.by = paste0(gene_target, "_binary"), test.use = "MAST", features = top_genes, latent.vars = valid_latent)
		  } else {
			FindMarkers(subset_cells, ident.1 = gene_pos, ident.2 = gene_neg, group.by = paste0(gene_target, "_binary"), test.use = "MAST", features = top_genes)
		  }
		},
		error = function(e) {
		  message("  ERROR in cluster ", ct, ": ", e$message)
		  return(NULL)
		}
	  )

	  res
	})

	names(DEG_list) <- paste0("cluster_", celltypes) #Into Results list obejct assign each analysis the cell type label it has.
	names(DEG_list) <- celltypes
	DEG_list <- DEG_list[!sapply(DEG_list, is.null)]
	DEG_list <- lapply(DEG_list, function(x){
	  x %>% rownames_to_column(var = "gene") # add gene column to each DEG result
	})
	DEG_results <- bind_rows(DEG_list, .id = "celltype") # combine into one table
	write.csv(DEG_results, paste0("DEGs_CT_MAST_", gene_target, ".csv"), row.names = FALSE) # save all results
	#####GFP get descriptive statistics
	resultdf <- data.frame(celltype = character(), mean_percentGene = numeric(), sd_percentGene = numeric(), max_percentGene = numeric(), mean_percentGene_Genepos = numeric(), sd_percentGene_Genepos = numeric(), pct_Genepos_cells = numeric(), N = integer(), stringsAsFactors = FALSE)
	clusters <- unique(SO@meta.data[['CT']]) # Get unique clusters
	# Loop through clusters
	for (ct in clusters) {
	  meta_sub <- SO@meta.data[SO@meta.data[['CT']] == ct, ] # Subset metadata for this cluster
	  percentgene <- SO[[paste0("percent.", gene_target)]]
	  genepos_vals <- percentgene[percentgene > 0] # GFP+ metadata
	  # Compute statistics (Same as described above)
	  mean_all <- mean(percentgene, na.rm = TRUE)
	  sd_all <- sd(percentgene, na.rm = TRUE)
	  max_all <- max(percentgene, na.rm = TRUE)
	  mean_Genepos <- if (nrow(genepos_vals) > 0) mean(genepos_vals, na.rm = TRUE) else NA
	  sd_Genepos <- if (nrow(genepos_vals) > 0) sd(genepos_vals, na.rm = TRUE) else NA
	  pct_Genepos <- (nrow(genepos_vals) / nrow(meta_sub)) * 100
	  N <- nrow(meta_sub)
	  # Append results
	  resultdf <- rbind(resultdf, data.frame(celltype = ct, mean_percentGene = mean_all, sd_percentGene = sd_all, max_percentGene = max_all, mean_percentGene_Genepos = mean_Genepos, sd_percentGene_Genepos = sd_Genepos, pct_Genepos_cells = pct_Genepos, N = N, stringsAsFactors = FALSE))
	  message("Finished cell type: ", ct, " (n=", N, " cells, Gene+ cells = ", nrow(genepos_vals), ")") # Optional progress message
	}
	# View and write results results
	resultdf
	write.csv(resultdf, paste0(gene_target, "_CT_stats.csv")
}
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
