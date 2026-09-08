################################################################################################
#The QC steps used in our analyses. This gets us to the raw counts tables used for analyses.
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
#sort(sapply(ls(), function(x) format(object.size(get(x)), unit = 'auto')))

#Setting up MT gene list
mt_file <- read.csv("../Alignement/Align/S25_105-3-Thalamus_THA_mt_genes.txt", header = FALSE)
mt_g <- as.character(mt_file[1,])
print(paste0(mt_g," MT Genes"))
set.seed(9578)

#Read count matrix into R
samples <- c("../Alignement/Align/${sample}/outs/S23_105-2-Thalamus_THA_seurat.h5","../Alignement/Align/${sample}/outs/S25_105-3-Thalamus_THA_seurat.h5","../Alignement/Align/${sample}/outs/S27_105-10-Thalamus_THA_seurat.h5","../Alignement/Align/${sample}/outs/S22_105-2-Striatum_STR_seurat.h5","../Alignement/Align/${sample}/outs/S24_105-3-Striatum_STR_seurat.h5","../Alignement/Align/${sample}/outs/S26_105-10-Striatum_STR_seurat.h5")

#Read Alignement Stats into R
otherqc <- c("../Alignement/Align/${sample}/outs/Contamintion/S23_105-2-Thalamus_THA_collapsed.tsv","../Alignement/Align/${sample}/outs/Contamintion/S25_105-3-Thalamus_THA_collapsed.tsv","../Alignement/Align/${sample}/outs/Contamintion/S27_105-10-Thalamus_THA_collapsed.tsv","../Alignement/Align/${sample}/outs/Contamintion/S22_105-2-Striatum_STR_collapsed.tsv","../Alignement/Align/${sample}/outs/Contamintion/S24_105-3-Striatum_STR_collapsed.tsv","../Alignement/Align/${sample}/outs/Contamintion/S26_105-10-Striatum_STR_collapsed.tsv")

#Sample IDs
ids <- c("S22_105-2-Striatum_STR", "S23_105-2-Thalamus_THA", "S24_105-3-Striatum_STR", "S25_105-3-Thalamus_THA", "S26_105-10-Striatum_STR", "S27_105-10-Thalamus_THA")
prefix <- "9823"

iter <- length(samples)
out1_list <- vector(mode='list', length=iter)
auc_results <- data.frame()

#MAD filtering Function
mad_filter <- function(x, nmads = 3, log_transform = FALSE) {
  if (log_transform) x <- log10(x + 1) 
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, center = med, constant = 1, na.rm = TRUE)
  lower <- med - nmads * mad_val
  upper <- med + nmads * mad_val
  return(list(lower = lower, upper = upper))
}

#Violin Plot Function because the Seurat one was doing funky stuff in this environement
VlnPlot_meta_fixed <- function(seurat_obj, features, group.by = NULL, cols = NULL, ncol = NULL, fill = "grey80") {
  # Verify requested variables
  missing <- setdiff(features, colnames(seurat_obj@meta.data))
  if (length(missing) > 0) {
    stop("Missing from meta.data: ", paste(missing, collapse = ", "))
  }
  df <- seurat_obj@meta.data[, features, drop = FALSE]
  if (!is.null(group.by)) {
    if (!(group.by %in% colnames(seurat_obj@meta.data))) {
      stop("group.by variable not found in meta.data")
    }
    df[[group.by]] <- seurat_obj@meta.data[[group.by]]
  } else {
    df[[group.by]] <- "All"
    group.by <- group.by %||% "Group"
  }
  df_long <- df |>
    tidyr::pivot_longer(cols = all_of(features), names_to = "Feature", values_to = "Value")
  p <- ggplot(df_long, aes(x = !!sym(group.by), y = Value, fill = !!sym(group.by))) + geom_violin(scale = "width", trim = TRUE) + facet_wrap(~Feature, scales = "free_y", ncol = ncol) + theme_bw(base_size = 12) + theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + labs(x = NULL, y = NULL)
  if (!is.null(cols)) p <- p + scale_fill_manual(values = cols)
  p
}

estimate_doublet_rate <- function(ncells){
  rates <- data.frame(
    cells = c(1000, 2000, 4000, 6000, 8000, 10000, 20000),
    rate  = c(0.004,0.008,0.02,0.04,0.06,0.08,0.15)
  )
  approx(x = rates$cells, y = rates$rate, xout = ncells, rule = 2)$y
}

for (s in 1:iter) {
	id <- Read10X_h5(filename = samples[s], use.names = TRUE) #Read .h5 file into memory
	id <- CreateSeuratObject(counts = id, project = ids[s], min.cells = 3, min.features = 100) #create seurat object (SO) adding sample ID and filtering genes found in less than 3 cells and cell with counts from less than 100 genes
	id <- RenameCells(object = id, add.cell.id = ids[s]) #Add sample ID to each cell
	set.seed(9578)
	id.sce <- as.SingleCellExperiment(id)
	# Estimate expected doublet rate (~1% per 1000 cells)
	ncells <- ncol(id)
	dbr_est <- estimate_doublet_rate(ncells)
	dbr_sd <- dbr_est * 0.1
	message(paste0("Sample: ", ids[s], " | Cells: ", ncells, " | Expected doublet rate: ", round(dbr_est, 4)))
	id.sce <- scDblFinder(id.sce, dbr = dbr_est, dbr.sd = dbr_sd)
	#Add doublet metrics to SO
	id$scDblFinder.score <- id.sce$scDblFinder.score
	id$scDblFinder.class <- id.sce$scDblFinder.class
	id$dbr_est <- dbr_est
	id$dbr_sd  <- dbr_sd
	id.sce <- NULL #memory clean
	id[["percent.mt"]] <- PercentageFeatureSet(id, features = mt_g) #Add MT metrics to SO 
	id[["percent.rp"]] <- PercentageFeatureSet(id, pattern = "^(RP[LS]|rp[ls])") #Add RPLS metrics to SO
	df <- read_tsv(otherqc[s]) #Read alignment metrics into memory
	df <- df %>%
		mutate(CB_full = paste0(ids[s], "_", CB)) #Add sample ID to the cell barcode in the table to match SO cell barcode
	#ensure barcodes match
	common_barcodes <- intersect(colnames(id), df$CB_full) 
	SO_metrics <- df %>%
		filter(CB_full %in% common_barcodes) %>%
		column_to_rownames("CB_full")
	SO_metrics <- SO_metrics[,-1]
	# Add to Seurat metadata
	id <- AddMetaData(id, SO_metrics)

	#QC violin plot
	png(paste0(ids[s],"_Pre_QC2_.png"), width=1000, height=700)
	print(VlnPlot_meta_fixed(id, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "scDblFinder.score", "pct_exonic", "pct_intronic", "pct_intergenic"), group.by = "orig.ident", ncol = 4))
	dev.off()
	# Apply MAD thresholds
	qc.df <- id@meta.data

	# nFeature_RNA (log-transformed, flag low and high)
	nfilt <- mad_filter(qc.df$nFeature_RNA, nmads = 3, log_transform = TRUE)
	qc.df$MAD_nFeature <- ifelse(log10(qc.df$nFeature_RNA+1) < nfilt$lower, "Low_nFeature_RNA", ifelse(log10(qc.df$nFeature_RNA+1) > nfilt$upper, "High_nFeature_RNA","Pass"))
	print(paste0(iter," Samples in the cohort"))

	# nCount_RNA (log-transformed)
	nc_filt <- mad_filter(qc.df$nCount_RNA, nmads = 3, log_transform = TRUE)
	qc.df$MAD_nCount <- ifelse(log10(qc.df$nCount_RNA+1) < nc_filt$lower, "Low_nCount", ifelse(log10(qc.df$nCount_RNA+1) > nc_filt$upper, "High_nCount","Pass"))

	# percent.mt (not log-transformed)
	mt_filt <- mad_filter(qc.df$percent.mt, nmads = 3)
	qc.df$MAD_MT <- ifelse(qc.df$percent.mt > mt_filt$upper, "High_MT","Pass")
					
	# percent.intergenic (not log-transformed)
	inter_filt <- mad_filter(qc.df$pct_intergenic, nmads = 3)
	qc.df$MAD_pct_inter <- ifelse(qc.df$pct_intergenic > inter_filt$upper, "High_inter","Pass")

	# doublet scores
	qc.df$MAD_dbl <- ifelse(qc.df$scDblFinder.class == "doublet", "Doublet", "Pass")

	# Combine into single call
	qc.df$QC_MAD <- apply(qc.df[,c("MAD_nFeature","MAD_nCount","MAD_MT","MAD_pct_inter","MAD_dbl")],1, function(x) ifelse(any(x!="Pass"), paste(x[x!="Pass"], collapse=";"), "Pass"))
	auc_results <- rbind(auc_results, data.frame(sample = ids[s], MAD_nFeature_low = nfilt$lower, MAD_nFeature_upper = nfilt$upper, MAD_nCount_RNA_low = nc_filt$lower, MAD_nCount_RNA_upper = nc_filt$upper, MAD_MT = mt_filt$upper, MAD_pct_inter = inter_filt$upper))
	# Add back to Seurat object
	id@meta.data <- qc.df
	out1_list[[s]] <- id
}
write.csv(auc_results, 'filter_metrics.csv', quote=F)
#Merge individual SOs together into one object.
combined <- merge(out1_list[[1]], y = c(out1_list[[2]], out1_list[[3]], out1_list[[4]], out1_list[[5]], out1_list[[6]]))

clus_res <- as.character(1.2)	#resolution parameter
iter <- length(samples)	#iteration value
out_list_int <- vector(mode='list', length=iter)	#list object of seurat objects for integration
#QC plots
png('Filter_metrics_Vln.png', width=1000, height=700) 
VlnPlot_meta_fixed(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "scDblFinder.score", "pct_exonic", "pct_intronic", "pct_intergenic"), group.by = "orig.ident", ncol = 4)
dev.off()
png('Filter_metrics_Pass_Vln.png', width=1000, height=700) 
VlnPlot_meta_fixed(subset(combined, subset = QC_MAD == 'Pass'), features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "scDblFinder.score", "pct_exonic", "pct_intronic", "pct_intergenic"), group.by = "orig.ident", ncol = 4)
dev.off()

#Normalization
combined <- SCTransform(combined, vars.to.regress = c("percent.mt", "pct_intergenic"), variable.features.n = 3000)
saveRDS(combined, paste0(prefix,"_pre_QC.rds")) #Happy save point
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


