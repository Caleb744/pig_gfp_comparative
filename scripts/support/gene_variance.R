SO <- readRDS("../Processing/9823_POST_INT_MAD.rds") #Filtered and integrated Seruat Object.
Labels <- read.csv("../Processing/Full_Labels.csv")
SO$Cluster <- Labels$Cluster
Labels$Final_celltype2 <- ifelse(Labels$Final_celltype == "Unclassified_Neuron", "Neuron_Unclassified", Labels$Final_celltype)
SO$CT <- Labels$Final_celltype2
SO$tissue <- Labels$Tissue
DefaultAssay(SO) <- "RNA" #Switch back to original counts 
SO <- JoinLayers(SO)
DefaultAssay(SO) <- "RNA"
SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000) #Normalize counts
SO <- ScaleData(SO) #Scale counts to account for library size differences


library(Seurat)
library(Matrix)

get_ubiquitous_genes <- function(obj, cluster_col = "seurat_clusters", assay = "RNA", layer = "data", min_pct = 0.1){
  expr <- GetAssayData(obj, assay = assay, layer = layer) #Pull log-normalized expression per cluster
  clusters <- obj[[cluster_col]][,1] #Extract cluster identities
  cluster_levels <- levels(factor(clusters)) #Set the cluster identities to a factor
  pct_list <- lapply(cluster_levels, function(ct){ #Function to get expression of selected gene per cell within a cluster
    cells <- colnames(obj)[clusters == ct] #list of Cells within the selected cluster
    mat <- expr[, cells, drop = FALSE] #for the cells in the cluster extract the expression matrix
    Matrix::rowMeans(mat > 0) #Identifiy the row means of the extracted matrix with expression.
  })
  pct_mat <- do.call(cbind, pct_list) #Bind the columns together
  colnames(pct_mat) <- cluster_levels #Set the matrix columns to cluster identities
  ubiquitous <- rownames(pct_mat)[rowSums(pct_mat >= min_pct) == length(cluster_levels)] #If the sum of percents across the clusters is greater than threshold add to vector of genes.
  list(genes = ubiquitous, pct_matrix = pct_mat) #Output the list of genes and the matrix of percents. 
}

res <- get_ubiquitous_genes(SO, min_pct = 0.10, cluster_col = 'CT')

length(res$genes)
head(res$genes)

variance_profile <- function(obj, genes, cluster_col = "seurat_clusters", assay = "RNA", layer = "data"){
  expr <- GetAssayData(obj, assay = assay, layer = layer)
  clusters <- factor(obj[[cluster_col]][,1])
  cluster_levels <- levels(clusters)
  expr <- expr[genes, , drop = FALSE]
  cluster_sizes <- table(clusters)
  # Cluster means
  cluster_means <- sapply(cluster_levels, function(ct){
    cells <- which(clusters == ct)
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  })
  # Between variance
  between_var <- apply(cluster_means, 1, var)
  # Total variance (per-cell)
  total_var <- apply(expr, 1, var)
  # Within variance (weighted)
  within_var <- sapply(seq_len(nrow(expr)), function(i){
    vals <- as.numeric(expr[i, ])
    cluster_vars <- tapply(vals, clusters, var)
    weighted.mean(cluster_vars, cluster_sizes[names(cluster_vars)], na.rm = TRUE)
  })
  # R²_cluster
  r2_cluster <- between_var / total_var
  tibble::tibble(gene = genes, total_var = total_var, within_var = within_var, between_var = between_var, r2_cluster = r2_cluster, stability_score = 1 - r2_cluster) |> 
    dplyr::arrange(stability_score)
}

var_tbl <- variance_profile(SO, res$genes, cluster_col = 'CT')
head(var_tbl)
write.csv(var_tbl, "Variance_Table_0.10.csv", quote = F, row.names=F)