#####################################################################################
library(pROC)
library(edgeR)
library(Seurat)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(data.table)
library(matrixStats)
setwd("./results/comparison")
################################################################################################
#PRINCIPAL COMPOENENT EXPRESSION
#####
PBA_obj <- readRDS("./PBA_INT_STR.rds")
PBA_obj$CT_clus <- paste0(PBA_obj$Cluster, "_", PBA_obj$CT)
HBA_obj <- readRDS("./HBA_INT_STR.rds")
HBA_obj$CT <- paste0(HBA_obj$species, "_", HBA_obj$CT)
MBA_obj <- readRDS("./MBA_INT_STR.rds")
MBA_obj$CT <- paste0(MBA_obj$species, "_", MBA_obj$CT)
get_ct_pseudobulk <- function(obj, ct_col="CT"){
  counts <- GetAssayData(obj, layer="counts")  #Pull expression data
  ct <- obj[[ct_col]][,1] #Extract metadata column defining cell types
  pb <- sapply(split(seq_len(ncol(counts)), ct), function(idx) Matrix::rowSums(counts[,idx,drop=FALSE])) #Pseudobulk
  return(pb)
}
avg_HBA <- get_ct_pseudobulk(HBA_obj) #HUMAN
avg_HBA <- DGEList(avg_HBA)
avg_HBA <- calcNormFactors(avg_HBA)
avg_HBA <- cpm(avg_HBA, log=TRUE, prior.count = 1)

avg_MBA <- get_ct_pseudobulk(MBA_obj) #MOUSE
avg_MBA <- DGEList(avg_MBA)
avg_MBA <- calcNormFactors(avg_MBA)
avg_MBA <- cpm(avg_MBA, log=TRUE, prior.count = 1)

HBAT_obj <- HBA_obj
HBAT_obj$CT <- HBAT_obj$CT_toplevel_species
avg_HBAT <- get_ct_pseudobulk(HBAT_obj) #HUMAN
avg_HBAT <- DGEList(avg_HBAT)
avg_HBAT <- calcNormFactors(avg_HBAT)
avg_HBAT <- cpm(avg_HBAT, log=TRUE, prior.count = 1)

MBAT_obj <- MBA_obj
MBAT_obj$CT <- MBAT_obj$CT_toplevel_species
avg_MBAT <- get_ct_pseudobulk(MBAT_obj) #MOUSE
avg_MBAT <- DGEList(avg_MBAT)
avg_MBAT <- calcNormFactors(avg_MBAT)
avg_MBAT <- cpm(avg_MBAT, log=TRUE, prior.count = 1)

avg_PBA <- get_ct_pseudobulk(PBA_obj, 'CT') #PIG
avg_PBA <- DGEList(avg_PBA)
avg_PBA <- calcNormFactors(avg_PBA)
avg_PBA <- cpm(avg_PBA, log=TRUE, prior.count = 1)

clus_obj <- PBA_obj
clus_obj$CT <- clus_obj$Cluster
clus_PBA <- get_ct_pseudobulk(clus_obj) #PIG
clus_PBA <- DGEList(clus_PBA)
clus_PBA <- calcNormFactors(clus_PBA)
clus_PBA <- cpm(clus_PBA, log=TRUE, prior.count = 1)



################################################################################################
#PRINCIPAL COMPOENENT EXPRESSION
#####
#Species Seurat Objects
objs <- list(HBA = HBA_obj, MBA = MBA_obj, PBA = PBA_obj, Cluster_PBA = clus_obj)
#Average expression matrices
avg_mats <- list(HBA = avg_HBA, MBA = avg_MBA, PBA = avg_PBA, Cluster_PBA = clus_PBA)

#Function to pull PCs and test their ability to deliniate cell types
pc_auc_by_ct <- function(obj, ct, ct_col = "CT_toplevel", reduction = "pca") {
  message("    Computing PC AUCs")
  scores <- Embeddings(obj, reduction) #Pulls Principal compenent scores from seurat object
  is_ct  <- obj[[ct_col]][, 1] == ct #logical vector of if a cell belongs to a cell type
  if (length(unique(is_ct)) < 2) {
    message("    :( CT ", ct, " has <2 classes — skipping")
    return(NULL)
  }
  auc <- apply(scores, 2, function(pc) { #Prediction of PC on cell type designation
    roc(is_ct, pc, quiet = TRUE)$auc #Calculate AUC (area under ROC curve)
  })
  sort(auc, decreasing = TRUE) #AUCs sorted from best to worst PC
}
#Function to extract top genes driving a set of PCs
pc_genes_by_ct <- function(obj, pcs, n_genes = 200, reduction = "pca") {
  message("    Extracting genes from ", length(pcs), " PCs")
  load <- Loadings(obj, reduction) #Pull gene loadings matrix from the PCA
  genes <- unique(unlist(lapply(pcs, function(pc) { #Collects top genes, combines across PCs and removes duplicates
    names(sort(abs(load[, pc]), decreasing = TRUE))[seq_len(min(n_genes, nrow(load)))] #
  }))) #Takes gene loadings, gets the absolute value, sorts by imporatnce, selects top genes and returns the names
  message("    Total unique genes: ", length(genes))
  genes
}
#Function builds a cell-type-specific gene program for one species
ct_program <- function(obj, ct, ct_col = "CT", auc_cutoff = 0.8, n_genes = 200, reduction = "pca") {
  message("  Cell type: ", ct)
  # Compute PC AUCs
  auc <- pc_auc_by_ct(obj, ct, ct_col, reduction)
  if (is.null(auc)) {
    message("    :( No valid AUCs — skipping")
    return(NULL)
  }
  # Select PCs passing the cutoff
  pcs <- names(auc)[auc >= auc_cutoff]
  message("    PCs passing AUC cutoff: ", length(pcs))
  if (length(pcs) == 0) {
    message("    :( No PCs passed cutoff — skipping")
    return(NULL)
  }
  # Pull top genes from selected PCs
  genes <- pc_genes_by_ct(obj, pcs, n_genes, reduction)
  if (length(genes) < 3) {
    message("    :( Too few genes — skipping")
    return(NULL)
  }
  # Compute **average absolute loading per gene across selected PCs**
  load_mat <- Loadings(obj, reduction)[genes, pcs, drop = FALSE]   # gene x PC
  loadings <- rowMeans(abs(load_mat))                                # collapse across PCs
  # Return a program object
  list(ct = ct, ct_raw = ct, pcs = pcs, auc = auc, genes = genes, loadings = loadings)
}

# Function to extract all programs per species
extract_all_programs <- function(objs, ct_col = "CT", auc_cutoff = 0.8, n_genes = 200, reduction = "pca") {
  programs <- list()
  for (sp in names(objs)) {
    message("\n=== Species: ", sp, " ===")
    obj <- objs[[sp]]
    # Only cell types present in this species
    cts <- unique(obj[[ct_col]][,1])
    message("Found ", length(cts), " cell types in ", sp)
    # Extract programs per cell type
    sp_programs <- setNames(
      lapply(cts, function(ct) {
        ct_program(obj, ct, ct_col = ct_col, auc_cutoff = auc_cutoff, n_genes = n_genes, reduction = reduction)
      }),
      cts
    )
    # Remove NULLs (cell types without valid programs)
    sp_programs <- sp_programs[!sapply(sp_programs, is.null)]
    message("Valid programs extracted: ", length(sp_programs))
    programs[[sp]] <- sp_programs
  }
  programs
}

# Extract programs per species
programs <- extract_all_programs(objs = objs, ct_col = "CT", auc_cutoff = 0.8, n_genes = 50)
#programs <- programs[-c(5,6)]
names(programs) <- names(objs)
canonical_ct <- function(x) {
  gsub("-", "_", x)
}

names(programs$HBA) <- canonical_ct(names(programs$HBA))
colnames(avg_HBA) <- canonical_ct(colnames(avg_HBA))
names(programs$MBA) <- canonical_ct(names(programs$MBA))
colnames(avg_MBA) <- canonical_ct(colnames(avg_MBA))
names(programs$PBA) <- canonical_ct(names(programs$PBA))
colnames(avg_PBA) <- canonical_ct(colnames(avg_PBA))
names(programs$Cluster_PBA) <- canonical_ct(names(programs$Cluster_PBA))
colnames(clus_PBA) <- canonical_ct(colnames(clus_PBA))

# Convert a program object into a named numeric vector: gene -> expression
extract_expression_vector <- function(prog, avg_mats, species, ct) {
  avg <- avg_mats[[species]]
  if (is.null(avg)) return(NULL)
  if (!ct %in% colnames(avg)) {
    return(NULL)
  }
  genes <- prog$genes
  common <- intersect(genes, rownames(avg))
  if (length(common) < 3) return(NULL)
  v <- avg[common, ct, drop = TRUE]
  names(v) <- common
  v
}


compare_programs <- function(prog1, prog2, sp1, ct1, sp2, ct2, avg_mats, method = "pearson", min_genes = 10) {
  if (is.null(prog1) || is.null(prog2)) {
    return(list(n1 = NA_integer_, n2 = NA_integer_, n_shared = 0L, n_used = 0L, jaccard = NA_real_, cor = NA_real_, valid = FALSE ))
  }
  e1 <- extract_expression_vector(prog1, avg_mats, sp1, ct1)
  e2 <- extract_expression_vector(prog2, avg_mats, sp2, ct2)
  empty_cmp <- function(reason = NA_character_) {
    list(n1 = NA_integer_, n2 = NA_integer_, n_shared = 0L, n_used = 0L, jaccard = NA_real_, cor = NA_real_, valid = FALSE, reason = reason)
  }
  if (is.null(e1) || is.null(e2)) {
    return(empty_cmp("expression vector missing"))
  }
  shared_expr <- intersect(names(e1), names(e2))
  n1 <- length(e1)
  n2 <- length(e2)
  n_shared <- length(shared_expr)
  out <- list(n1 = n1, n2 = n2, n_shared = n_shared, n_used = 0L, jaccard = n_shared / (n1 + n2 - n_shared), cor = NA_real_, valid = FALSE)
  if (n_shared < min_genes) return(out)
  if (length(shared_expr) >= min_genes) {
    out$cor <- suppressWarnings(cor(e1[shared_expr], e2[shared_expr], method = method))
    out$n_used <- length(shared_expr)
    out$valid <- TRUE
  }
  out
}
compare_programs <- function(prog1, prog2, sp1, ct1, sp2, ct2, avg_mats, method = "pearson", min_genes = 10) {

  if (is.null(prog1) || is.null(prog2)) {
    return(list(n1 = NA_integer_, n2 = NA_integer_, n_shared = 0L, n_used = 0L, jaccard = NA_real_, cor = NA_real_, valid = FALSE, shared_genes = character(0)))
  }

  e1 <- extract_expression_vector(prog1, avg_mats, sp1, ct1)
  e2 <- extract_expression_vector(prog2, avg_mats, sp2, ct2)

  if (is.null(e1) || is.null(e2)) {
    return(list(n1 = NA_integer_, n2 = NA_integer_, n_shared = 0L, n_used = 0L, jaccard = NA_real_, cor = NA_real_, valid = FALSE, shared_genes = character(0)))
  }
  shared_expr <- intersect(names(e1), names(e2))
  n1 <- length(e1)
  n2 <- length(e2)
  n_shared <- length(shared_expr)
  out <- list(n1 = n1, n2 = n2, n_shared = n_shared, n_used = 0L, jaccard = n_shared / (n1 + n2 - n_shared), cor = NA_real_, valid = FALSE, shared_genes = shared_expr)
  if (n_shared < min_genes) return(out)
  if (length(shared_expr) >= min_genes) {
    out$cor <- suppressWarnings(cor(e1[shared_expr], e2[shared_expr], method = method))
    out$n_used <- length(shared_expr)
    out$valid <- TRUE
  }
  out
}


#Main pipeline function
run_all_program_comparisons <- function(programs, method = "pearson", min_genes = 10, avg_mats = avg_mats) {
  species <- names(programs)
  results <- list()
  for (sp1 in species) {
    for (sp2 in species) {
      message("Comparing ", sp1, " vs ", sp2)
      cts1 <- names(programs[[sp1]])
      cts2 <- names(programs[[sp2]])
      res <- expand.grid(sp1 = sp1, ct1 = cts1, sp2 = sp2, ct2 = cts2, stringsAsFactors = FALSE)
      # initialize columns
      res$n1 <- NA_integer_
      res$n2 <- NA_integer_
      res$n_shared <- NA_integer_
      res$n_used <- NA_integer_
      res$jaccard <- NA_real_
      res$cor <- NA_real_
      res$valid <- FALSE
      res$shared_genes <- vector("list", nrow(res))
      message("....calculating correlation and jaccard")
      for (i in seq_len(nrow(res))) {
        cmp <- compare_programs(
          programs[[res$sp1[i]]][[res$ct1[i]]],
          programs[[res$sp2[i]]][[res$ct2[i]]],
          sp1 = res$sp1[i],
          ct1 = res$ct1[i],
          sp2 = res$sp2[i],
          ct2 = res$ct2[i],
          avg_mats = avg_mats
        )
        res$n1[i] <- cmp$n1
        res$n2[i] <- cmp$n2
        res$n_shared[i] <- cmp$n_shared
        res$n_used[i] <- cmp$n_used
        res$jaccard[i] <- cmp$jaccard
        res$cor[i] <- cmp$cor
        res$valid[i] <- cmp$valid
        res$shared_genes[[i]] <- cmp$shared_genes
      }
      results[[paste(sp1, sp2, sep = "_vs_")]] <- res
    }
  }
  results
}

#Pipeline call
results <- run_all_program_comparisons(programs = programs, method = "pearson", min_genes = 5, avg_mats = avg_mats)
saveRDS(results, "PC_gene_CT_corr_AUC080_FEAT50_STR.rds")
saveRDS(programs, "PC_gene_CT_programs_AUC080_FEAT50_STR.rds")
#Need to check out https://github.com/jiang-junyao/CACIMAR

H2M <- results[["HBA_vs_MBA"]]
H2Msubset <- H2M[H2M$ct1 == "H_Astrocyte" & H2M$ct2 == "M_Astroependymal NN", ]
H2Msubset$shared_genes[[1]]

H2P <- results[["HBA_vs_PBA"]]
H2Psubset <- H2P[H2P$ct1 == "H_Astrocyte" & H2P$ct2 == "PBA_Astrocytes", ]
H2Psubset$shared_genes[[1]]

H2C <- results[["HBA_vs_Cluster_PBA"]]
H2Csubset <- H2C[H2C$ct1 == "H_Astrocyte" & H2C$ct2 == "19", ]
H2Csubset$shared_genes[[1]]

genes <- unique(unlist(list(
  H2Msubset$shared_genes[[1]],
  H2Psubset$shared_genes[[1]],
  H2Csubset$shared_genes[[1]]
)))

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

build_ct_matrix <- function(corr_list, value_col = "cor", use_valid = TRUE) {
  df <- bind_rows(corr_list, .id = "comparison")
  if(use_valid){
    df <- df %>% filter(valid)
  }
  mat <- df %>%
    select(ct1, ct2, value = all_of(value_col)) %>%
    distinct() %>%                # protects against duplicated pairs
    pivot_wider(names_from = ct2, values_from = value) %>%
    column_to_rownames("ct1") %>%
    as.matrix()
  # enforce symmetric ordering
  all_ct <- union(rownames(mat), colnames(mat))
  mat_full <- matrix(NA, nrow = length(all_ct), ncol = length(all_ct), dimnames = list(all_ct, all_ct))
  mat_full[rownames(mat), colnames(mat)] <- mat
  # fill the opposite triangle
  mat_full[lower.tri(mat_full)] <- t(mat_full)[lower.tri(mat_full)]
  return(mat_full)
}
cor_matrix <- build_ct_matrix(results, "cor")
jaccard_matrix <- build_ct_matrix(results, "jaccard")
combined <- (cor_matrix * jaccard_matrix)
write.csv(cor_matrix, "cor_matrix_AUC085_FEAT50_STR.csv", row.names=T, quote=F)
write.csv(jaccard_matrix, "jaccard_matrix_AUC080_FEAT50_STR.csv", row.names=T, quote=F)
write.csv(combined, "Corr_Jaccard_matrix_AUC080_FEAT50_STR.csv", row.names=T, quote=F)

png('cor_matrix_AUC080_FEAT50_heat_STR.png', width = 2600, height = 2600)
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
png('jaccard_matrix_AUC080_FEAT50_heat_STR.png', width = 2600, height = 2600)
pheatmap(jaccard_matrix, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
png('corr_Jaccard_matrix_AUC080_FEAT50_heat_STR.png', width = 2600, height = 2600)
pheatmap(combined, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
