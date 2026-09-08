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
#We First must make objects to do analyses with. 

#####CONVERSION table
ensembl <- useEnsembl("ensembl", version = 112) #Alignement consistent version.
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl) #Pull human mart
pig <- useDataset("sscrofa_gene_ensembl", mart = ensembl) #Pull pig mart
mouse <- useDataset("mmusculus_gene_ensembl", mart = ensembl) #Pull mouse mart
# Human <-> Pig
orth_hp <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "sscrofa_homolog_ensembl_gene", "sscrofa_homolog_associated_gene_name", "sscrofa_homolog_orthology_type"), mart = human)
# Human <-> Mouse
orth_hm <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"), mart = human)
orth_hp_1to1 <- subset(orth_hp, sscrofa_homolog_orthology_type == "ortholog_one2one") #Only keep true one to one orthologs
orth_hm_1to1 <- subset(orth_hm, mmusculus_homolog_orthology_type == "ortholog_one2one") #Only keep true one to one orthologs
orth_merge <- merge(orth_hp_1to1, orth_hm_1to1, by = c("ensembl_gene_id", "external_gene_name")) #Merged list of 121s
orth3 <- orth_merge[, c("external_gene_name", "ensembl_gene_id", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene")] #Only relevant columns
colnames(orth3) <- c("Human_Symbol", "Human_Ensembl", "Pig_Ensembl", "Mouse_Ensembl") #Column names.

#####PIG
SO <- readRDS("../Processing/9823_POST_INT_MAD.rds") #Filtered and integrated Seruat Object.
Labels <- read.csv("../Processing/Full_Labels.csv")
SO$Cluster <- Labels$Cluster
Labels$Final_celltype2 <- ifelse(Labels$Final_celltype == "Unclassified_Neuron", "Neuron_Unclassified", Labels$Final_celltype)
SO$CT <- Labels$Final_celltype2
SO$tissue <- Labels$Tissue
DefaultAssay(SO) <- "RNA" #Switch back to original counts 
SO <- JoinLayers(SO)

SOSTR <- subset(SO, tissue == 'Striatum')
expr_matrix <- GetAssayData(SOSTR, layer = 'counts') #extract raw counts
symbol2ensembl <- getBM(attributes=c("external_gene_name","ensembl_gene_id"), mart=pig) #Get the gene ID that matches with the gene names
gene2ensembl <- symbol2ensembl$ensembl_gene_id[match(rownames(expr_matrix), symbol2ensembl$external_gene_name)] #Match pig sybols to gene ID
orth_idx <- match(gene2ensembl, orth3$Pig_Ensembl) #Make an index of the gene IDs that match the human conversion table
final_gene_ids <- ifelse(!is.na(orth_idx), orth3$Human_Ensembl[orth_idx], gene2ensembl) #if missing a human ortholog set as the pig gene ID
rownames(expr_matrix) <- make.unique(final_gene_ids) #replace and make unique the rownames of the pig count matrix
outmat <- as.data.frame(expr_matrix) #output class is data.frame
data.table::fwrite(outmat, file="PIG_counts_STR.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write new count matrix
cell_labels <- colnames(expr_matrix)
metadata <- data.frame(cell_label = cell_labels, orig.ident = "PIG", CT = paste0("PBA_", SOSTR$CT), Cluster = SOSTR$Cluster) #Make simple metadata file with cell annotations marked as pig
data.table::fwrite(metadata, file="PIG_meta_STR.csv", quote=F, row.names=FALSE, col.names=TRUE) #Write metadata


#####SET UP COUNT MATRICIES AND METADATA

#####BUILD COUNT MATRICIES AND SEURAT OBJECTS
PBA <- as.data.frame(fread("./PIG_counts_STR.csv"))
PBAmeta <- as.data.frame(fread("./PIG_meta_STR.csv"))
rownames(PBA) <- PBA$V1
PBA <- PBA[,-1]
rownames(PBAmeta) <- PBAmeta$cell_label
PBAmeta <- PBAmeta[,-1]
PBA[1:5,1:5]
head(PBAmeta)

#HUMAN
HBA_Counts <- as.data.frame(fread("../Processing/AnnotationRef/HBA_counts_STR.csv"))
HBA_meta <- as.data.frame(fread("../Processing/AnnotationRef/HBA_meta_STR.csv"))
rownames(HBA_Counts) <- HBA_Counts$V1 #Make gene ID the rownames
HBA_Counts <- HBA_Counts[,-1] #Remove gene ID column
rownames(HBA_meta) <- HBA_meta[,1] #Make Cell ID the rownames
colnames(HBA_meta) <- HBA_meta[1,] #make first row the column names 
HBA_meta <- HBA_meta[-1,] #Remove if nessicary
HBA_meta <- HBA_meta[,-1] #Remove if nessicary
HBA_Counts[1:5,1:5] #check orientation and structure
head(HBA_meta) #check orientation and structure
#MOUSE
MBA_Counts <- as.data.frame(fread("../Processing/AnnotationRef/MBA_counts.csv"))
MBA_meta <- as.data.frame(fread("../Processing/AnnotationRef/MBA_meta.csv"))
rownames(MBA_Counts) <- MBA_Counts$V1 #Make gene ID the rownames
MBA_Counts <- MBA_Counts[,-1] #Remove gene ID column
rownames(MBA_meta) <- MBA_meta[,1] #Make Cell ID the rownames
colnames(MBA_meta) <- MBA_meta[1,] #make first row the column names 
MBA_meta <- MBA_meta[-1,] #Remove if nessicary
MBA_meta <- MBA_meta[,-1] #Remove if nessicary
MBA_Counts[1:5,1:5] #check orientation and structure
head(MBA_meta) #check orientation and structure
#
common_genes <- Reduce(intersect, list(rownames(PBA), rownames(HBA_Counts), rownames(MBA_Counts))) #Identify the one-to-one orthologs
#Filter housekeeping genes
HKG <- fread("../Processing/AnnotationRef/Human_Mouse_Common.csv") #Identify where this comes from https://housekeeping.unicamp.br/?download
orth_idx <- match(HKG$Human, orth3$Human_Symbol) #Match gene symbol to human gene ID
HKG$HEID <- orth3$Human_Ensembl[orth_idx] #Pull out those with a match
HKG <- na.omit(HKG$HEID) #Make a vector of only the matching genes

#filter house keeping genes and identify MT and RSP genes
common_genes <- common_genes[-which(common_genes %in% HKG)] #Remove the Housekeeping genes from the ortholog list
ensembl <- useEnsembl("ensembl", version = 112)
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
mt_ribo <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "gene_biotype", "external_gene_name"), mart = human) #Get conversion table to identify MT and RSP genes
mt_genes <- mt_ribo$ensembl_gene_id[mt_ribo$chromosome_name == "MT"] # Identify mitochondrial genes
ribo_genes <- mt_ribo$ensembl_gene_id[grepl("^RPL|^RPS", mt_ribo$external_gene_name)] # Identify ribosomal genes (RPL / RPS)
bad_genes <- unique(c(mt_genes, ribo_genes)) #Single list of genes
common_genes <- common_genes[-which(common_genes %in% bad_genes)] #final set of genes to use
length(common_genes)

#Remove genes from count matrix to get the final analysis ready count matrix.
PBA <- PBA[rownames(PBA) %in% common_genes, ]
HBA_Counts <- HBA_Counts[rownames(HBA_Counts) %in% common_genes, ]
MBA_Counts <- MBA_Counts[rownames(MBA_Counts) %in% common_genes, ]

PBA_obj <- CreateSeuratObject(counts = PBA, meta.data = PBAmeta) #Generate Seurat Obeject
PBA_obj <- NormalizeData(PBA_obj, normalization.method = "LogNormalize", scale.factor = 10000) #Normalize counts
PBA_obj <- ScaleData(PBA_obj) #Scale counts to account for library size differences
PBA_obj <- FindVariableFeatures(PBA_obj) #Identify variable genes
PBA_obj <- RunPCA(PBA_obj, npcs=100) #Determine dimensionality
HBA_obj <- CreateSeuratObject(counts = HBA_Counts, meta.data = HBA_meta)
HBA_obj <- NormalizeData(HBA_obj, normalization.method = "LogNormalize", scale.factor = 10000)
HBA_obj <- ScaleData(HBA_obj)
HBA_obj <- FindVariableFeatures(HBA_obj)
HBA_obj <- RunPCA(HBA_obj, npcs=100)
MBA_obj <- CreateSeuratObject(counts = MBA_Counts, meta.data = MBA_meta)
MBA_obj <- NormalizeData(MBA_obj, normalization.method = "LogNormalize", scale.factor = 10000)
MBA_obj <- ScaleData(MBA_obj)
MBA_obj <- FindVariableFeatures(MBA_obj)
MBA_obj <- RunPCA(MBA_obj, npcs=100)

# Add species labels
PBA_obj$species <- "P"
HBA_obj$species  <- "H"
MBA_obj$species  <- "M"


HBA_obj$CT_toplevel_species <- paste0(HBA_obj$species, "_", HBA_obj$CT_toplevel) #Add species back to annotation at top level
MBA_obj$CT_toplevel_species <- paste0(MBA_obj$species, "_", MBA_obj$CT_toplevel) #Add species back to annotation at top level

(table(PBA_obj$CT)/sum(table(PBA_obj$CT)))*100
(table(HBA_obj$CT_toplevel)/sum(table(HBA_obj$CT_toplevel)))*100
(table(MBA_obj$CT_toplevel)/sum(table(MBA_obj$CT_toplevel)))*100

saveRDS(PBA_obj, file = "PBA_INT_STR.rds")
saveRDS(HBA_obj, file = "HBA_INT_STR.rds")
saveRDS(MBA_obj, file = "MBA_INT_STR.rds")
combined <- merge(PBA_obj, y = list(HBA_obj, MBA_obj), add.cell.ids = c("PBA", "HBA", "MBA"))
saveRDS(combined, file = "ALL_INT_STR.rds") #Save for later integration
#These species objects are the starting point for the CELL TYPE EXPRESSION, MARKER GENE EXPRESSION, and PRINCIPAL COMPOENENT EXPRESSION analyses.
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
PBA_obj <- readRDS("./PBA_INT_STR.rds")
HBA_obj <- readRDS("./HBA_INT_STR.rds")
MBA_obj <- readRDS("./MBA_INT_STR.rds")

HBA_obj$CT <- paste0(HBA_obj$species, "_", HBA_obj$CT) #Add species back to annotation at top level
MBA_obj$CT <- paste0(MBA_obj$species, "_", MBA_obj$CT) #Add species back to annotation at top level
PBA_obj$CT <- paste0(PBA_obj$species, "_", PBA_obj$CT) 


################################################################################################
#CELL TYPE EXPRESSION COMPARIOSN. So this is the average expression of most all to most expressed genes per cell grouping. 
#####BULK EXPRESSION COMPARIOSN AND FILTERING
get_avg_expr <- function(obj) rowSums(GetAssayData(obj, layer = "counts")) #Function to pull normalized counts per gene
avg_PBA <- get_avg_expr(PBA_obj) #Pig
avg_HBA <- get_avg_expr(HBA_obj) #Human
avg_MBA <- get_avg_expr(MBA_obj) #Mouse
shared <- Reduce(intersect, list(names(avg_PBA), names(avg_HBA), names(avg_MBA))) # Use rownames (genes)

P <- avg_PBA[shared]
H <- avg_HBA[shared]
M <- avg_MBA[shared]

Pdge2 <- DGEList(P)
Pdge2 <- calcNormFactors(Pdge2)
Plogcpm <- cpm(Pdge2, log=TRUE, prior.count = 1)
Hdge2 <- DGEList(H)
Hdge2 <- calcNormFactors(Hdge2)
Hlogcpm <- cpm(Hdge2, log=TRUE, prior.count = 1)
Mdge2 <- DGEList(M)
Mdge2 <- calcNormFactors(Mdge2)
Mlogcpm <- cpm(Mdge2, log=TRUE, prior.count = 1)

cors <- c(
  H_P = cor(Hlogcpm, Plogcpm, method = "pearson"),
  H_M = cor(Hlogcpm, Mlogcpm, method = "pearson"),
  P_M = cor(Plogcpm, Mlogcpm, method = "pearson")
)

#PLOTTING to pick low expressed threshold
m_PBA <- data.frame(PBAExpress = Plogcpm) #Data frame of the expression values per gene
m_PBA$Gene <- rownames(m_PBA) #Rownames become a new column in this new matrix
m_HBA <- data.frame(HBAExpress = Hlogcpm)
m_HBA$Gene <- rownames(m_HBA)
m_MBA <- data.frame(MBAExpress = Mlogcpm)
m_MBA$Gene <- rownames(m_MBA)

combinede <- merge(m_PBA, m_HBA, by = "Gene") #Merge the matricies by gene ID: Human with pig
combinede <- merge(combinede, m_MBA, by = "Gene") #Then add mouse
rownames(combinede) <- combinede$Gene #Make the rownames the gene list
combinede <- combinede[,-1] #Remove the gene column
combinede <- as.matrix(combinede) #Convert to matrix for easy math
GeneMean <- rowMeans(combinede) #Mean across the three species
gene_sd <- rowSds(combinede) #Standard deviation across the three species

df_expr <- data.frame(mean_expr = GeneMean) #Set up plotable device
#Density plot of expression
p_volcano <- ggplot(df_expr, aes(x = mean_expr)) + geom_density(fill = "steelblue", alpha = 0.4) + geom_vline(xintercept = c(1, 2, 3), linetype = "dashed") + labs(title = NULL, x = "Max log-normalized expression", y = "Density") + theme_classic()
ggsave("Bulk_express_density_STR.png", p_volcano)

#Expression vs variability
df_var <- data.frame(mean_expr = GeneMean, sd_expr   = gene_sd)
p_volcano <- ggplot(df_var, aes(x = mean_expr, y = sd_expr)) + geom_point(alpha = 0.15, size = 0.5) + geom_smooth(method = "loess", color = "red") + geom_vline(xintercept = c(1, 2, 3), linetype = "dashed") + labs(title = NULL, x = "Mean log-normalized expression", y = "SD across species (log scale)") + theme_classic()
ggsave("Bulk_express_v_varV2_STR.png", p_volcano)

combinede <- cbind(Plogcpm, Hlogcpm, Mlogcpm)

gene_mean <- rowMeans(combinede)
gene_var  <- matrixStats::rowVars(combinede)

df_mv <- data.frame(
  mean = gene_mean,
  variance = gene_var
)

p <- ggplot(df_mv, aes(mean, variance)) +
  geom_point(alpha = 0.2, size = 0.4) +
  geom_smooth(method = "loess", color = "red") + geom_vline(xintercept = c(1, 2, 3), linetype = "dashed") + 
  theme_classic()
ggsave("mean_variance_trend_STR.png", p, width = 6, height = 5)

#Two filtering levels
bulk_keep3 <- GeneMean >= 3 

#####MAJOR CELL TYPE EXPRESSION FILTERING
get_ct_pseudobulk <- function(obj, ct_col="CT_toplevel"){
  counts <- GetAssayData(obj, layer="counts")  #Pull expression data
  ct <- obj[[ct_col]][,1] #Extract metadata column defining cell types
  pb <- sapply(split(seq_len(ncol(counts)), ct), function(idx) Matrix::rowSums(counts[,idx,drop=FALSE])) #Pseudobulk
  return(pb)
}
ct_PBA  <- get_ct_pseudobulk(PBA_obj, ct_col="CT") #PIG
ct_PBA2 <- DGEList(ct_PBA)
ct_PBA2 <- calcNormFactors(ct_PBA2)
ct_PBAlogcpm2 <- cpm(ct_PBA2, log=TRUE, prior.count = 1)
ct_HBA  <- get_ct_pseudobulk(HBA_obj) #HUMAN
ct_HBA2 <- DGEList(ct_HBA)
ct_HBA2 <- calcNormFactors(ct_HBA2)
ct_HBAlogcpm2 <- cpm(ct_HBA2, log=TRUE, prior.count = 1)
ct_MBA  <- get_ct_pseudobulk(MBA_obj) #MOUSE
ct_MBA2 <- DGEList(ct_MBA)
ct_MBA2 <- calcNormFactors(ct_MBA2)
ct_MBAlogcpm2 <- cpm(ct_MBA2, log=TRUE, prior.count = 1)

ct_max_expr <- pmax(
  rowMaxs(ct_PBAlogcpm2),
  rowMaxs(ct_HBAlogcpm2),
  rowMaxs(ct_MBAlogcpm2)
)

CT_keep <- ct_max_expr >= 1


get_ct_detection <- function(obj, ct_col = "CT_toplevel", thresh = 0.7) {
  expr <- GetAssayData(obj, layer = "data") #Pull expression data
  ct   <- obj[[ct_col]][,1]  #Extract metadata column defining cell types
  t(vapply(
    split(seq_len(ncol(expr)), ct), #Create a list of cell indices grouped by cell type
    function(idx)
      Matrix::rowMeans(expr[, idx, drop = FALSE] > thresh), #Subset expression from the cells targeted and calucate average gene expression if the expression is greater than 0.1
    numeric(nrow(expr))
  ))
}

ct_max_detect <- pmax(
  apply(get_ct_detection(PBA_obj, ct_col="CT"),  2, max),
  apply(get_ct_detection(HBA_obj),  2, max),
  apply(get_ct_detection(MBA_obj),  2, max)
)
hist(ct_max_expr, breaks=100) #the maximum mean log-normalized expression of a gene across all cell types and species
hist(ct_max_detect, breaks=100) #maximum fraction of cells (within any one cell type) in which the gene is detected
png("CT_detect_express.png")
plot(ct_max_detect, ct_max_expr, pch = 16, cex = 0.4, col = rgb(0,0,0,0.1), xlab = "ct_max_detect", ylab = "ct_max_expr")
abline(v = 0.3, h = 5, col = "red", lty = 2)
dev.off()

CT_keep <- ct_max_expr >= 5 & ct_max_detect >= 0.3 #List of genes to keep
#This means: Gene is expressed at ≥5 In 30% of cells have log-normalized expression > 0.7 In at least one CT in at least one species

#Repeat but with more stringent filtering on expression
keep_genes3 <- bulk_keep3 & CT_keep
keep_genes3 <- data.frame(genes = keep_genes3)
keep_genes3$names <- rownames(keep_genes3)
analysis_genes3 <- ifelse(keep_genes3$genes == FALSE, NA, keep_genes3$names)
analysis_genes3 <- na.omit(analysis_genes3)
length(analysis_genes3)

PBA_obj3  <- subset(PBA_obj,  features = analysis_genes3)
HBA_obj3  <- subset(HBA_obj,  features = analysis_genes3)
MBA_obj3  <- subset(MBA_obj,  features = analysis_genes3)
combined3 <- merge(PBA_obj3, y = list(HBA_obj3, MBA_obj3))

get_avg_expr <- function(obj) rowSums(GetAssayData(obj, layer = "counts")) #Function to pull normalized counts per gene
avg_PBA <- get_avg_expr(PBA_obj3) #Pig
avg_HBA <- get_avg_expr(HBA_obj3) #Human
avg_MBA <- get_avg_expr(MBA_obj3) #Mouse
shared <- Reduce(intersect, list(names(avg_PBA), names(avg_HBA), names(avg_MBA))) # Use rownames (genes)

P <- avg_PBA[shared]
H <- avg_HBA[shared]
M <- avg_MBA[shared]

Pdge2 <- DGEList(P)
Pdge2 <- calcNormFactors(Pdge2)
Plogcpm <- cpm(Pdge2, log=TRUE, prior.count = 1)
Hdge2 <- DGEList(H)
Hdge2 <- calcNormFactors(Hdge2)
Hlogcpm <- cpm(Hdge2, log=TRUE, prior.count = 1)
Mdge2 <- DGEList(M)
Mdge2 <- calcNormFactors(Mdge2)
Mlogcpm <- cpm(Mdge2, log=TRUE, prior.count = 1)

cors <- c(
  H_P = cor(Hlogcpm, Plogcpm, method = "pearson"),
  H_M = cor(Hlogcpm, Mlogcpm, method = "pearson"),
  P_M = cor(Plogcpm, Mlogcpm, method = "pearson")
)


#Save the objects
saveRDS(combined, file = "PREINT_ALL_INTV_STR.rds")
saveRDS(combined3, file = "PREINT_HIGHEXP3_INTV2_STR.rds")

#####PER CELL TYPE COMPARISON
get_ct_pseudobulk <- function(obj, ct_col="CT"){
  counts <- GetAssayData(obj, layer="counts")  #Pull expression data
  ct <- obj[[ct_col]][,1] #Extract metadata column defining cell types
  pb <- sapply(split(seq_len(ncol(counts)), ct), function(idx) Matrix::rowSums(counts[,idx,drop=FALSE])) #Pseudobulk
  return(pb)
}
avg_HBA  <- get_ct_pseudobulk(HBA_obj3) #HUMAN
avg_HBA <- DGEList(avg_HBA)
avg_HBA <- calcNormFactors(avg_HBA)
avg_HBA <- cpm(avg_HBA, log=TRUE, prior.count = 1)
avg_MBA  <- get_ct_pseudobulk(MBA_obj3) #MOUSE
avg_MBA <- DGEList(avg_MBA)
avg_MBA <- calcNormFactors(avg_MBA)
avg_MBA <- cpm(avg_MBA, log=TRUE, prior.count = 1)
avg_HBAT  <- get_ct_pseudobulk(HBA_obj3, ct_col="CT_toplevel_species") #HUMAN
avg_HBAT <- DGEList(avg_HBAT)
avg_HBAT <- calcNormFactors(avg_HBAT)
avg_HBAT <- cpm(avg_HBAT, log=TRUE, prior.count = 1)
avg_MBAT  <- get_ct_pseudobulk(MBA_obj3, ct_col="CT_toplevel_species") #MOUSE
avg_MBAT <- DGEList(avg_MBAT)
avg_MBAT <- calcNormFactors(avg_MBAT)
avg_MBAT <- cpm(avg_MBAT, log=TRUE, prior.count = 1)

avg_PBA  <- get_ct_pseudobulk(PBA_obj3, 'CT') #PIG
avg_PBA <- DGEList(avg_PBA)
avg_PBA <- calcNormFactors(avg_PBA)
avg_PBA <- cpm(avg_PBA, log=TRUE, prior.count = 1)
clus_PBA  <- get_ct_pseudobulk(PBA_obj3, 'Cluster') #PIG
clus_PBA <- DGEList(clus_PBA)
clus_PBA <- calcNormFactors(clus_PBA)
clus_PBA <- cpm(clus_PBA, log=TRUE, prior.count = 1)


#Function to compare all pairwise counts of cell grouping to cell grouping
align_by_genes <- function(mat1, mat2) {
  common <- intersect(rownames(mat1), rownames(mat2)) #Ensure matching genes
  mat1a <- mat1[common, , drop = FALSE] #Subset counts to shared genes
  mat2a <- mat2[common, , drop = FALSE] #Subset counts to shared genes
  mat2a <- mat2a[rownames(mat1a), , drop = FALSE] # enforce identical order
  stopifnot(identical(rownames(mat1a), rownames(mat2a))) #Sanity check
  list(mat1 = mat1a, mat2 = mat2a) #Output count matricies
}
#Function to do correlations of count matrcies
ct_similarity <- function(mat1, mat2, method = "pearson") {
  aligned <- align_by_genes(mat1, mat2) #Call function to preprare the matcies.
  cor(as.matrix(aligned$mat1), as.matrix(aligned$mat2), method = method) #Correlation analysies 
}
#Correlation analyses
common_genes <- Reduce(intersect, list(rownames(avg_HBA), rownames(avg_MBA), rownames(avg_PBA), rownames(clus_PBA), rownames(avg_HBAT), rownames(avg_MBAT))) #Get list of orthologous genes
avg_all <- cbind(HBA = avg_HBA[common_genes, , drop = FALSE], MBA = avg_MBA[common_genes, , drop = FALSE], PBA = avg_PBA[common_genes, , drop = FALSE], Cluster = clus_PBA[common_genes, , drop = FALSE], HBATOP = avg_HBAT[common_genes, , drop = FALSE], MBATOP = avg_MBAT[common_genes, , drop = FALSE]) #bind the objects together.
#Correlation of all possible cell groupings in all datasets to each other
sim_avg_all <- ct_similarity(avg_all, avg_all)
#Heat map of correlations cell grouping to cell grouping
png('agg_pb_ct_all_heat_STR3.png', width = 2600, height = 2600)
pheatmap(sim_avg_all, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
write.csv(sim_avg_all, "agg_pb_ct_all_corr_STR3.csv", quote=F, row.names=T) #NOTE THAT THIS WAY THE CELL TPYE ASSIGNER FOR THE PIG IS UNKNOWN AND WILL HAVE TO BE ADDED MANUALLY
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
