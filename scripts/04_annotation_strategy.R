#Package list in Seurat conda environement.
library(hdf5r)
library(rhdf5)
library(Seurat)
library(biomaRt)
library(tidyverse)
library(data.table)
library(BPCells) #I had to do a remotes install on this one
library(Azimuth) #I had to do a remotes install on this one 

#####NORMALIZE CELL ANNOTATIONS
#Function which makes all annotations lower cases relatces the species label with nothing, subs dilimiters with spaces. and removes plurals and misspellings.
normalize_ct_labels <- function(x) {
  x <- tolower(x)
  x <- gsub("^(mba|hba|pba)_", "", x)
  x <- gsub("[_-]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- gsub("oligodendrocytes?", "oligodendrocyte", x)
  x <- gsub("astrocytes?", "astrocyte", x)
  x <- gsub("microglia?", "microglia", x)
  x
}
#Major cell type list fitting the various types in each reference into one of these main ones.
TOPLEVEL_CT <- list(
  Astrocytes = c("astrocyte", "astro-", "astro te", "astro nt", "bergmann"),
  Oligodendrocytes_OPC = c("oligodendrocyte", "opc", "cop", "committed oligodendrocyte", "mol", "mfol", "nfol"),
  Microglia = c("microglia", "macrophage", "bam"),
  Endothelial_Vascular = c("endothelial", "endo", "vascular", "vlmc", "pericyte", "peri", "smc", "fibroblast", "abc", "arachnoid"),
  Ependymal_Tanycyte = c("ependymal", "tanycyte", "hypendymal", "choroid"),
  Glutamatergic_Neurons = c("excitatory", "glut", "ex imn", "intratelencephalic", "corticothalamic", "pyramidal", "hippocampal", "thalamic", "amygdala", "dg", "pir"),
  GABAergic_Neurons = c("gaba", "inhibitory", "inh imn", "interneuron", "mge", "cge", "lamp5", "chandelier", "medium spiny", "msn", "str"),
  Hypothalamic_Neurons = c("hypothalamic", "mammillary", "vmh", "arh", "dmh", "pvh", "lha"),
  Brainstem_Nuclei = c("rhombic", "snc", "vta", "ramb", "pag", "prn", "pcg", "ldt", "ppn", "dopa", "foxa1", "splatter"),
  Cerebellar_Neurons = c("cerebellar"),
  Immune = c("immune", "t cell", "tcell", "dendritic", "dc", "b cell", "nk", "monocyte"),
  Neuron_Unclassified = c("Neuron_Unclassified")
)
#Normalize the cell types into the main cell types.
collapse_to_toplevel <- function(x, ontology = TOPLEVEL_CT) {
  x <- normalize_ct_labels(x) #Call the Normalization function
  vapply(x, function(label) {
    hits <- names(ontology)[ #Search in the cell type ontology list (TOPLEVEL_CT) and return the value if matching any part
      vapply(ontology, function(keys) {
        any(grepl(paste(keys, collapse = "|"), label))
      }, logical(1))
    ]
    if (length(hits) == 0) { #If nothing matches or is already NA it is unassigned.
      return("Unassigned")
    }
    if (length(hits) > 1) {
      return(hits[1])  # deterministic, ordered priority
    }
    hits
  }, character(1))
}
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


#####MOUSE ANNOTATION REFERENCE PREP
#Processesed data sets.
data.list <- c() #Empty objects to fill
metadata.list <- c()
#Get large dataset as a Seurat Object written to disk
path <- "./Zeng-Aging-Mouse-10Xv3-raw.h5ad"
data <- open_matrix_anndata_hdf5(path) #Open tunnel to counts on disk not into memory.
data <- convert_matrix_type(data, "uint32_t") #Convert counts to appropriate class
genes <- h5read(path, "var/gene_identifier") #Pull out specific gene list. NOTE: there was an issue with this dataset with additional columns for some genes and not others. So steps are added to correct this
any(grepl("dtype|Name:|gene-symbol", genes)) #Finding the extra column rows which lead to off counts. 
clean_genes <- sub(".*(ENSMUSG[0-9]+).*","\\1",genes) #ensure uniqness if needed
any(grepl("dtype|Name:|gene-symbol", clean_genes)) #Finding the extra column rows which lead to off counts in the cleaned gene list as a sanity check (should be 0 now)
clean_genes_unique <- make.unique(clean_genes) #ensure uniqness
dimnames(data)[[1]] <- clean_genes_unique #Add gene names back to count matrix
write_matrix_dir(mat = data, dir = "./ZENG_MBA_RAW_FIXED") #Write to disk in manner readable for Seurat. 
mat <- open_matrix_dir("./ZENG_MBA_RAW_FIXED") #Stream the count Matrix
any(grepl("dtype|Name:|gene-symbol", rownames(mat))) #Sanity check that cleaning worked.
length(unique(rownames(mat))) == nrow(mat) #Sanity check that cleaning worked.
#Prepare the metadata
metadata.list[[1]] <- LoadH5ADobs(path = path) #collect metadata from H5 file.
zeng_clus <- read.csv("./cell_cluster_annotations.csv") #Collect metadata from annotation file provided by Allen Institute.
order <- c(1:1162565) #Order of cells in full dataset.
zeng_clus$order <- c(1:1162565) #add to metadata.
metadata.list[[2]] <- zeng_clus #add to metadata list object. 
metadata <- merge(metadata.list[[1]], metadata.list[[2]], by="cell_label") #Merge objects by cell ID.
rownames(metadata) = metadata$cell_label
metadata <- metadata[order(metadata$order),] #Reorder based on original order so it matches the counts matrix. 
#Generate seurat Object.
neuron.object <- CreateSeuratObject(counts = mat, meta.data = metadata) 
cluster_name <- as.data.frame(neuron.object$cluster_name) #Extract the cluster name which is a combination of major clustering tissue, cell type, and specific cluster name.
colnames(cluster_name) <- 'cluster' #Set the column name
cluster_name <- cluster_name %>% separate(cluster, c('clus', 'type', 'name'), sep = '_') #Split the information out to select the cell type annotation. 
neuron.object$cell_type <- cluster_name$type #add the cell type back to the seurat object. 
neuron.object <- NormalizeData(neuron.object, normalization.method = "LogNormalize") #Perfrom normalization on the whole dataset. 
saveRDS(neuron.object, "./MBA_FIXED.rds") #Save this point.
#Subset a sufficent number to capture all cell types but still be managable as a reference set in memory. 
neuron_sub <- subset(x = neuron.object, downsample = 50000)
##NOTE: I am ending the Mouse here because the steps to get to each of the reference types is the same for each species once they get to this point. So I will do the other species to this point then do the steps to generate the references once.

#####HUMAN ANNOTATION REFERENCE PREP


library(Seurat)
library(BPCells) #I had to do a remotes install on this one
library(Azimuth) #I had to do a remotes install on this one 
data.list <- c() #Empty objects to fill
metadata.list <- c()
##NEURON
#Get large dataset as a Seurat Object written to disk
path <- "./c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad"
data <- open_matrix_anndata_hdf5(path) #Open tunnel to counts on disk not into memory.
data <- convert_matrix_type(data, "uint32_t") #Convert counts to appropriate class
write_matrix_dir(mat = data, dir = "./Neuron") #Write to disk in manner readable for Seurat.
mat <- open_matrix_dir(dir = "./Neuron") #Stream the count Matrix
mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human") #Convert to symbol. 
metadata.list[[1]] <- LoadH5ADobs(path = path) #collect metadata from H5 file.
data.list[[1]] <- mat #add to metadata.
metadata <- Reduce(rbind, metadata.list)
neuron.object <- CreateSeuratObject(counts = data.list, meta.data = metadata) #Create Neuron Seurat Object
saveRDS(neuron.object, "NEURON_OG.rds") #Save it
##NONNEURON
#Get large dataset as a Seurat Object written to disk
path <- ".f/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad"
data <- open_matrix_anndata_hdf5(path) #Open tunnel to counts on disk not into memory.
data <- convert_matrix_type(data, "uint32_t") #Convert counts to appropriate class
write_matrix_dir(mat = data, dir = "./NoNNeuron")  #Write to disk in manner readable for Seurat.
mat <- open_matrix_dir(dir = "./NoNNeuron") #Stream the count Matrix
mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human") #Convert to symbol.
metadata.list[[1]] <- LoadH5ADobs(path = path) #collect metadata from H5 file.
data.list[[1]] <- mat #add to metadata.
metadata <- Reduce(rbind, metadata.list)
neuron.object <- CreateSeuratObject(counts = data.list, meta.data = metadata) #Create NON Neuron Seurat Object
saveRDS(neuron.object, "NONNEURON_OG.rds") #Save it
#Clear memory and stream the two objects into memory
non_neuron <- readRDS("./NONNEURON_OG.rds")
neuron <- readRDS("./NEURON_OG.rds")
#combine the objects together and denote Neuron and Non Neuron.
combined_thal <- merge(non_neuron, y = neuron, add.cell.ids = c("NonN", "N"), project = "HBA_REF")
#Subset a sufficent number to capture all cell types but still be managable as a reference set in memory. 
neuron_sub <- subset(combined_thal, downsample = 500)
##NOTE: I am ending the Human here because the steps to get to each of the reference types is the same for each species once they get to this point. So I will do the other species to this point then do the steps to generate the references once.


################################################################################################
#Here we set up cross species use. Only a subset of genes are orthologs. We do assess below that the subset is sufficent for partitioning major cell types with the recognition of subtype or neuron types being clouded.
#####ONE-TO-ONE ortholog set up
library(Seurat)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(data.table)
library(matrixStats)

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


################################################################################################
neuron.object <- readRDS("./MBA_FIXED.rds")
neuron.object <- subset(neuron.object, anatomical_division_label == 'STR')
neuron.object$CT_toplevel <- collapse_to_toplevel(neuron.object$cell_type)

cell_types <- neuron.object@meta.data[, c('cell_type','CT_toplevel','anatomical_division_label')]
colnames(cell_types) <- c('CT', 'CT_toplevel','Tissue')

em <- GetAssayData(neuron.object, layer = 'counts')
emcol <- colnames(em)
emrow <- rownames(em)
em <- as.matrix(em)
em <- as.data.table(em)

#MOUSE
orth_idx <- match(emrow, orth3$Mouse_Ensembl) #Mouse was already gene IDS so matching 121 IDs to each other.
final_gene_ids <- ifelse(!is.na(orth_idx), orth3$Human_Ensembl[orth_idx], emrow) #If no human ID keep mouse ID
rownames(em) <- make.unique(final_gene_ids) #ensure uniqness and add back to count matrix.
fwrite(em, file="./MBA_counts.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write new count matrix
cell_types$orig.ident <- 'MBA'
fwrite(cell_types, file="./MBA_meta.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write metadata
#These would be the steps. but since I already did them I shall just use what was already made.


hnon_neuron <- readRDS("./NONNEURON_OG.rds")
hneuron <- readRDS("./NEURON_OG.rds")
#combine the objects together and denote Neuron and Non Neuron.
hcombined <- merge(hnon_neuron, y = hneuron, add.cell.ids = c("NonN", "N"), project = "HBA_REF")
roi_of_STR <- c(
  "Human CaB", "Human NAC", "Human Pu"
)
roi_of_THAL <- c(
  "Human ANC", "Human CM", "Human CM-Pf", "Human ETH", "Human LG", 
  "Human LP", "Human LP-VPL", "Human MD", "Human MD-Re", "Human MG", 
  "Human Pul", "Human STH", "Human VA", "Human VLN", "Human VPL"  # Thalamus
)


# Subset the Seurat object
hcombined_sub_STR <- subset(hcombined, roi %in% roi_of_STR)
hcombined_sub_STR <- JoinLayers(hcombined_sub_STR)

hcombined_sub_THAL <- subset(hcombined, roi %in% roi_of_THAL)
hcombined_sub_THAL <- JoinLayers(hcombined_sub_THAL)

set.seed(123)  # for reproducibility
hcombined_sub_STR$CT_toplevel <- collapse_to_toplevel(hcombined_sub_STR$supercluster_term)

set.seed(123)  # for reproducibility
hcombined_sub_THAL_sub <- hcombined_sub_THAL[, unlist(lapply(unique(hcombined_sub_THAL$roi), function(r) {
  cells <- WhichCells(hcombined_sub_THAL, expression = roi == r)
  if(length(cells) > 5000) {
    sample(cells, 5000)
  } else {
    cells
  }
}))]
hcombined_sub_THAL_sub$CT_toplevel <- collapse_to_toplevel(hcombined_sub_THAL_sub$supercluster_term)


em <- GetAssayData(hcombined_sub_STR, layer = 'counts')
emcol <- colnames(em)
emrow <- rownames(em)
em <- as.matrix(em)
em <- as.data.table(em)

symbol2ensembl_human <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = human) #Get the gene ID that matches with the gene names
orth_idx <- match(emrow, symbol2ensembl_human$external_gene_name)#Mouse was already gene IDS so matching 121 IDs to each other.
final_gene_ids <- ifelse(!is.na(orth_idx), symbol2ensembl_human$ensembl_gene_id[orth_idx], emrow) #If no human ID keep mouse ID
rownames(em) <- make.unique(final_gene_ids) #ensure uniqness and add back to count matrix.
fwrite(em, file="./HBA_counts_STR.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write new count matrix
cell_types <- hcombined_sub_STR@meta.data[, c('supercluster_term','CT_toplevel', 'roi')]
colnames(cell_types) <- c('CT', 'CT_toplevel','Tissue')
cell_types$orig.ident <- 'HBA'
fwrite(cell_types, file="./HBA_meta_STR.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write metadata


em <- GetAssayData(hcombined_sub_THAL_sub, layer = 'counts')
emcol <- colnames(em)
emrow <- rownames(em)
em <- as.matrix(em)
em <- as.data.table(em)

symbol2ensembl_human <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = human) #Get the gene ID that matches with the gene names
orth_idx <- match(emrow, symbol2ensembl_human$external_gene_name)#Mouse was already gene IDS so matching 121 IDs to each other.
final_gene_ids <- ifelse(!is.na(orth_idx), symbol2ensembl_human$ensembl_gene_id[orth_idx], emrow) #If no human ID keep mouse ID
rownames(em) <- make.unique(final_gene_ids) #ensure uniqness and add back to count matrix.
fwrite(em, file="./HBA_counts_THAL.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write new count matrix
cell_types <- hcombined_sub_THAL_sub@meta.data[, c('supercluster_term','CT_toplevel', 'roi')]
colnames(cell_types) <- c('CT', 'CT_toplevel','Tissue')
cell_types$orig.ident <- 'HBA'
fwrite(cell_types, file="./HBA_meta_THAL.csv", quote=F, row.names = TRUE, col.names = TRUE) #Write metadata

HBAS_Counts <- as.data.frame(fread("./HBA_counts_STR.csv"))
HBAS_meta <- as.data.frame(fread("./HBA_meta_STR.csv"))
rownames(HBAS_Counts) <- HBAS_Counts$V1 #Make gene ID the rownames
HBAS_Counts <- HBAS_Counts[,-1] #Remove gene ID column
rownames(HBAS_meta) <- HBAS_meta[,1] #Make Cell ID the rownames
colnames(HBAS_meta) <- HBAS_meta[1,] #make first row the column names 
HBAS_meta <- HBAS_meta[-1,] #Remove if nessicary
HBAS_meta <- HBAS_meta[,-1] #Remove if nessicary
HBAS_Counts[1:5,1:5] #check orientation and structure
head(HBAS_meta) #check orientation and structure
HBAT_Counts <- as.data.frame(fread("./HBA_counts_THAL.csv"))
HBAT_meta <- as.data.frame(fread("./HBA_meta_THAL.csv"))
rownames(HBAT_Counts) <- HBAT_Counts$V1 #Make gene ID the rownames
HBAT_Counts <- HBAT_Counts[,-1] #Remove gene ID column
rownames(HBAT_meta) <- HBAT_meta[,1] #Make Cell ID the rownames
colnames(HBAT_meta) <- HBAT_meta[1,] #make first row the column names 
HBAT_meta <- HBAT_meta[-1,] #Remove if nessicary
HBAT_meta <- HBAT_meta[,-1] #Remove if nessicary
HBAT_Counts[1:5,1:5] #check orientation and structure
head(HBAT_meta) #check orientation and structure
#MOUSE
MBA_Counts <- as.data.frame(fread("./MBA_counts.csv"))
MBA_meta <- as.data.frame(fread("./MBA_meta.csv"))
rownames(MBA_Counts) <- MBA_Counts$V1 #Make gene ID the rownames
MBA_Counts <- MBA_Counts[,-1] #Remove gene ID column
rownames(MBA_meta) <- MBA_meta[,1] #Make Cell ID the rownames
colnames(MBA_meta) <- MBA_meta[1,] #make first row the column names 
MBA_meta <- MBA_meta[-1,] #Remove if nessicary
MBA_meta <- MBA_meta[,-1] #Remove if nessicary
MBA_Counts[1:5,1:5] #check orientation and structure
head(MBA_meta) #check orientation and structure
#Reduce to ONE-TO-ONE genes shared by all.
common_genes <- Reduce(intersect, list(rownames(HBAS_Counts), rownames(MBA_Counts), rownames(HBAT_Counts)))
#Hounkpe, B. W., Chenou, F., de Lima, F., & De Paula, E. V. (2021). HRT Atlas v1.0 database: Redefining human and mouse housekeeping genes and candidate reference transcripts by mining massive RNA-seq datasets. Nucleic Acids Research, 49(D1), D947–D955. https://doi.org/10.1093/nar/gkaa609
HKG <- fread("./Human_Mouse_Common.csv")

#Match with human symbol
orth_idx <- match(HKG$Human, orth3$Human_Symbol) #Match symbols to get index
HKG$HEID <- orth3$Human_Ensembl[orth_idx] #Use index to find IDs and make a vector of them
HKG <- na.omit(HKG$HEID) #Remove NAs
common_genes <- common_genes[-which(common_genes %in% HKG)] #Subtract from common set of shared genes

#Remove mitochondrial and ribosomal protein genes.
mt_ribo <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "gene_biotype", "external_gene_name"), mart = human)
# Identify mitochondrial genes
mt_genes <- mt_ribo$ensembl_gene_id[mt_ribo$chromosome_name == "MT"]
# Identify ribosomal genes (RPL / RPS)
ribo_genes <- mt_ribo$ensembl_gene_id[grepl("^RPL|^RPS", mt_ribo$external_gene_name)]
bad_genes <- unique(c(mt_genes, ribo_genes))
#final set of genes to use
common_genes <- common_genes[-which(common_genes %in% bad_genes)]

#Subset count matricies down to shared genes.
HBAT_Counts <- HBAT_Counts[rownames(HBAT_Counts) %in% common_genes, ]
HBAS_Counts <- HBAS_Counts[rownames(HBAS_Counts) %in% common_genes, ]
MBA_Counts <- MBA_Counts[rownames(MBA_Counts) %in% common_genes, ]
#Make seurat objects and normalize
HBAT_obj  <- CreateSeuratObject(counts = HBAT_Counts, meta.data = HBAT_meta)
HBAT_obj <- NormalizeData(HBAT_obj, normalization.method = "LogNormalize", scale.factor = 10000)
HBAS_obj  <- CreateSeuratObject(counts = HBAS_Counts, meta.data = HBAS_meta)
HBAS_obj <- NormalizeData(HBAS_obj, normalization.method = "LogNormalize", scale.factor = 10000)
MBA_obj  <- CreateSeuratObject(counts = MBA_Counts, meta.data = MBA_meta)
MBA_obj <- NormalizeData(MBA_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#####NORMALIZE CELL ANNOTATIONS
combined <- merge(HBAT_obj, y = c(MBA_obj, HBAS_obj), add.cell.ids = c("HBAT", "HBAS", "MBA"))
#####Build the References
combined <- JoinLayers(combined) #Combine all counts together
combined <- NormalizeData(combined, normalization.method = "LogNormalize") #Normalization of counts
expr_matrix <- GetAssayData(combined, layer = 'data') #Extract normalized counts
cell_types <- combined@meta.data$CT #Pull out original annotation 
cell_types_top <- combined@meta.data$CT_toplevel #Pull out clean annotation 
#Build SingleR annotation object.
HBA_data <- SummarizedExperiment(
  assays = list(logcounts = as.matrix(expr_matrix)),
  colData = DataFrame(labels = cell_types, top_labels = cell_types_top)
)
#Save it
saveRDS(HBA_data, "./HBAMBA_STRTHL_SINR_REF.rds")
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+



################################################################################################
#####AUTOMATED ANNOTATION FILTER APPLIED
SO <- readRDS("./9823_POST_INT_MAD.rds") #Filtered and integrated Seruat Object.
SO <- FindNeighbors(SO, reduction = "harmony", dims = 1:40)	#neighbor detection
SO <- FindClusters(object = SO, resolution = c(0.8)) #Selected for annotation
Idents(object = SO) <- "SCT_snn_res.0.8" #Set selected group identity

#SO <- combinedMAD
DefaultAssay(SO) <- 'RNA'
SO <- JoinLayers(SO)
SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000) #normalizatioin
SO <- ScaleData(SO) #scaling counts
sinR_in <- LayerData(SO)

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

symbol2ensembl <- getBM(attributes=c("external_gene_name","ensembl_gene_id"), mart=pig) #Get the gene ID that matches with the gene names
gene2ensembl <- symbol2ensembl$ensembl_gene_id[match(rownames(sinR_in), symbol2ensembl$external_gene_name)] #Match pig sybols to gene ID
orth_idx <- match(gene2ensembl, orth3$Pig_Ensembl) #Make an index of the gene IDs that match the human conversion table
final_gene_ids <- ifelse(!is.na(orth_idx), orth3$Human_Ensembl[orth_idx], gene2ensembl) #if missing a human ortholog set as the pig gene ID
rownames(sinR_in) <- make.unique(final_gene_ids) #replace and make unique the rownames of the pig count matrix

#This metadata structre is what is needed for the analyses below. 
PCA <- readRDS("./AnnotationRef/HBAMBA_STRTHL_SINR_REF.rds") #Load all species reference
common_genes <- rownames(sinR_in)[which(rownames(sinR_in) %in% PCA@NAMES)] #identify only genes in the reference
sinR_in <- sinR_in[rownames(sinR_in) %in% common_genes, ] #Keep only genes shared between the target and reference
Idents(object = SO) <- "SCT_snn_res.0.8" #set grouping identity
#Run annotation
ann <- SingleR(test = sinR_in, de.method = "wilcox", clusters = Idents(SO), de.n = 20, ref = PCA, labels = PCA$labels, aggr.ref = TRUE, BNPARAM = RandomParam(deferred = TRUE, rank = 20))
#extract labels and apply them to the cells.
cluster_annotations <- data.frame(cluster = levels(Idents(SO)), cell_type = ann$pruned.labels, score = ann$delta.next)
ident_vec_num <- SO[["SCT_snn_res.0.8"]] # Extract cluster IDs for each cell in the Seurat object
ident_vec_num$cell_type <- cluster_annotations$cell_type[match(ident_vec_num$SCT_snn_res.0.8, cluster_annotations$cluster)] #match cell type with cell type grouping
write.csv(ident_vec_num, file = "./ANN_SinR_Labels.csv", quote=F, row.names=T) #save annotation
SO$CT <- ident_vec_num$cell_type #Add annotation to seurat Object
DimPlot(SO, reduction = "umap", group.by = "CT", raster=FALSE)
png('ANN_SinR_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "CT", label = TRUE, label.size = 6, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + NoLegend()
dev.off()
SO$CT_toplevel <- collapse_to_toplevel(SO$CT)
png('ANN_SinR_toplevel_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "CT_toplevel", label = TRUE, label.size = 6, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + NoLegend()
dev.off()

SO_STR <- subset(SO, tissue == 'Striatum')

# Make sure we are using RNA assay
DefaultAssay(SO) <- "RNA"

# Define Marker Sets
keep_genes <- function(g) intersect(g, rownames(SO))

neuronal_genes <- keep_genes(c("SYT1", "MYT1L", "RBFOX1", "SYNPR", "NSF"))
gaba_genes <- keep_genes(c("GAD2", "DLGAP2", "GRIP1", "GABRG3"))
glut_genes <- keep_genes(c("ARPP21", "NELL2", "SH3GL2", "SV2B"))
msn_genes <- keep_genes(c("GNAL", "ANO3", "ADCY5"))
thal_genes <- keep_genes(c("MLLT11", "TUBB2A", "SCG5", "STMN2"))

oligo_genes <- keep_genes(c("ST18", "CTNNA3", "RNF220", "MBP"))
opc_genes <- keep_genes(c("PTPRZ1", "TNR", "LHFPL3", "EPN2"))
astro_genes <- keep_genes(c("SLC1A2", "GPC5", "PITPNC1", "BMPR1B"))
microglia_genes <- keep_genes(c("APBB1IP", "ARHGAP24", "DOCK8", "CSF1R", "DOCK2"))
endo_genes <- keep_genes(c("FLT1", "VWF", "CLDN5", "EPAS1"))
glut_markers <- c("SLC17A7","SLC17A6")
gaba_markers <- c("GAD1","GAD2","SLC32A1")

# Add Module Scores
SO <- AddModuleScore(SO, list(neuronal_genes), name="NEURON")
SO <- AddModuleScore(SO, list(gaba_genes), name="GABA")
SO <- AddModuleScore(SO, list(glut_genes), name="GLUT")
SO <- AddModuleScore(SO, list(msn_genes), name="MSN")
SO <- AddModuleScore(SO, list(thal_genes), name="THAL")

SO <- AddModuleScore(SO, list(oligo_genes), name="OLIGO")
SO <- AddModuleScore(SO, list(opc_genes), name="OPC")
SO <- AddModuleScore(SO, list(astro_genes), name="ASTRO")
SO <- AddModuleScore(SO, list(microglia_genes), name="MICRO")
SO <- AddModuleScore(SO, list(endo_genes), name="ENDO")

cluster_scores <- SO@meta.data %>%
  mutate(cluster = Idents(SO)) %>%
  group_by(cluster) %>%
  summarise(NEURON = mean(NEURON1),GABA = mean(GABA1), GLUT = mean(GLUT1), MSN = mean(MSN1), THAL = mean(THAL1), OLIGO = mean(OLIGO1), OPC = mean(OPC1), ASTRO = mean(ASTRO1), MICRO = mean(MICRO1), ENDO = mean(ENDO1), n_cells = n())

cluster_scores <- cluster_scores %>%
  mutate(
    neuron_vs_glia =
      NEURON - pmax(OLIGO, OPC, ASTRO, MICRO, ENDO),
    gaba_vs_glut =
      GABA - GLUT,
    msn_vs_thal =
      MSN - THAL
  )


cluster_tissue <- SO@meta.data %>%
  mutate(cluster = Idents(SO)) %>%
  group_by(cluster, tissue) %>%
  summarise(n = n(), .groups="drop") %>%
  tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

cluster_scores <- left_join(cluster_scores, cluster_tissue, by="cluster")

cluster_scores <- cluster_scores %>%
  mutate(
    total_cells = n_cells,
    striatum_prop = Striatum / total_cells,
    thalamus_prop = Thalamus / total_cells
  )

cluster_scores <- cluster_scores %>%
  mutate(
    neuron_confidence = abs(neuron_vs_glia),
    neurotransmitter_confidence = abs(gaba_vs_glut),
    regional_confidence = abs(msn_vs_thal),
    tissue_support =
      pmax(striatum_prop, thalamus_prop)
  )


cluster_scores$Major_class <- "Unassigned"
cluster_scores$Neuron_subtype <- "None"

for (i in 1:nrow(cluster_scores)) {
  row <- cluster_scores[i,]
  glia_max <- max(row$OLIGO, row$OPC, row$ASTRO, row$MICRO, row$ENDO)
  # Neuron detection
  if (row$NEURON > glia_max) {
    cluster_scores$Major_class[i] <- "Neuron"
    # Neurotransmitter
    if (row$gaba_vs_glut > 0) {
      # GABA neuron
      if (row$msn_vs_thal > 0 | row$striatum_prop > 0.6) {
        cluster_scores$Neuron_subtype[i] <- "MSN_GABA"
      } else {
        cluster_scores$Neuron_subtype[i] <- "Other_GABA"
      }
    } else {
      # Glut neuron
      if (row$msn_vs_thal < 0 | row$thalamus_prop > 0.6) {
        cluster_scores$Neuron_subtype[i] <- "Thalamic_Glut"
      } else {
        cluster_scores$Neuron_subtype[i] <- "Other_Glut"
      }
    }
  } else {
    # Glia
    if (row$OLIGO == glia_max) cluster_scores$Major_class[i] <- "Oligodendrocytes"
    if (row$OPC   == glia_max) cluster_scores$Major_class[i] <- "OPC"
    if (row$ASTRO == glia_max) cluster_scores$Major_class[i] <- "Astrocytes"
    if (row$MICRO == glia_max) cluster_scores$Major_class[i] <- "Microglia"
    if (row$ENDO  == glia_max) cluster_scores$Major_class[i] <- "Endothelial"
  }
}

  

cluster_map <- cluster_scores$Major_class
names(cluster_map) <- cluster_scores$cluster

cluster_annotations <- data.frame(cluster = names(cluster_map), cell_type = cluster_map)
ident_vec_num <- SO[["SCT_snn_res.0.8"]] # Extract cluster IDs for each cell in the Seurat object
ident_vec_num$cell_type <- cluster_annotations$cell_type[match(ident_vec_num$SCT_snn_res.0.8, cluster_annotations$cluster)] #match cell type with cell type grouping
write.csv(ident_vec_num, file = "./ANN_Marker.csv", quote=F, row.names=T) #save annotation
SO$Major_class_cluster <- ident_vec_num$cell_type #Add annotation to seurat Object


subtype_map <- cluster_scores$Neuron_subtype
names(subtype_map) <- cluster_scores$cluster

cluster_annotations <- data.frame(cluster = names(subtype_map), cell_type = subtype_map)
ident_vec_num <- SO[["SCT_snn_res.0.8"]] # Extract cluster IDs for each cell in the Seurat object
ident_vec_num$cell_type <- cluster_annotations$cell_type[match(ident_vec_num$SCT_snn_res.0.8, cluster_annotations$cluster)] #match cell type with cell type grouping
write.csv(ident_vec_num, file = "./ANN_Marker_NeurSub.csv", quote=F, row.names=T) #save annotation
SO$Neuron_subtype_cluster <- ident_vec_num$cell_type #Add annotation to seurat Object



prop.table(table(SO$tissue, SO$Major_class_cluster), margin = 1) * 100

prop.table(table(SO$tissue, SO$Neuron_subtype_cluster),1)*100

cluster_scores$dominant_tissue <- ifelse(cluster_scores$striatum_prop > cluster_scores$thalamus_prop, "Striatum", "Thalamus")

cluster_scores$CT <- ifelse(cluster_scores$Neuron_subtype == "None", cluster_scores$Major_class, cluster_scores$Neuron_subtype)
subtype_map <- cluster_scores$CT
names(subtype_map) <- cluster_scores$cluster

cluster_annotations <- data.frame(cluster = names(subtype_map), cell_type = subtype_map)
ident_vec_num <- SO[["SCT_snn_res.0.8"]] # Extract cluster IDs for each cell in the Seurat object
ident_vec_num$cell_type <- cluster_annotations$cell_type[match(ident_vec_num$SCT_snn_res.0.8, cluster_annotations$cluster)] #match cell type with cell type grouping
write.csv(ident_vec_num, file = "./ANN_Marker_CT.csv", quote=F, row.names=T) #save annotation
SO$CT <- ident_vec_num$cell_type #Add annotation to seurat Object

png('ANN_Marker_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "Major_class_cluster", label = TRUE, label.size = 6, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + NoLegend()
dev.off()
png('ANN_Marker_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "CT", label = TRUE, label.size = 6, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2") + NoLegend()
dev.off()


SinR <- read.csv("./ANN_SinR_Labels.csv")
SO$SinR <- SinR$cell
SO$SinR_toplevel <- collapse_to_toplevel(SO$CT)

DefaultAssay(SO) <- "RNA"
Idents(SO) <- "SCT_snn_res.0.8"

# Cluster-level majority vote
meta <- SO@meta.data

cluster_vote <- meta %>%
  group_by(SCT_snn_res.0.8) %>%
  summarise(
    CT_major = names(sort(table(CT), decreasing = TRUE))[1],
    SinR_major = names(sort(table(SinR_toplevel), decreasing = TRUE))[1]
  )

# combine annotations
cluster_vote$Major_class <- ifelse(
  cluster_vote$CT_major == cluster_vote$SinR_major,
  cluster_vote$CT_major,
  cluster_vote$CT_major
)

# attach to object
SO$Major_class_cluster <- cluster_vote$Major_class[match(
  SO$SCT_snn_res.0.8,
  cluster_vote$SCT_snn_res.0.8
)]

#glial identities
glia_types <- c("Astrocytes","Microglia","Oligodendrocytes","OPC","Endothelial")

SO$Major_class_cluster <- ifelse(
  SO$Major_class_cluster %in% glia_types,
  SO$Major_class_cluster,
  "Neuron"
)
# Score neuronal markers
glut_markers <- c("SLC17A7","SLC17A6")
gaba_markers <- c("GAD1","GAD2","SLC32A1")

SO <- AddModuleScore(SO, features=list(glut_markers), name="GlutScore")
SO <- AddModuleScore(SO, features=list(gaba_markers), name="GabaScore")

# Tissue-aware neuron class
gaba_thresh <- 0.1
glut_thresh <- 0.1

SO$Neuron_subtype <- "None"

# Striatal projection neurons
SO$Neuron_subtype[
  SO$Major_class_cluster == "Neuron" &
  SO$GabaScore1 > gaba_thresh &
  SO$tissue == "Striatum"
] <- "Striatal_GABA"

# Thalamic glut neurons
SO$Neuron_subtype[
  SO$Major_class_cluster == "Neuron" &
  SO$GlutScore1 > glut_thresh &
  SO$tissue == "Thalamus"
] <- "Thalamic_Glut"

# Other glut neurons
SO$Neuron_subtype[
  SO$Major_class_cluster == "Neuron" &
  SO$GlutScore1 > glut_thresh &
  SO$Neuron_subtype == "None"
] <- "Other_Glut"

# Other GABA neurons
SO$Neuron_subtype[
  SO$Major_class_cluster == "Neuron" &
  SO$GabaScore1 > gaba_thresh &
  SO$Neuron_subtype == "None"
] <- "Other_GABA"
# check results
SO$Neuron_subtype[SO$Neuron_subtype == "None" &
                  SO$Major_class_cluster == "Neuron"] <- "Unclassified_Neuron"
# Final label
SO$Final_celltype <- ifelse(
  SO$Major_class_cluster == "Neuron",
  SO$Neuron_subtype,
  SO$Major_class_cluster
)

# check results
table(SO$tissue, SO$Final_celltype)
DimPlot(SO, group.by="Final_celltype", label=TRUE)
png('ANN_Final_celltype_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "Final_celltype", label = TRUE, label.size = 6, repel = TRUE) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()
png('ANN_Final_celltype_Tiss_UMAP.png', width = 1200, height = 700)
DimPlot(SO, reduction = "umap", group.by = "Final_celltype", split.by = 'tissue', label = F) + theme(text = element_text(size = 16)) + labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
dev.off()

OutLabels <- data.frame(Cluster = SO$SCT_snn_res.0.8, Tissue = SO$tissue, SinR = SO$SinR, SinR_toplevel = SO$SinR_toplevel, Major_class_cluster = SO$Major_class_cluster, Final_celltype = SO$Final_celltype)
write.csv(OutLabels, 'Full_Labels.csv', quote=F, row.names=T)


gaba_genes <- keep_genes(c("GAD2", "DLGAP2", "GRIP1", "GABRG3", "GAD1","SLC32A1"))
glut_genes <- keep_genes(c("ARPP21", "NELL2", "SH3GL2", "SV2B", "SLC17A7","SLC17A6"))

oligo_genes <- keep_genes(c("ST18", "CTNNA3", "RNF220", "MBP"))
opc_genes <- keep_genes(c("PTPRZ1", "TNR", "LHFPL3", "EPN2"))
astro_genes <- keep_genes(c("SLC1A2", "GPC5", "PITPNC1", "BMPR1B"))
microglia_genes <- keep_genes(c("APBB1IP", "ARHGAP24", "DOCK8", "CSF1R", "DOCK2"))
endo_genes <- keep_genes(c("FLT1", "VWF", "CLDN5", "EPAS1"))

SO$CT <- factor(SO$CT, levels = c("Striatal_GABA", "Other_GABA", "Thalamic_Glut", "Other_Glut", "Neuron_Unclassified", "Oligodendrocytes", "OPC", "Astrocytes", "Microglia", "Endothelial"))

markers <- c("GAD1", "GAD2", "DLGAP2", "SLC32A1", "NELL2", "SV2B", "SLC17A7", "SLC17A6", "SOX2", "NES", "DCX", "NEUROD1", "ST18", "CTNNA3", "RNF220", "MBP", "PTPRZ1", "TNR", "LHFPL3", "EPN2", "SLC1A2", "GPC5", "PITPNC1", "BMPR1B", "ARHGAP24", "DOCK8", "CSF1R", "DOCK2", "FLT1", "VWF", "CLDN5", "EPAS1", 'GFP')
png('ANN_celltype_Dot.png', width = 1200, height = 700)
DotPlot(SO, features = markers, group.by = 'CT') + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), axis.text.y = element_text(size = 18), text = element_text(size = 18))
dev.off()
