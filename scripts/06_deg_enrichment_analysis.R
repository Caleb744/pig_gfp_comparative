setwd("./results/GFP")
#Performing enrichment analysis. We did get genes that were differentially expressed but this is expected. But does putting some expectation on happenchance and difference (FC) lead to a list of genes that is meaningful?
library(fgsea)
library(org.Ss.eg.db)  # Pig annotation; use org.Mm.eg.db for mouse, org.Hs.eg.db for human
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(biomaRt)
#####Function to convert genes to the right ID for use by enrichment tools.
convert_ids <- function(genes) {
  cat("[ID] Converting genes...\n")
  genes <- genes[!genes %in% c(gene_target, "(Intercept)")] #We expect GFP and the intercept to be differnet so remove them
  ens_ids <- grep("^ENS", genes, value = TRUE) #Find Gene ID genes
  sym_ids <- setdiff(genes, ens_ids) #Find Gene Symbol genes
  df <- data.frame(original = character(), entrez = character()) #Empty dataframe for the conversion
  if (length(ens_ids) > 0) {
    mart <- biomaRt::useMart("ensembl", dataset = "sscrofa_gene_ensembl") #Select bio mart with relevant information
    ens_map <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = "ensembl_gene_id", values = ens_ids, mart = mart) #Pull out needed ID information
    df <- rbind(df, data.frame(original = ens_map$ensembl_gene_id, entrez = as.character(ens_map$entrezgene_id))) #Bind IDs together
  }
  if (length(sym_ids) > 0) {
    sym_map <- bitr(sym_ids, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Ss.eg.db) #Convert the symbol to an ID
    df <- rbind(df, data.frame(original = sym_map$SYMBOL, entrez = sym_map$ENTREZID)) #Bind IDs together
  }
  df <- df[!is.na(df$entrez), ] #Output only those with ID readable by the tools
  df
}

#####GENE LIST PRIOTIZATION
for (s in c("GFP", "TUBB", "NUCKS1")) {
	gene_target <- s
	gene_pos <- paste0(gene_target, "+")
	gene_neg <- paste0(gene_target, "-")
	DEG <- read.csv(paste0("DEGs_", gene_target, "_MAST.csv")) #Read results in
	DEG$p_val_adj[DEG$p_val_adj == 0] <- 1e-320  # Prevent Inf before volcano
	hist(DEG$avg_log2FC, breaks = 100) #Histogram of results
	DEG$Significant <- DEG$p_val_adj < 0.001 & abs(DEG$avg_log2FC) > 0.5 #More libral filtering easing difference of expression
	write.csv(DEG[DEG$Significant, ], paste0("DEGs_", gene_target, "_SIG_GENES.csv"), quote = FALSE, row.names = FALSE) #Write signficant genes
	DEG <- DEG[DEG$X != gene_target, ]
	# Compute ranking metric
	DEG$rank_stat <- sign(DEG$avg_log2FC) * (-log10(DEG$p_val_adj))

	#Volcano plot of all dataset results
	p_volcano <- ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significant)) + geom_point(alpha = 0.7) + scale_color_manual(values = c("grey", "red")) + geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + labs(x = print0("log2 Fold Change (", gene_target, "+ vs ", gene_target ,"−)"), y = "-log10(FDR)", color = "Status") + theme_minimal()
	ggsave(paste0("DEGs_", gene_target, "_volcano.png"), p_volcano, width = 6, height = 5)


	# Rank and convert IDs
	gene_ranks <- DEG[!is.na(DEG$rank_stat), c("X", "rank_stat", "Significant")] #Pull gene, ranking variable and signficance designation into a table
	gene_ranks <- gene_ranks[order(gene_ranks$rank_stat, decreasing = TRUE), ] #Order based on ranking variable
	id_map <- convert_ids(gene_ranks$X) #Use the new function to convert the gene names/IDs to ENTREZID
	merged <- merge(gene_ranks, id_map, by.x = "X", by.y = "original") #Merge ranked gene table with ID table
	merged <- merged[order(merged$rank_stat, decreasing = TRUE), ] #Just to be sure rerank
	ranks <- merged$rank_stat #Pull out the ranking variable
	names(ranks) <- merged$entrez #Assign the gene to the ranking variable

	#####ENRICHMENT ANLAYSIS on whole dataset
	gsea_res <- tryCatch({gseKEGG(geneList = ranks, organism = "ssc", pvalueCutoff = 0.05)}, error = function(e) NULL) #Gene Set Enrichment analysis with ranked list of all tested genes
	write.csv(gsea_res@result, paste0("DEGs_", gene_target, "_gsea.csv"), row.names = FALSE)

	eg <- merged[which(merged$Significant == "TRUE"),] #Subset to "signficant genes"
	go <- enrichGO(gene = eg$entrez, OrgDb = org.Ss.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05) #Enrichment analysis with Gene Ontology Terms
	kegg <- enrichKEGG(gene = eg$entrez, organism = "ssc", keyType = "kegg", pAdjustMethod = "BH", qvalueCutoff = 0.05) #Enrichment analysis with KEGG pathway Terms
	res <- list() #Empty results obejct
	res$go <- if (!is.null(go)) as.data.frame(go) else NULL #Add results to results object
	res$kegg <- if (!is.null(kegg)) as.data.frame(kegg) else NULL #Add results to results object
	write.csv(res$kegg, file = paste0("DEGs_", gene_target, "_kegg.csv"), quote=F, row.names=F)
	write.csv(res$go, file = paste0("DEGs_", gene_target, "_go.csv"), quote=F, row.names=F)

	#####CELL TYPE ENRICHMENT ANLAYSIS
	DEG <- read.csv(paste0("DEGs_CT_MAST_", gene_target, ".csv"))
	if (!dir.exists(paste0("GSEA_",gene_target))) dir.create(paste0("GSEA_",gene_target)) #Create output directory
	DEG <- DEG[DEG$gene != gene_target, ] #Remove GFP from results
	celltypes <- unique(DEG$celltype) #Get list of cell types (cell groupings)
	#Loop to perform signficance and enrichment analyses
	for (ct in celltypes) {
	  message("Processing ", ct, " ...")
	  df <- DEG[DEG$celltype == ct, ] #Subset to only the target cell type
	  if (nrow(df) < 50) next #If less than 50 gene skip
	  df$Significant <- df$p_val_adj < 0.001 & abs(df$avg_log2FC) > 0.5 #Signficance designation. 
	  plot_name <- gsub(" / ", "_", ct) #Cleaning up cell type name
	  plot_name <- gsub("/", "_", plot_name) #Cleaning up cell type name
	  df$p_val_adj[df$p_val_adj == 0] <- 1e-320  # Prevent Inf before volcano
	  write.csv(df[df$Significant, ], paste0("GSEA_",gene_target,"/", plot_name, "_SIG_GENES.csv"), quote = FALSE, row.names = FALSE) #Write signficant genes
	  # Volcano
	  p_volcano <- ggplot(df, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(color = Significant), alpha = 0.7) + scale_color_manual(values = c("grey", "red"))
	  ggsave(paste0("GSEA_",gene_target,"/", plot_name, "_volcano.png"), p_volcano)
	  # Compute ranking metric
	  df$rank_stat <- sign(df$avg_log2FC) * (-log10(df$p_val_adj))
	  if (all(df$rank_stat == 0, na.rm = T)) {
		message("Skipping ", ct, " — no signal for ranking / GSEA_", gene_target, " .")
		next
	  }
	  gene_ranks <- df[order(df$rank_stat, decreasing = TRUE), c("gene", "rank_stat", "Significant")] #Extract and rank the following colums by the ranking metric
	  id_map <- convert_ids(gene_ranks$gene) #Convert the genes to the right ID
	  merged <- merge(gene_ranks, id_map, by.x = "gene", by.y = "original") #Merge ranked table with id table
	  if (nrow(merged) < 20) next
	  ranks <- merged$rank_stat #Pull the ranking stat
	  names(ranks) <- merged$entrez #And name them with the gene they correlate with
	  # GSEA
	  gsea_res <- tryCatch(gseKEGG(ranks, organism = "ssc", pvalueCutoff = 0.05), error = function(e) NULL)
	  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
		write.csv(gsea_res@result, paste0("GSEA_",gene_target,"/", plot_name, "_GSEA_results.csv"), row.names = FALSE)
	  }
	  # For GO/KEGG over-representation: only significant gene Entrez
	  sig_entrez <- merged$entrez[merged$Significant]
	  if (length(sig_entrez) > 0) {
		go <- tryCatch(enrichGO(sig_entrez, OrgDb = org.Ss.eg.db, keyType = "ENTREZID"), error = function(e) NULL)
		kegg <- tryCatch(enrichKEGG(sig_entrez, organism = "ssc", pAdjustMethod = "BH"), error = function(e) NULL)
		if (!is.null(go)) write.csv(as.data.frame(go), paste0("GSEA_",gene_target,"/", plot_name, "_go_results.csv"), row.names = FALSE)
		if (!is.null(kegg)) write.csv(as.data.frame(kegg), paste0("GSEA_",gene_target,"/", plot_name, "_kegg_results.csv"), row.names = FALSE)
	  } else {
		message("No significant genes for ", ct, " — skipping ORA (GO/KEGG).")
	  }
	}
}

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
