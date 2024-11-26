checkMart <- function(species, updateMart = FALSE){
  if(species == "Human"){
    speciesDataset <- "hsapiens_gene_ensembl"
    martFile <- "./data/humanMart.rds"
  }
  else if(species == "Mouse"){
    speciesDataset <- "mmusculus_gene_ensembl"
    martFile <- "./data/mouseMart.rds"
  }
  if (file.exists(martFile) & updateMart == FALSE) {
    mart <- readRDS(martFile)
  } 
  else if (!file.exists(martFile) | updateMart == TRUE){
    mart <- useMart("ensembl", dataset = speciesDataset)
    saveRDS(mart, martFile)
  }
  return(mart)
}

loadAndPreprocess <- function(folderCellRangerOut, gene_expression_cutoff, spot_gene_cutoff, species){
  withProgress(message = "Loading data...", value = 0, {
    incProgress(0.1, detail = "Preparing data...")
    seuratObj <- Load10X_Spatial(folderCellRangerOut)
    
    incProgress(0.1, detail = "Filtering genes...")
    # Filter genes based on minimum expression in % of cells
    percent_expressed <- rowSums(seuratObj[["Spatial"]]$counts > 0) / ncol(seuratObj) * 100
    genes_to_keep <- names(percent_expressed[percent_expressed >= gene_expression_cutoff])
    filtered_genes <- setdiff(rownames(seuratObj), genes_to_keep)
    seuratObj <- subset(seuratObj, features = genes_to_keep)
    
    incProgress(0.1, detail = "Filtering spots...")
    # Filter spots based on minimum number of genes expressed per spot
    expressed_genes_per_spot <- colSums(seuratObj[["Spatial"]]$counts > 0)
    spots_to_keep <- names(expressed_genes_per_spot[expressed_genes_per_spot >= spot_gene_cutoff])
    filtered_spots <- setdiff(colnames(seuratObj), spots_to_keep)
    seuratObj <- subset(seuratObj, cells = spots_to_keep)
    
    incProgress(0.4, detail = "Normalizing and scaling the data...")
    # Normalize and Scale the data
    seuratObj <- NormalizeData(seuratObj, normalization.method = "RC", scale.factor = 1e6)
    seuratObj <- ScaleData(seuratObj)
    
    incProgress(0.1, detail = "Loading appropriate mart...")
    mart <- checkMart(species)

    # Return the filtered seurat object and the counts of filtered genes and spots
    return(list(seuratObj = seuratObj, filtered_genes = length(filtered_genes), 
                filtered_spots = length(filtered_spots), mart = mart))
  })
}

loadClusterMat <- function(filenameCluster, seuratObj) {
  clusterMat <- read.csv2(filenameCluster, sep = ",", row.names = 1)
  
  barcodes <- colnames(seuratObj[["Spatial"]]$data)
  
  common_barcodes <- intersect(barcodes, rownames(clusterMat))
  clusterMat <- clusterMat[common_barcodes, , drop = FALSE]
  
  matched_indices <- match(barcodes, rownames(clusterMat))
  
  valid_indices <- !is.na(matched_indices)
  clusterMat <- clusterMat[matched_indices[valid_indices], , drop = FALSE]
  
  return(clusterMat)
}


prepare_gene_data <- function(gene, data_loaded) {
  countMatrix <- data_loaded$seuratObj[["Spatial"]]$data
  
  gene_data <- data.frame(
    Expression = countMatrix[gene, ], 
    Barcode = colnames(countMatrix),
    Cluster = data_loaded$clusterMat[,1]
  )
  
  return(gene_data)
}


create_violin_plot <- function(gene_data, gene) {
  ggplot(gene_data, aes(x = factor(Cluster), y = Expression)) +
    geom_violin(aes(fill = factor(Cluster)), trim = TRUE) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "gray") +
    theme_minimal() +
    labs(title = paste("Violin Plot for", gene), x = "Cluster", y = "Expression Level")
}

create_beeswarm_plot <- function(gene_data, gene) {
  ggplot(gene_data, aes(x = factor(Cluster), y = Expression)) +
    geom_quasirandom(aes(color = factor(Cluster)), size = 0.8, stroke = 0.3) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "gray") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_minimal() +
    labs(title = paste("Beeswarm Plot for", gene), x = "Cluster", y = "Expression Level")
}

create_plot_with_stats <- function(plot_func, gene_data, gene, comparisons, display_pval) {
  p <- plot_func(gene_data, gene)
  
  if (length(comparisons) > 0) {
    # Toggle label based on display_pval value
    label_format <- if (display_pval) "p.value" else "p.signif"
    
    p <- p + stat_compare_means(comparisons = comparisons, 
                                method = "wilcox.test",
                                label = label_format)
  }
  
  p <- p + stat_compare_means(
    method = "kruskal.test",
    label.x = 0.5,  
    label.y = Inf,
    vjust = 1.2,
    hjust = 0
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

format_pval <- function(pval, threshold = 1e-6) {
  ifelse(pval < threshold, format(pval, scientific = TRUE, digits = 3), round(pval, 3))
}

requestGeneTable <- function(mart, species){
  full_gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                         mart = mart)
  saveRDS(full_gene_map, paste0("./data/", tolower(species), "GeneTable.rds"))
}

convertGeneMap <- function(genes, mart, species){
  genes <- rownames(genes)
  
  tryCatch({
    full_gene_map <- readRDS(paste0("./data/", tolower(species), "GeneTable.rds"))
  }, error = function(e) {
    message("Gene table not found. Creating the table...")
    requestGeneTable(mart, species)
    full_gene_map <- readRDS(paste0("./data/", tolower(species), "GeneTable.rds"))
  })
  
  gene_map <- full_gene_map |>
    dplyr::filter(hgnc_symbol %in% genes) |>
    dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
  
  return(gene_map)
}

runPathwayAnalysis <- function(genes, method, database, species, mart) {
  incProgress(0.1, detail = "Running Pathway Analysis")
  
  # Convert genes to Ensembl format
  gene_map <- convertGeneMap(genes, mart, species)
  genes <- merge(genes, gene_map, by.x = "row.names", by.y = "hgnc_symbol")

  
  # Set species database for GO terms
  go_species <- ifelse(species == "Human", "org.Hs.eg.db", "org.Mm.eg.db")
  if(database == "HALLMARK"){
    hallmark_gene_sets <- msigdbr(species = tolower(species), category = "H")
    hallmark_gene_list <- hallmark_gene_sets |>
      dplyr::select(gs_name, entrez_gene)
  }
  
  # ORA method
  if (method == "ORA") {
    gene_list <- genes |>
      filter(p_val_adj <= 0.05)
    print(gene_list)
    
    # Enrichment based on the selected database
    result <- switch(database,
                     GO = enrichGO(gene = gene_list$ensembl_gene_id,
                                   OrgDb = go_species,
                                   keyType = "ENSEMBL",
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.2),
                     KEGG = enrichKEGG(gene = gene_map$entrezgene_id,
                                       organism = ifelse(species == "Human", "hsa", "mmu"),
                                       pvalueCutoff = 0.05),
                     HALLMARK = enricher(gene = gene_map$entrezgene_id,
                                         TERM2GENE = hallmark_gene_list,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05),
                     stop("Unsupported database"))
    
    # FGSEA method
  } else if (method == "FGSEA") {
    genes_sorted <- genes |>
      dplyr::arrange(desc(avg_log2FC))
    genes_sorted$hgnc_symbol <- rownames(genes_sorted)
    genes_sorted <- merge(genes_sorted, gene_map, by = "hgnc_symbol", by.y = "entrezgene_id")
    genes_sorted <- genes_sorted[!is.na(genes_sorted$entrezgene_id),]

    ranks <- as.numeric(genes_sorted$avg_log2FC)
    names(ranks) <- genes_sorted$entrezgene_id
    ranks <- sort(ranks, decreasing = TRUE)
    ranks <- ranks[!duplicated(names(ranks))]
    
    result <- switch(database,
                         GO = gseGO(geneList = ranks,
                           OrgDb = go_species,
                           ont = "BP",
                           pvalueCutoff = 1),
                         KEGG = gseKEGG(geneList = ranks,
                                        keyType = "ncbi-geneid",
                                        organism = ifelse(species == "Human", "hsa", "mmu"),
                                        pvalueCutoff = 1),
                         HALLMARK = GSEA(ranks, 
                                         TERM2GENE = hallmark_gene_list,
                                         pvalueCutoff = 1),
                         stop("Unsupported database"))
    
  } else {
    stop("Unsupported method")
  }
  
  return(result)
}
