loadAndPreprocess <- function(folderCellRangerOut, gene_expression_cutoff, spot_gene_cutoff){
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
    
    incProgress(0.3, detail = "Normalizing the data...")
    # Normalize and Scale the data
    seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize")
    
    incProgress(0.2, detail = "Scaling the data...")
    seuratObj <- ScaleData(seuratObj)
    
    # Return the filtered seurat object and the counts of filtered genes and spots
    list(seuratObj = seuratObj, filtered_genes = length(filtered_genes), filtered_spots = length(filtered_spots))
  })
}

loadClusterMat <- function(filenameCluster, seuratObj) {
  clusterMat <- read.csv2(filenameCluster, sep = ",", row.names = 1)
  
  barcodes <- colnames(seuratObj[["Spatial"]]$data)
  
  common_barcodes <- intersect(barcodes, rownames(clusterMat))
  clusterMat <- clusterMat[common_barcodes, , drop = FALSE]
  clusterMat <- clusterMat[match(barcodes, rownames(clusterMat)), drop = FALSE]
  
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

convert_to_ensembl <- function(genes) {
  incProgress(0.2, detail = "Query of Ensembl IDs")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  incProgress(0.2, detail = "Mapping Ensembl IDs to Symbol")
  gene_map <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = genes, mart = mart)
  gene_map <- gene_map %>% distinct(hgnc_symbol, .keep_all = TRUE)
  
  return(gene_map)
}

runPathwayAnalysis <- function(genes, method = "clusterProfiler") {
  incProgress(0.2, detail = "Running Pathway Analysis")
  
  ensemblGenes <- convert_to_ensembl(genes)
  
  if (method == "clusterProfiler") {
    incProgress(0.2, detail = "Enrichment analysis using ClusterProfiler")
    result <- enrichGO(gene = ensemblGenes$ensembl_gene_id, OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
  } else if (method == "fgsea") {
    incProgress(0.2, detail = "Enrichment analysis using fgsea")
    pathways <- fgsea::examplePathways
    ranks <- stats::rnorm(length(ensemblGenes))
    names(ranks) <- ensemblGenes
    result <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  }
  
  incProgress(0.2, detail = "Rendering plots & data table")
  
  return(result)
}
