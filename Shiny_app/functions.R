loadAndPreprocess <- function(folderCellRangerOut){
  seuratObj <- Load10X_Spatial(folderCellRangerOut)
  percent_expressed <- rowSums(seuratObj[["Spatial"]]$counts > 0) / ncol(seuratObj) * 100
  # Keep genes expressed in at least 1% of cells
  genes_to_keep <- names(percent_expressed[percent_expressed >= 1])
  seuratObj <- subset(seuratObj, features = genes_to_keep)
  
  seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize")
  
  return(seuratObj)
}

loadClusterMat <- function(filenameCluster, seuratObj) {
  clusterMat <- read.csv2(filenameCluster, sep = ",", row.names = 1)
  
  barcodes <- colnames(seuratObj[["Spatial"]]$data)
  
  common_barcodes <- intersect(barcodes, rownames(clusterMat))
  clusterMat <- clusterMat[common_barcodes, , drop = FALSE]
  clusterMat <- clusterMat[match(barcodes, rownames(clusterMat)), , drop = FALSE]
  
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
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.15)))
  
  return(p)
}