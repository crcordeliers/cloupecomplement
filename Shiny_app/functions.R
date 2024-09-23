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
  clusterMat <- read.csv2(filenameCluster, sep = ",")
  
  # Ensure alignment between barcodes in clusterMat and seuratObj
  barcodes <- colnames(seuratObj[["Spatial"]]$data)
  clusterMat <- clusterMat[clusterMat$Barcode %in% barcodes, ]
  
  # Print the number of rows in clusterMat after filtering
  print(paste("Number of rows in clusterMat after filtering:", nrow(clusterMat)))
  
  return(clusterMat)
}


prepare_gene_data <- function(gene, data_loaded) {
  countMatrix <- data_loaded$seuratObj[["Spatial"]]$data
  
  gene_data <- data.frame(
    Expression = countMatrix[gene, ], 
    Barcode = colnames(countMatrix),
    Cluster = data_loaded$clusterMat[,2]
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