checkMart <- function(species, updateMart = FALSE){
  if (species == "Human") {
    speciesDataset <- "hsapiens_gene_ensembl"
    martFile <- "./data/humanMart.rds"
  }
  else if (species == "Mouse") {
    speciesDataset <- "mmusculus_gene_ensembl"
    martFile <- "./data/mouseMart.rds"
  }
  if (file.exists(martFile) & updateMart == FALSE) {
    mart <- readRDS(martFile)
  } 
  else if (!file.exists(martFile) | updateMart == TRUE) {
    mart <- useMart("ensembl", dataset = speciesDataset)
    saveRDS(mart, martFile)
  }
  return(mart)
}

loadAndPreprocess <- function(h5FilePath, gene_expression_cutoff, spot_gene_cutoff, species, normalisation_method){
  withProgress(message = "Loading data...", value = 0, {
    incProgress(0.1, detail = "Preparing data...")
    data <- Read10X_h5(h5FilePath)
    seuratObj <- CreateSeuratObject(counts = data)
    seuratObj$sample <- sub(".*-", "", colnames(seuratObj))
    
    incProgress(0.1, detail = "Filtering genes...")
    # Filter genes based on minimum expression in % of cells
    percent_expressed <- rowSums(GetAssayData(seuratObj, layer = "counts") > 0) / ncol(seuratObj) * 100
    genes_to_keep <- names(percent_expressed[percent_expressed >= gene_expression_cutoff])
    filtered_genes <- setdiff(rownames(seuratObj), genes_to_keep)
    seuratObj <- subset(seuratObj, features = genes_to_keep)
    
    incProgress(0.1, detail = "Filtering spots...")
    # Filter spots based on minimum number of genes expressed per spot
    expressed_genes_per_spot <- colSums(GetAssayData(seuratObj, layer = "counts") > 0)
    spots_to_keep <- names(expressed_genes_per_spot[expressed_genes_per_spot >= spot_gene_cutoff])
    filtered_spots <- setdiff(colnames(seuratObj), spots_to_keep)
    seuratObj <- subset(seuratObj, cells = spots_to_keep)
    
    incProgress(0.4, detail = "Normalizing and scaling the data...")
    # Normalize and Scale the data
    if (normalisation_method == "LogNormalize") {
      seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize")
      seuratObj <- ScaleData(seuratObj)
      seuratObj <- FindVariableFeatures(seuratObj)
      if (length(unique(seuratObj$sample)) > 1) {
        incProgress(0.2, detail = "Correcting batch effect...")
        seuratObj <- RunPCA(seuratObj)
        seuratObj <- RunHarmony(seuratObj, group.by.vars = "sample")
      }
    } else if (normalisation_method == "SCTransform") {
      options(future.globals.maxSize = 2 * 1024^3)
      if (length(unique(seuratObj$sample)) > 1) {
        seuratObj <- SCTransform(seuratObj, vars.to.regress = "sample")
        seuratObj <- RunPCA(seuratObj)
        seuratObj <- RunHarmony(seuratObj, group.by.vars = "sample")
      } else {
        seuratObj <- SCTransform(seuratObj)
      }
    }
    
    incProgress(0.2, detail = "Loading appropriate mart...")
    mart <- checkMart(species)
    
    # Return the filtered seurat object and the counts of filtered genes and spots
    return(list(seuratObj = seuratObj, filtered_genes = length(filtered_genes), 
                filtered_spots = length(filtered_spots), mart = mart))
  })
}

loadClusterMat <- function(filenameCluster, seuratObj) {
  clusterMat <- read.csv2(filenameCluster, sep = ",", row.names = 1)
  
  barcodes <- colnames(GetAssayData(seuratObj, layer = "data"))
  
  common_barcodes <- intersect(barcodes, rownames(clusterMat))
  clusterMat <- clusterMat[common_barcodes, , drop = FALSE]
  
  matched_indices <- match(barcodes, rownames(clusterMat))
  
  valid_indices <- !is.na(matched_indices)
  clusterMat <- clusterMat[matched_indices[valid_indices], , drop = FALSE]
  
  # Order by cluster number to avoid lexicographic order
  clusterMat[,1] <- factor(clusterMat[,1], levels = mixedsort(unique(clusterMat[,1])))
  
  return(clusterMat)
}


prepare_gene_data <- function(gene, data_loaded) {
  countMatrix <- GetAssayData(data_loaded$seuratObj, layer = "data")
  
  barcodes <- colnames(countMatrix)
  
  gene_data <- data.frame(
    Expression = countMatrix[gene, ], 
    Barcode = colnames(countMatrix),
    Cluster = data_loaded$clusterMat[barcodes,1]
  )
  
  return(gene_data)
}


create_violin_plot <- function(gene_data, gene) {
  ggplot(gene_data, aes(x = factor(Cluster), y = Expression)) +
    geom_violin(aes(fill = factor(Cluster)), trim = TRUE) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "gray") +
    theme_minimal() +
    labs(title = paste("Violin Plot for", gene), x = NULL, y = "Expression Level", fill = "Cluster")
}

create_beeswarm_plot <- function(gene_data, gene) {
  ggplot(gene_data, aes(x = factor(Cluster), y = Expression)) +
    geom_quasirandom(aes(color = factor(Cluster)), size = 0.8, stroke = 0.3) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "gray") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_minimal() +
    labs(title = paste("Beeswarm Plot for", gene), x = NULL, y = "Expression Level", color = "Cluster")
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
  if (tolower(species) == "mouse") {
    attrs <- c("ensembl_gene_id", "mgi_symbol", "entrezgene_id")
  } else {
    attrs <- c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id")
  }
  
  full_gene_map <- getBM(attributes = attrs, mart = mart)
  saveRDS(full_gene_map, paste0("./data/", tolower(species), "GeneTable.rds"))
}

convertGeneMap <- function(genes, mart, species){
  genes <- rownames(genes)
  
  symbol_col <- switch(tolower(species),
                       "human" = "hgnc_symbol",
                       "mouse" = "mgi_symbol",
                       stop("Unsupported species"))
  
  tryCatch({
    full_gene_map <- readRDS(paste0("./data/", tolower(species), "GeneTable.rds"))
  }, error = function(e) {
    message("Gene table not found. Creating the table...")
    requestGeneTable(mart, species)
    full_gene_map <- readRDS(paste0("./data/", tolower(species), "GeneTable.rds"))
  })
  
  gene_map <- full_gene_map |>
    dplyr::filter(.data[[symbol_col]] %in% genes) |>
    dplyr::distinct(.data[[symbol_col]], .keep_all = TRUE)
  
  return(gene_map)
}

runPathwayAnalysis <- function(genes, method, database, species, mart) {
  incProgress(0.1, detail = "Running Pathway Analysis")

  # Set gene symbol column based on species
  gene_symbol_col <- ifelse(tolower(species) == "mouse", "mgi_symbol", "hgnc_symbol")

  # Convert genes to Ensembl/Entrez
  gene_map <- convertGeneMap(genes, mart, species)
  genes[[gene_symbol_col]] <- rownames(genes)
  genes <- merge(genes, gene_map, by.x = gene_symbol_col, by.y = gene_symbol_col)

  # Set OrgDb for GO terms
  go_species <- switch(tolower(species),
                       human = "org.Hs.eg.db",
                       mouse = "org.Mm.eg.db",
                       stop("Unsupported species for GO analysis"))

  # Prepare hallmark gene sets if needed
  if (database == "HALLMARK") {
    hallmark_gene_sets <- msigdbr(species = tolower(species), category = "H")
    hallmark_gene_list <- hallmark_gene_sets |>
      dplyr::select(gs_name, entrez_gene)
  }

  # ORA method
  if (method == "ORA") {
    gene_list <- genes |>
      dplyr::filter(p_val_adj <= 0.05) |>
      dplyr::filter(avg_log2FC >= 0)

    result <- switch(database,
                     GO = enrichGO(gene = gene_list$ensembl_gene_id,
                                   OrgDb = get(go_species),
                                   keyType = "ENSEMBL",
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.2),
                     KEGG = enrichKEGG(gene = gene_list$entrezgene_id,
                                       organism = ifelse(tolower(species) == "human", "hsa", "mmu"),
                                       pvalueCutoff = 0.05),
                     HALLMARK = enricher(gene = gene_list$entrezgene_id,
                                         TERM2GENE = hallmark_gene_list,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05),
                     stop("Unsupported database"))

    # FGSEA method
  } else if (method == "FGSEA") {
    genes_sorted <- genes |>
      dplyr::arrange(desc(avg_log2FC)) |>
      dplyr::filter(!is.na(entrezgene_id))

    ranks <- as.numeric(genes_sorted$avg_log2FC)
    names(ranks) <- genes_sorted$entrezgene_id
    ranks <- sort(ranks, decreasing = TRUE)
    ranks <- ranks[!duplicated(names(ranks))]

    result <- switch(database,
                     GO = gseGO(geneList = ranks,
                                OrgDb = get(go_species),
                                ont = "BP",
                                keyType = "ENTREZID",
                                pvalueCutoff = 1),
                     KEGG = gseKEGG(geneList = ranks,
                                    keyType = "ncbi-geneid",
                                    organism = ifelse(tolower(species) == "human", "hsa", "mmu"),
                                    pvalueCutoff = 1),
                     HALLMARK = GSEA(ranks,
                                     TERM2GENE = hallmark_gene_list,
                                     pvalueCutoff = 1),
                     stop("Unsupported database"))

  } else {
    stop("Unsupported method")
  }

  # Convert geneID from Entrez/ENSEMBL to common gene symbols
  if (!is.null(result) && nrow(result@result) > 0) {
    # Create case-insensitive lookup for gene symbols
    gene_map_lower <- gene_map
    gene_map_lower[[gene_symbol_col]] <- tolower(gene_map[[gene_symbol_col]])

    result@result$geneID <- sapply(strsplit(as.character(result@result$geneID), "/"), function(gene_ids) {
      # Convert each ID to gene symbol
      converted <- sapply(gene_ids, function(id) {
        original_id <- id  # Keep original ID in case we can't convert

        # Try exact match first (Entrez or ENSEMBL)
        match_row <- gene_map[gene_map$entrezgene_id == id | gene_map$ensembl_gene_id == id, ]

        if (nrow(match_row) > 0 && !is.na(match_row[[gene_symbol_col]][1]) && match_row[[gene_symbol_col]][1] != "") {
          return(as.character(match_row[[gene_symbol_col]][1]))
        }

        # Try case-insensitive match on gene symbol itself (in case ID is already a symbol)
        id_lower <- tolower(id)
        match_row_lower <- gene_map_lower[tolower(gene_map_lower[[gene_symbol_col]]) == id_lower, ]

        if (nrow(match_row_lower) > 0 && !is.na(gene_map[[gene_symbol_col]][match_row_lower[1, ]]) && gene_map[[gene_symbol_col]][match_row_lower[1, ]] != "") {
          return(as.character(gene_map[[gene_symbol_col]][which(tolower(gene_map[[gene_symbol_col]]) == id_lower)[1]]))
        }

        # If all else fails, return the original ID (better than NA)
        return(original_id)
      })
      paste(converted, collapse = "/")
    })

    # Remove ID column as it's redundant with rownames
    if ("ID" %in% colnames(result@result)) {
      result@result <- result@result[, !colnames(result@result) %in% "ID"]
    }
  }

  return(result)
}

