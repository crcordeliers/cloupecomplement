server <- function(input, output, session) {
  
  source("functions.R")
  
  # Reactive values for loaded data and comparisons
  data_loaded <- reactiveValues(seuratObj = NULL, clusterMat = NULL)
  comparisons <- reactiveVal(list())
  diffexp_results <- reactiveVal(data.frame())
  
  # Event to load data when the user clicks the load button
  observeEvent(input$load_data, {
    req(input$cellranger_out, input$cluster_csv, input$gene_expression_cutoff, input$spot_gene_cutoff, input$species)
    
    folderCellRangerOut <- input$cellranger_out
    filenameCluster <- input$cluster_csv$datapath
    
    filter_results <- loadAndPreprocess(folderCellRangerOut, input$gene_expression_cutoff, input$spot_gene_cutoff, input$species)
    data_loaded$seuratObj <- filter_results$seuratObj
    data_loaded$mart <- filter_results$mart
    data_loaded$clusterMat <- loadClusterMat(filenameCluster, data_loaded$seuratObj)
    data_loaded$seuratObj[[]]["clusterMat"] <- data_loaded$clusterMat
    
    Idents(data_loaded$seuratObj) <- data_loaded$seuratObj[[]]["clusterMat"][[1]]
    
    gene_expression_sums <- Matrix::rowSums(data_loaded$seuratObj[["Spatial"]]$counts)
    ordered_genes <- names(sort(gene_expression_sums, decreasing = TRUE))
    
    sorted_clusters <- sort(unique(data_loaded$clusterMat[,1]))
    updateSelectizeInput(session, "gene_select", choices = ordered_genes, server = TRUE)
    updateSelectizeInput(session, "gene_select_dotplot", choices = ordered_genes, server = TRUE)
    updateSelectizeInput(session, "comparison_select", choices = sorted_clusters, server = TRUE)
    updateSelectizeInput(session, "selected_cluster", choices = sorted_clusters, selected = sorted_clusters[1])
    
    # Update the filtered out information
    output$data_info <- renderPrint({
      cat("Seurat Object Dimensions:", dim(data_loaded$seuratObj), "\n")
      cat("Cluster Matrix Dimensions:", dim(data_loaded$clusterMat), "\n")
      cat("Filtered out genes:", filter_results$filtered_genes, "\n")
      cat("Filtered out spots:", filter_results$filtered_spots, "\n")
    })
  })
  
  output$data_info <- renderPrint({
    if (is.null(data_loaded$seuratObj) || is.null(data_loaded$clusterMat)) {
      cat("No data loaded yet.")
    } else {
      cat("Seurat Object Dimensions:", dim(data_loaded$seuratObj), "\n")
      cat("Cluster Matrix Dimensions:", dim(data_loaded$clusterMat), "\n")
    }
  })
  
  # Update the saved Mart file
  observeEvent(input$update_mart, {
    req(input$species)
    withProgress(message = "Updating mart...", value = 0, {
      incProgress(0.4, detail = "Loading mart...")

      checkMart(input$species, updateMart = TRUE)
      incProgress(0.4, detail = "Saving mart...")
      Sys.sleep(1)
    })
  })
  
  # Add a comparison to the list
  observeEvent(input$add_comparison, {
    req(input$comparison_select)
    
    if (length(input$comparison_select) == 2) {
      current <- comparisons()
      current[[length(current) + 1]] <- as.character(input$comparison_select)
      comparisons(current)
      
      # Clear comparison select input after adding comparison
      updateSelectizeInput(session, "comparison_select", selected = character(0))
    }
  })
  
  # Remove the most recent comparison
  observeEvent(input$remove_comparison, {
    current <- comparisons()
    if (length(current) > 0) {
      comparisons(current[-length(current)])
    }
  })
  
  # Display current comparisons
  output$current_comparisons <- renderPrint({
    cat("Current comparisons:\n")
    for (i in seq_along(comparisons())) {
      cat(paste0(i, ": ", paste(comparisons()[[i]], collapse = " vs "), "\n"))
    }
  })
  
  # Render the violin plot with comparisons
  output$violinPlot <- renderPlot({
    req(input$gene_select, data_loaded$seuratObj, data_loaded$clusterMat)
    gene_data <- prepare_gene_data(input$gene_select, data_loaded)
    create_plot_with_stats(create_violin_plot, gene_data, input$gene_select, comparisons(), input$display_pval)
  })
  
  # Render the beeswarm plot with comparisons
  output$beeswarmPlot <- renderPlot({
    req(input$gene_select, data_loaded$seuratObj, data_loaded$clusterMat)
    gene_data <- prepare_gene_data(input$gene_select, data_loaded)
    create_plot_with_stats(create_beeswarm_plot, gene_data, input$gene_select, comparisons(), input$display_pval)
  })
  
  # Download PDF plot (violin + beeswarm)
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      gene_data <- prepare_gene_data(input$gene_select, data_loaded)
      
      pdf(file, width = 8, height = 12)
      par(mfrow = c(2, 1))
      plot(create_plot_with_stats(create_violin_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      plot(create_plot_with_stats(create_beeswarm_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      dev.off()
    }
  )
  
  # Heatmap & Dotplot
  output$heatmapPlot <- renderPlot({
    req(input$gene_select_dotplot, data_loaded$seuratObj)
    
    DoHeatmap(
      object = data_loaded$seuratObj,
      features = input$gene_select_dotplot,
      group.by = "clusterMat"
    ) + scale_fill_viridis(option = "plasma")
  })
  
  output$dotPlot <- renderPlot({
    req(input$gene_select_dotplot, data_loaded$seuratObj)
    
    DotPlot(
      object = data_loaded$seuratObj,
      features = input$gene_select_dotplot,
      group.by = "clusterMat"
    ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_viridis_c(option = "plasma")
  })
  
  output$download_combined_pdf <- downloadHandler(
    filename = function() {
      paste("heatmap_dotplot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 8, height = 12)
      
      DoHeatmap(
        object = data_loaded$seuratObj,
        features = input$gene_select_dotplot,
        group.by = "clusterMat"
      ) + scale_fill_viridis(option = "plasma")
      
      dev.flush()
      
      DotPlot(
        object = data_loaded$seuratObj,
        features = input$gene_select_dotplot,
        group.by = "clusterMat"
      ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_viridis_c(option = "plasma")
      
      dev.off()
    }
  )
  
  # diffexp
  observeEvent(input$selected_cluster, {
    req(data_loaded$seuratObj, input$selected_cluster)
    
    selected_cluster <- input$selected_cluster
    diffexp <- FindMarkers(data_loaded$seuratObj, ident.1 = selected_cluster, ident.2 = NULL)
    
    if ("p_val" %in% colnames(diffexp)) {
      diffexp$p_val <- format_pval(diffexp$p_val)
    }
    if ("p_val_adj" %in% colnames(diffexp)) {
      diffexp$p_val_adj <- format_pval(diffexp$p_val_adj)
    }
    if ("avg_log2FC" %in% colnames(diffexp)) {
      diffexp$avg_log2FC <- format_pval(diffexp$avg_log2FC)
    }
    
    diffexp_results(diffexp)
    
    output$diffexp_table <- DT::renderDataTable({
      DT::datatable(diffexp_results(), options = list(pageLength = 10, autoWidth = TRUE))
    })
  })
  
  output$download_diffexp <- downloadHandler(
    filename = function() {
      paste("DiffExp_Cluster", input$selected_cluster, "_vs_all_others.csv", sep = "")
    },
    content = function(file) {
      write.csv(diffexp_results(), file)
    }
  )
  
  diffexp_data <- reactive({
    if (input$use_custom_diffexp) {
      req(input$diffexp_file)
      df <- read.csv(input$diffexp_file$datapath, header = TRUE, sep = ",")
      return(df)
    } else {
      req(diffexp_results)
      return(diffexp_results())
    }
  })
  
  # Run pathway analysis
  observeEvent(input$run_pathway, {
    req(diffexp_data(), input$pathway_method, input$species, data_loaded$mart)
    
    genes <- rownames(diffexp_data())
    method <- input$pathway_method
    species <- input$species
    mart <- data_loaded$mart
    
    if (is.null(genes) || length(genes) == 0) {
      showNotification("No genes found for pathway analysis", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Running Pathway Analysis', value = 0, {
      result <- runPathwayAnalysis(genes, method, species, mart)
      
      # Check if result is NULL
      if (is.null(result)) {
        showNotification("Pathway analysis failed, please check your input or connection.", type = "error")
        return(NULL)
      }
      
      output$pathway_results <- DT::renderDataTable({
        resultDt <- as.data.frame(result)
        DT::datatable(resultDt, options = list(pageLength = 20))
      })
      
      output$pathway_plot <- renderPlot({
        if (method == "clusterProfiler") {
          if (inherits(result, "enrichResult")) {
            dotplot(result, showCategory = 20)
          } else {
            plot(1, 1, main = "Error: Result not compatible")
          }
        } else if (method == "fgsea") {
          if ("NES" %in% colnames(result)) {
            ggplot(result, aes(x = reorder(pathway, NES), y = NES)) +
              geom_bar(stat = "identity") +
              coord_flip()
          } else {
            plot(1, 1, main = "Error: Result not compatible with barplot")
          }
        }
      })
    })
  })
  
  # Download Pathway Plot as a PDF
  output$download_pathway_plot_pdf <- downloadHandler(
    filename = function() {
      paste("pathway_analysis_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 8, height = 6)
      
      if (input$pathway_method == "clusterProfiler") {
        result <- runPathwayAnalysis(rownames(diffexp_data()), input$pathway_method)
        if (inherits(result, "enrichResult")) {
          barplot(result, showCategory = 20)
        } else {
          plot(1, 1, main = "Error: Result not compatible with barplot")
        }
      } else if (input$pathway_method == "fgsea") {
        result <- runPathwayAnalysis(rownames(diffexp_data()), input$pathway_method)
        if ("NES" %in% colnames(result)) {
          ggplot(result, aes(x = reorder(pathway, NES), y = NES)) +
            geom_bar(stat = "identity") +
            coord_flip()
        } else {
          plot(1, 1, main = "Error: Result not compatible with barplot")
        }
      }
      
      dev.off()
    }
  )
  
  # Download Pathway Data Table as CSV
  output$download_pathway_data_csv <- downloadHandler(
    filename = function() {
      paste("pathway_analysis_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      result <- runPathwayAnalysis(rownames(diffexp_data()), input$pathway_method)
      write.csv(as.data.frame(result), file, row.names = FALSE)
    }
  )
}