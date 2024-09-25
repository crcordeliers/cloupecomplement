server <- function(input, output, session) {
  # Load functions from functions.R
  source("functions.R")
  
  # Reactive values for loaded data and comparisons
  data_loaded <- reactiveValues(seuratObj = NULL, clusterMat = NULL)
  comparisons <- reactiveVal(list())
  
  # Event to load data when the user clicks the load button
  observeEvent(input$load_data, {
    req(input$cellranger_out, input$cluster_csv)
    
    folderCellRangerOut <- input$cellranger_out
    filenameCluster <- input$cluster_csv$datapath
    
    data_loaded$seuratObj <- loadAndPreprocess(folderCellRangerOut)
    data_loaded$clusterMat <- loadClusterMat(filenameCluster, data_loaded$seuratObj)
    
    sorted_clusters <- sort(unique(data_loaded$clusterMat[,1]))
    updateSelectizeInput(session, "gene_select", choices = rownames(data_loaded$seuratObj), server = TRUE)
    updateSelectizeInput(session, "gene_select_dotplot", choices = rownames(data_loaded$seuratObj), server = TRUE)
    updateSelectizeInput(session, "comparison_select", choices = sorted_clusters, server = TRUE)
  })
  
  output$data_info <- renderPrint({
    if (is.null(data_loaded$seuratObj) || is.null(data_loaded$clusterMat)) {
      cat("No data loaded yet.")
    } else {
      cat("Seurat Object Dimensions:", dim(data_loaded$seuratObj), "\n")
      cat("Cluster Matrix Dimensions:", dim(data_loaded$clusterMat), "\n")
    }
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
  
  # Download PNG plot (violin + beeswarm)
  output$download_png <- downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      gene_data <- prepare_gene_data(input$gene_select, data_loaded)
      
      png(file, width = 800, height = 1200)
      par(mfrow = c(2, 1))
      plot(create_plot_with_stats(create_violin_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      plot(create_plot_with_stats(create_beeswarm_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      dev.off()
    }
  )
  
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
}
