server <- function(input, output, session) {
  source("functions.R")
  
  data_loaded <- reactiveValues(seuratObj = NULL, clusterMat = NULL)
  
  observeEvent(input$load_data, {
  req(input$cellranger_out, input$cluster_csv)
  
  folderCellRangerOut <- input$cellranger_out
  filenameCluster <- input$cluster_csv$datapath
  
  # Load and preprocess CellRanger output
  data_loaded$seuratObj <- loadAndPreprocess(folderCellRangerOut)
  
  # Print dimensions of the count matrix to ensure it's loaded correctly
  countMatrix <- data_loaded$seuratObj[["Spatial"]]$data
  print(paste("Dimensions of countMatrix:", dim(countMatrix)))
  
  # Load cluster matrix from CSV
  data_loaded$clusterMat <- loadClusterMat(filenameCluster, data_loaded$seuratObj)
  
  # Print dimensions of the cluster matrix to ensure it's loaded correctly
  print(paste("Dimensions of clusterMat:", dim(data_loaded$clusterMat)))
  
  # Check for common barcodes between countMatrix and clusterMat
  print(paste("Number of common barcodes:", length(intersect(colnames(countMatrix), data_loaded$clusterMat$Barcode))))
  
  # Update gene select input with loaded genes
  updateSelectizeInput(session, "gene_select", choices = rownames(data_loaded$seuratObj), server = TRUE)
})

  
  
  # Display data information in Data Loading tab
  output$data_info <- renderPrint({
    if (is.null(data_loaded$seuratObj) || is.null(data_loaded$clusterMat)) {
      cat("No data loaded yet.")
    } else {
      cat("Seurat Object Dimensions:", dim(data_loaded$seuratObj), "\n")
      cat("Cluster Matrix Dimensions:", dim(data_loaded$clusterMat), "\n")
    }
  })
  
  # Generate Violin plot
  output$violinPlot <- renderPlot({
    req(input$gene_select, data_loaded$seuratObj, data_loaded$clusterMat)  
    
    gene_data <- prepare_gene_data(input$gene_select, data_loaded)
    create_violin_plot(gene_data, input$gene_select)
  })
  
  # Generate Beeswarm plot
  output$beeswarmPlot <- renderPlot({
    req(input$gene_select, data_loaded$seuratObj, data_loaded$clusterMat)  # Ensure inputs are available
    
    gene_data <- prepare_gene_data(input$gene_select, data_loaded)
    create_beeswarm_plot(gene_data, input$gene_select)
  })
  
  # Download as PNG
  output$download_png <- downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      gene_data <- prepare_gene_data(input$gene_select, data_loaded)
      
      png(file, width = 800, height = 600)
      plot(create_violin_plot(gene_data, input$gene_select))
      plot(create_beeswarm_plot(gene_data, input$gene_select))
      dev.off()
    }
  )
  
  # Download as PDF
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      gene_data <- prepare_gene_data(input$gene_select, data_loaded)
      
      pdf(file, width = 8, height = 6)
      plot(create_violin_plot(gene_data, input$gene_select))
      plot(create_beeswarm_plot(gene_data, input$gene_select))
      dev.off()
    }
  )
  
}