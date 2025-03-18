server <- function(input, output, session) {
  
  source("functions.R")
  disable("use_custom_diffexp")
  
  # Reactive values for loaded data and comparisons
  data_loaded <- reactiveValues(seuratObj = NULL, clusterMat = NULL)
  comparisons <- reactiveVal(list())
  diffexp_all <- reactiveVal(NULL)
  diffexp_results <- reactiveVal(NULL)
  diffexp_status <- reactiveVal(NULL)
  diffexp_message <- reactiveVal("Waiting for input...")
  pathway_result <- reactiveVal(NULL)
  dotplot_pathway <- reactiveVal(NULL)
  timer <- reactiveTimer(1000)

  # Event to load data when the user clicks the load button
  observeEvent(input$load_data, {
    req(input$cellranger_out, input$cluster_csv, input$gene_expression_cutoff,
        input$spot_gene_cutoff, input$species, input$normalisation_method)
    
    folderCellRangerOut <- input$cellranger_out
    filenameCluster <- input$cluster_csv$datapath
    
    filter_results <- loadAndPreprocess(folderCellRangerOut, input$gene_expression_cutoff,
                                        input$spot_gene_cutoff, input$species,
                                        input$normalisation_method)
    data_loaded$seuratObj <- filter_results$seuratObj
    data_loaded$mart <- filter_results$mart
    data_loaded$clusterMat <- loadClusterMat(filenameCluster, data_loaded$seuratObj)
    data_loaded$seuratObj[[]]["clusterMat"] <- data_loaded$clusterMat
    
    Idents(data_loaded$seuratObj) <- data_loaded$seuratObj[[]]["clusterMat"][[1]]
    
    gene_expression_sums <- Matrix::rowSums(data_loaded$seuratObj[[Assays(data_loaded$seuratObj)]]$counts)
    ordered_genes <- names(sort(gene_expression_sums, decreasing = TRUE))
    
    sorted_clusters <- sort(unique(data_loaded$clusterMat[,1]))

    updateSelectizeInput(session, "gene_select", choices = ordered_genes, server = TRUE)
    updateSelectizeInput(session, "gene_select_dotplot", choices = ordered_genes, server = TRUE)
    updateSelectizeInput(session, "comparison_select", choices = sorted_clusters, server = TRUE)
    updateSelectizeInput(session, "selected_cluster", choices = sorted_clusters, selected = sorted_clusters[1])
    
    observe({
      session$sendCustomMessage("enhanceSelectize", "gene_select_dotplot")
    })
    
    # Update the filtered out information
    output$data_info <- renderPrint({
      cat("Seurat Object Dimensions:", dim(data_loaded$seuratObj), "\n")
      cat("Cluster Matrix Dimensions:", dim(data_loaded$clusterMat), "\n")
      cat("Filtered out genes:", filter_results$filtered_genes, "\n")
      cat("Filtered out spots:", filter_results$filtered_spots, "\n")
    })
    
    # Pre-calculate diffexp results while user is busy looking at something else
    process <- callr::r_bg(FindAllMarkers, 
      args = list(data_loaded$seuratObj),
      package = "Seurat")
    diffexp_status(process)
  })
  
  # Timer to observe when the diffexp is done
  observe({
    timer()
    process <- diffexp_status()
    if (!is.null(process) && process$is_alive()) {
      diffexp_message("Calculation of differential expression is running, this may take some time depending on hardware...")
    } else if (!is.null(process) && process$is_alive() == FALSE) {
      diffexp_all(process$get_result())
      diffexp_status(NULL)
      shinyjs::toggle("diffexp_message_box", anim = TRUE)
    }
    output$diffexp_message <- renderText({
      diffexp_message()
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
  
  #Update the saved Gene Table
  observeEvent(input$update_gene_table, {
    req(input$species)
    
    withProgress(message = "Updating gene table...", value = 0, {
      incProgress(0.4, detail = "Querying gene table...")
      mart <- checkMart(input$species)
      
      if (!is.null(mart)) {
        requestGeneTable(mart, input$species)
      } else {
        showNotification("Failed to load the mart", type = "error")
      }
      
      incProgress(0.4, detail = "Saving gene table...")
      Sys.sleep(0.2)
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
      paste0(Sys.Date(), "_violin_plot_", input$gene_select, ".pdf")
    },
    content = function(file) {
      gene_data <- prepare_gene_data(input$gene_select, data_loaded)
      
      pdf(file, width = 8, height = 8)
      par(mfrow = c(2, 1))
      plot(create_plot_with_stats(create_violin_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      plot(create_plot_with_stats(create_beeswarm_plot, gene_data, input$gene_select, comparisons(), input$display_pval))
      dev.off()
    }
  )
  
  # Heatmap & Dotplot
  heatmap_data <- eventReactive(input$run_heatmap, {
    req(input$gene_select_dotplot, data_loaded$seuratObj)
    
    withProgress(message = "Generating Heatmap...", value = 0, {
      
      incProgress(0.2, detail = "Extracting data")
      
      genes <- input$gene_select_dotplot
      gexp <- GetAssayData(data_loaded$seuratObj, slot = "scale.data")[genes, , drop = FALSE]
      
      sample_annot <- data_loaded$seuratObj[[]] %>% rownames_to_column("sample")
      
      incProgress(0.3, detail = "Processing data")
      
      dataHm <- as.data.frame(t(gexp)) %>%
        rownames_to_column("sample") %>%
        left_join(sample_annot, by = "sample") %>%
        tibble() %>%
        group_by(clusterMat)
      
      color_limit <- max(quantile(dataHm[, genes, drop = FALSE], 0.99, na.rm = TRUE),
                         -quantile(dataHm[, genes, drop = FALSE], 0.01, na.rm = TRUE))
      
      incProgress(0.2, detail = "Generating heatmap")
      
      hmplot <- ggheatmap(dataHm,
                          colv = "sample",
                          rowv = genes,
                          hm_colors = "RdBu",
                          scale = TRUE,
                          center = TRUE, 
                          hm_color_limits = c(-color_limit, color_limit),
                          show_dend_col = FALSE,
                          show_dend_row = FALSE,
                          show_colnames = FALSE,
                          show_rownames = TRUE,
                          colors_title = "Scaled expression (log2 UQ)") +
        plot_layout(guides = 'collect')
      
      incProgress(0.3, detail = "Done")
      hmplot
    })
  })
  
  
  dotplot_data <- eventReactive(input$run_heatmap, {
    req(input$gene_select_dotplot, data_loaded$seuratObj)
    
    withProgress(message = "Generating DotPlot...", value = 0, {
      incProgress(0.5, detail = "Processing data")
      
      plot <- DotPlot(
        object = data_loaded$seuratObj,
        features = input$gene_select_dotplot,
        group.by = "clusterMat"
      ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_viridis_c(option = "plasma")
      
      incProgress(0.5, detail = "Rendering plot")
      plot
    })
  })
  
  
  output$heatmapPlot <- renderPlot({
    heatmap_data()
  })
  
  output$dotPlot <- renderPlot({
    dotplot_data()
  })
  
  output$download_combined_pdf <- downloadHandler(
    filename = function() {
      first_genes <- input$gene_select_dotplot[!is.na(input$gene_select_dotplot[1:3])][1:3]
      paste0(Sys.Date(), "_heatmap_dotplot_", paste(first_genes, collapse = "_"), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 8, height = 12)
      par(mfrow = c(2, 1))
      heatmap <- DoHeatmap(
        object = data_loaded$seuratObj,
        features = input$gene_select_dotplot,
        group.by = "clusterMat"
      ) + scale_fill_viridis(option = "plasma")
      
      dotplot <- DotPlot(
        object = data_loaded$seuratObj,
        features = input$gene_select_dotplot,
        group.by = "clusterMat"
      ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_viridis_c(option = "plasma")
      
      print(heatmap)
      print(dotplot)
      
      dev.off()
    }
  )
  
  observe_trigger <- reactive({
    req(data_loaded$seuratObj)
    list(
      selected_cluster = input$selected_cluster,
      diffexp = diffexp_all()
    )
  })
  
  # diffexp
  observe({
    trigger_values <- observe_trigger()
    selected_cluster <- trigger_values$selected_cluster
    diffexp <- trigger_values$diffexp
    
    selected_cluster <- input$selected_cluster

    if (!is.null(diffexp_all())) {
      diffexp <- diffexp_all() |>
        filter(str_detect(cluster, selected_cluster)) |>
        dplyr::select(-contains("cluster"))
      
      if ("gene" %in% colnames(diffexp)) {
        rownames(diffexp) <- diffexp$gene
      } else {
        warning("The 'gene' column is missing in the differential expression results.")
      }
      
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
    }
  })
  
  output$download_diffexp <- downloadHandler(
    filename = function() {
      paste0("DiffExp_Cluster", input$selected_cluster, "_vs_all_others.csv")
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
  
  observeEvent(input$run_ct_enrichment, {
    withProgress(message = "Running Cell Type Enrichment Analysis...", value = 0, {
      req(diffexp_all(), data_loaded$mart, input$species, input$celltype_method, input$celltype_db)
      
      # Load and prepare data
      genes <- diffexp_all()
      gene_map <- convertGeneMap(genes, data_loaded$mart, input$species)
      Cell_marker <- read_xlsx("./data/Cell_marker_All.xlsx")
      
      pathways  <- Cell_marker |>
        filter(str_detect(species, input$species)) |>
        filter(!is.na(Symbol)) |>
        distinct() |>
        pull(Symbol, cell_name)
      pathways <- split(pathways, names(pathways))
      
      genes_sorted <- genes |>
        dplyr::mutate(cluster = factor(cluster, levels = paste("Cluster", sort(as.numeric(gsub("Cluster ", "", unique(cluster))))))) |>
        dplyr::arrange(cluster, desc(avg_log2FC))
      
      enrichment_results <- list()
      barplots_celltype <- list()
      
      # Number of clusters for progress calculation
      total_clusters <- length(unique(genes_sorted$cluster))
      cluster_progress_step <- 1 / total_clusters  # Calculate progress step per cluster
      
      for (clust in sort(unique(genes_sorted$cluster))) {
        # Update progress for each cluster iteration
        incProgress(cluster_progress_step, message = paste("Processing Cluster:", clust))
        
        cluster_genes <- genes_sorted[genes_sorted$cluster == clust, ]
        clust_ranks <- setNames(as.numeric(cluster_genes$avg_log2FC), rownames(cluster_genes))
        clust_ranks <- sort(clust_ranks[!duplicated(names(clust_ranks))], decreasing = TRUE)
        
        if (input$celltype_method == "FGSEA") {
          enrichment_output <- fgsea::fgsea(pathways = pathways, stats = clust_ranks)
          enrichment_results[[as.character(clust)]] <- enrichment_output |>
            as.data.frame() |>
            dplyr::filter(NES >= 0) |>
            dplyr::arrange(padj)
          barplots_celltype[[as.character(clust)]] <- ggplot(head(enrichment_results[[as.character(clust)]], 10), 
                                                             aes(x = reorder(pathway, -padj), y = NES, fill = padj))
        } else if (input$celltype_method == "Enrichr Web Query") {
          enrichment_output <- enrichR::enrichr(names(clust_ranks), databases = input$celltype_db)
          enrichment_results[[as.character(clust)]] <- enrichment_output[[1]] |>
            as.data.frame() |>
            dplyr::arrange(Adjusted.P.value)
            
          barplots_celltype[[as.character(clust)]] <- ggplot(head(enrichment_results[[as.character(clust)]], 10), 
                                                             aes(x = reorder(Term, -Adjusted.P.value), y = Combined.Score, fill = Adjusted.P.value))
        }
        
        barplots_celltype[[as.character(clust)]] <- barplots_celltype[[as.character(clust)]] +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_gradient(low = "midnightblue", high = "tan1", name = "Adjusted p-value", limits = c(0, 1)) +
          labs(x = "Cell Type", y = "Enrichment Score",
               title = paste0(input$celltype_method, " cell type enrichment on ", input$celltype_db, " database in ", input$species, " : ", clust)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, face = "bold"),
                axis.text.y = element_text(size = 10),
                axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12),
                panel.background = element_blank(),
                panel.grid.major = element_line(colour = "gray90")) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))
      }
      
      output$fgsea_plots <- renderUI({
        plotTabs <- lapply(names(barplots_celltype), function(clust) {
          tabPanel(
            title = clust,
            plotOutput(outputId = paste0(input$celltype_method, "_plot_", clust), height = "400px")
          )
        })
        do.call(tabsetPanel, c(id = "plot_tabs", plotTabs))
      })
      
      output$fgsea_tables <- renderUI({
        tableTabs <- lapply(names(enrichment_results), function(clust) {
          tabPanel(
            title = clust,
            DT::dataTableOutput(outputId = paste0(input$celltype_method, "_table_", clust))
          )
        })
        do.call(tabsetPanel, c(id = "table_tabs", tableTabs))
      })
      
      # Render each plot and data table for each cluster
      observe({
        for (clust in names(barplots_celltype)) {
          local({
            cluster_name <- clust
            
            # Render the plot for the current cluster
            output[[paste0(input$celltype_method, "_plot_", cluster_name)]] <- renderPlot({
              barplots_celltype[[cluster_name]]
            })
            
            # Render the data table for the current cluster
            output[[paste0(input$celltype_method, "_table_", cluster_name)]] <- DT::renderDataTable({
              enrichment_results[[cluster_name]]
            })
          })
        }
      })
      output$download_ct_plot_pdf <- downloadHandler(
        filename = function() {
          paste0(Sys.Date(), "_", input$plot_tabs, "_Cell_Type_Enrichment_Plot.pdf")
        },
        content = function(file) {
          selected_cluster <- input$plot_tabs
          selected_plot <- barplots_celltype[[selected_cluster]]
          pdf(file, width = 12, height = 6)
          print(selected_plot)
          dev.off()
        }
      )
      output$download_ct_data_csv <- downloadHandler(
        filename = function() {
          paste0(Sys.Date(), "_", input$table_tabs, "_Cell_Type_Enrichment_table.csv")
        },
        content = function(file) {
          selected_cluster <- input$table_tabs
          selected_table <- enrichment_results[[selected_cluster]]
          selected_table <- dplyr::mutate(selected_table, across(where(is.list), ~sapply(.x, paste, collapse = ", ")))
          write.csv(selected_table, file, row.names = FALSE)
        }
      )
    })
  })
  
  # Run pathway analysis
  observeEvent(input$run_pathway, {
    req(diffexp_data(), input$pathway_method, input$species, data_loaded$mart,
        input$selected_cluster, input$pathway_database)
    
    genes <- diffexp_data()
    method <- input$pathway_method
    species <- input$species
    mart <- data_loaded$mart
    chosenCluster <- input$selected_cluster
    database <- input$pathway_database
    
    withProgress(message = 'Running Pathway Analysis', value = 0, {
      result <- runPathwayAnalysis(genes, method, database, species, mart)
      pathway_result(result) 
      
      incProgress(0.1, detail = "Rendering Data table")
      output$pathway_results <- DT::renderDataTable({
        resultDt <- as.data.frame(result)
        DT::datatable(resultDt, options = list(pageLength = 20))
      })
      
      incProgress(0.1, detail = "Rendering plots")
      dotplot_pathway(dotplot(result, showCategory = 30) +
                        ggtitle(paste0(method, " method on ", database, " database in ", species, " : ", chosenCluster, " vs all")) +
                        theme(
                          plot.title = element_text(size = 12, face = "bold"),
                          axis.text.y = element_text(size = 10),
                          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                          axis.title.x = element_text(size = 12),
                          axis.title.y = element_text(size = 12),
                          panel.background = element_blank(),
                          panel.grid.major = element_line(colour = "gray90")
                        ) +
                        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)))
      output$pathway_plot <- renderPlot({
        dotplot_pathway()
      })
    })
  })
  
  # Download Pathway Plot as a PDF
  output$download_pathway_plot_pdf <- downloadHandler(
    filename = function() {
      paste("pathway_analysis_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 6, height = 9)
      
      print(dotplot_pathway())
      
      dev.off()
    }
  )
  
  
  # Download Pathway Data Table as CSV
  output$download_pathway_data_csv <- downloadHandler(
    filename = function() {
      paste("pathway_analysis_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      result <- pathway_result()
      if (!is.null(result)) {
        write.csv(as.data.frame(result), file, row.names = FALSE)
      }
    }
  )
}