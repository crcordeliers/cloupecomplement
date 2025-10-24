if (!require("BiocManager")) install.packages("BiocManager", quiet = TRUE)
if (!require("devtools")) install.packages("devtools", quiet = TRUE)
if (!require("pacman")) install.packages("pacman", quiet = TRUE)
if (!require("enrichR")) remotes::install_github("wjawaid/enrichR")
if (!require("presto")) remotes::install_github("immunogenomics/presto")
if (!require("ggheatmapper")) remotes::install_github("csgroen/ggheatmapper")

options(shiny.maxRequestSize=100*1024^2)

pacman::p_load(shiny, shinydashboard, ggplot2, shinyWidgets, dplyr, ggbeeswarm,
               Seurat, reshape2, ggpubr, ggheatmapper, viridis, clusterProfiler,
               org.Hs.eg.db, org.Mm.eg.db, biomaRt, fgsea, msigdbr, tidyverse, readxl, devtools,
               enrichR, callr, shinyjs, gtools, WriteXLS, harmony)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ui <- dashboardPage(
  dashboardHeader(title = "cLoupeComplement"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Loading", tabName = "data_loading", icon = icon("upload")),
      menuItem("Violin & Beeswarm Plots", tabName = "VB_plots", icon = icon("chart-simple")),
      menuItem("Heatmap & Dotplots", tabName = "hmap_dotplot", icon = icon("chart-bar")),
      menuItem("Diffexp", tabName = "diffexp", icon = icon("table")),
      menuItem("Celltype Enrichment", tabName = "celltype_enrichment", icon = icon("magnifying-glass")),
      menuItem("Pathway Analysis", tabName = "pathway_analysis", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(
      # Data loading tab
      tabItem(
        tabName = "data_loading",
        h2("Data Loading"),
        
        # data loading inputs
        fluidRow(
          box(
            width = 12,
            title = "Data Input",
            solidHeader = TRUE,
            status = "primary",

            fluidRow(
              column(6,
                     fileInput("h5_file", "Upload filtered_feature_bc_matrix.h5 file (max 100MB)",
                               accept = c(".h5", ".hdf5"), placeholder = "No file selected"),
                     helpText("For Visium HD: use binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
              ),
              column(6,
                     fileInput("cluster_csv", "Choose Cluster CSV File (max 100MB)",
                               accept = ".csv", placeholder = "No file selected")
              )
            )
          )
        ),
        
        # preprocessing options
        fluidRow(
          box(
            width = 12,
            title = "Preprocessing",
            solidHeader = TRUE,
            status = "warning",
            
            selectInput("species", "Select Species",
                        choices = c("Human", "Mouse"),
                        selected = "Human"),
            selectInput("normalisation_method", "Select normalisation method",
                        choices = c("LogNormalize", "SCTransform"),
                        selected = "LogNormalize"),
            numericInput("gene_expression_cutoff",
                         "Minimum % of cells expressing gene:",
                         value = 1, min = 0, max = 100, step = 1),
            numericInput("spot_gene_cutoff",
                         "Minimum number of genes expressed per spot:",
                         value = 100, min = 0, step = 10),

            # Visium HD options
            checkboxInput("is_visium_hd", "Visium HD Dataset (Use Sketching)", value = FALSE),
            conditionalPanel(
              condition = "input.is_visium_hd == true",
              numericInput("sketch_size",
                           "Number of spots to sketch (for memory optimization):",
                           value = 5000, min = 1000, max = 50000, step = 1000),
              helpText("Sketching randomly samples spots to reduce memory usage for large Visium HD datasets. 5,000-10,000 spots typically provides good balance between performance and quality.")
            )
          )
        ),
        actionButton("load_data", "Load Data"),
        br(),
        verbatimTextOutput("data_info"),

        # Info box for Visium HD sketching
        conditionalPanel(
          condition = "input.is_visium_hd == true",
          fluidRow(
            box(
              width = 12,
              title = "Visium HD Sketching Information",
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,

              HTML("<p><strong>Operations using SKETCHED data (faster, memory efficient):</strong></p>
                   <ul>
                     <li>Violin & Beeswarm plots</li>
                     <li>Heatmaps</li>
                     <li>DotPlots</li>
                     <li>Data visualization and exploration</li>
                   </ul>
                   <p><strong>Operations using FULL data (accurate, slower):</strong></p>
                   <ul>
                     <li>Differential expression analysis (FindAllMarkers)</li>
                     <li>Statistical testing</li>
                   </ul>
                   <p><em>This approach ensures visualizations are responsive while maintaining statistical accuracy for differential expression.</em></p>")
            )
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Advanced",
            solidHeader = TRUE,
            status = "danger",
            collapsible = TRUE,
            collapsed = TRUE,
            
            fluidRow(
              column(6,
                     actionButton("update_mart", "Update Mart")
                     ),
              column(6,
                     actionButton("update_gene_table", "Update Gene Table")
              )
            )
          )
        )
      ),
      
      # Violin & Beeswarm Plots tab
      tabItem(tabName = "VB_plots",
              h2("Violin & Beeswarm Plots"),
              
              # download buttons
              fluidRow(
                column(6, 
                       downloadButton("download_pdf", "Download as PDF"))
              ),
              br(),
              # gene selection and options
              fluidRow(
                box(
                  width = 12,
                  title = "Plot Options",
                  solidHeader = TRUE,
                  status = "info",
                  
                  selectizeInput("gene_select", "Select Gene of Interest:",
                                 choices = NULL, multiple = FALSE),
                  checkboxInput("show_comparisons", "Show pairwise comparisons", FALSE)
                )
              ),
              
              # Conditional panel for comparisons
              conditionalPanel(
                condition = "input.show_comparisons == true",
                fluidRow(
                  box(
                    width = 12,
                    title = "Comparisons",
                    solidHeader = TRUE,
                    status = "warning",
                    
                    fluidRow(
                      column(4, 
                             selectizeInput("comparison_select", "Select clusters to compare:",
                                            choices = NULL, multiple = TRUE, 
                                            options = list(maxItems = 2))),
                      column(4, 
                             checkboxInput("display_pval", "Show comparisons as pval", FALSE)),
                      column(4, 
                             actionButton("add_comparison", "Add Comparison", class = "btn-success"),
                             actionButton("remove_comparison", "Remove Last Comparison", class = "btn-danger"))
                    ),
                    verbatimTextOutput("current_comparisons")
                  )
                )
              ),
              
              # display plots
              fluidRow(
                box(
                  width = 6, 
                  title = "Violin Plot", 
                  solidHeader = TRUE, 
                  status = "primary",
                  plotOutput("violinPlot")
                ),
                box(
                  width = 6, 
                  title = "Beeswarm Plot", 
                  solidHeader = TRUE, 
                  status = "primary",
                  plotOutput("beeswarmPlot")
                )
              )
      ),
      
      # Heatmap & Dotplot tab
      tabItem(tabName = "hmap_dotplot",
              h2("Heatmap & Dotplots"),
              
              downloadButton("download_combined_pdf", "Download as PDF"),
              br(), br(),
              
              # Gene selection
              fluidRow(
                box(
                  width = 12,
                  title = "Plot Options",
                  solidHeader = TRUE,
                  status = "info",
                  
                  selectizeInput("gene_select_dotplot", "Select Genes of Interest:",
                                 choices = NULL, multiple = TRUE),
                  actionButton("run_heatmap", "Run Heatmap", class = "btn-primary")
                )
              ),
              
              tags$script(HTML("
      Shiny.addCustomMessageHandler('enhanceSelectize', function(inputId) {
        var $select = $('#' + inputId);
        var selectize = $select[0].selectize;

        if (selectize) {
          // Remove any existing paste handler
          selectize.$control_input.off('paste');

          // Add paste event handler
          selectize.$control_input.on('paste', function(e) {
            e.preventDefault();
            var pastedData = (e.originalEvent || e).clipboardData.getData('text/plain');

            // Split by newlines, commas, semicolons, or spaces
            var genes = pastedData.split(/[\\n\\r,;\\s]+/).filter(function(item) {
              return item.trim().length > 0;
            });

            // Add each gene as an item
            genes.forEach(function(gene) {
              var trimmedGene = gene.trim();
              // Check if the gene exists in available options
              if (selectize.options[trimmedGene]) {
                selectize.addItem(trimmedGene, true);
              } else {
                // Try to add it anyway (will work if server-side selectize allows it)
                selectize.addOption({value: trimmedGene, text: trimmedGene});
                selectize.addItem(trimmedGene, true);
              }
            });

            selectize.refreshOptions(false);
          });
        }
      });
    ")),
              
              # Results
              fluidRow(
                box(
                  title = "Heatmap & Dotplot Results", 
                  width = 12, 
                  solidHeader = TRUE, 
                  status = "primary",
                  
                  tabsetPanel(
                    id = "hmap_dotplot_tabs",
                    
                    tabPanel(
                      title = "Heatmap", 
                      value = "heatmap",
                      br(), br(),
                      plotOutput("heatmapPlot", height = "700px")
                    ),
                    
                    tabPanel(
                      title = "Dotplot", 
                      value = "dotplot",
                      br(), br(),
                      plotOutput("dotPlot", height = "700px")
                    )
                  )
                )
              )
      ),
      
      
      # Diffexp tab
      tabItem(
        useShinyjs(),
        tabName = "diffexp",
        h2("Differential Expression Analysis"),
        
        downloadButton("download_all_diffexp", "Download all Results"),
        br(), br(),
        
        # cluster selection
        fluidRow(
          box(
            width = 12,
            title = "Analysis Options",
            solidHeader = TRUE,
            status = "info",
            
            selectInput("selected_cluster", "Select Cluster:", choices = NULL, selected = NULL)
          )
        ),
        
        # Status message
        fluidRow(
          
          box(
            id = "diffexp_message_box",
            width = 12,
            status = "warning",
            solidHeader = TRUE,
            textOutput("diffexp_message")
          )
        ),
        
        # diffexp table
        fluidRow(
          box(
            width = 12,
            title = "Results",
            solidHeader = TRUE,
            status = "primary",
            
            downloadButton("download_diffexp", "Download this cluster's results"),
            br(), br(),
            
            DT::dataTableOutput("diffexp_table")
          )
        )
      ),
      
      tabItem(
        tabName = "celltype_enrichment",
        
        fluidRow(
          box(
            title = "Cell Type Enrichment Settings", 
            width = 12, 
            solidHeader = TRUE, 
            status = "info",
            
            selectInput("celltype_method", "Select Method:", choices = c("FGSEA", "Enrichr Web Query"), selected = "FGSEA"),
            selectInput("celltype_db", "Select Database:", choices = c("CellMarker_2024"), selected = "CellMarker_2024"),
            actionButton("run_ct_enrichment", "Run Cell Type Enrichment", class = "btn-primary")
          )
        ),
        
        fluidRow(
          box(
            title = "Cell Type Enrichment Results", 
            width = 12, 
            solidHeader = TRUE, 
            status = "primary",
            
            tabsetPanel(
              id = "ct_results_tabs",
              
              tabPanel(
                title = "Plot", 
                value = "plot",
                br(),
                downloadButton("download_ct_plot_pdf", "Download Cell Type Plot as PDF"),
                br(), br(),
                uiOutput("fgsea_plots")
              ),
              
              tabPanel(
                title = "Data Table", 
                value = "table",
                br(),
                downloadButton("download_ct_data_csv", "Download Cell Type Data Table as CSV"),
                br(), br(),
                uiOutput("fgsea_tables")
              )
            )
          )
        )
      ),
      
      
      # Pathway Analysis tab
      tabItem(tabName = "pathway_analysis",
              fluidRow(
                box(
                  title = "Pathway Analysis Settings", 
                  width = 12, 
                  solidHeader = TRUE, 
                  status = "info",
                  
                  # custom diffexp checkbox & file input
                  checkboxInput("use_custom_diffexp", "Use Custom Differential Expression Results", value = FALSE),
                  conditionalPanel(
                    condition = "input.use_custom_diffexp == true",
                    fileInput("diffexp_file", "Upload Differential Expression CSV", accept = ".csv")
                  ),
                  
                  selectInput("pathway_method", "Select Method:", choices = c("ORA", "FGSEA"),
                              selected = "ORA"),
                  selectInput("pathway_database", "Select Database:", choices = c("GO", "KEGG", "HALLMARK"),
                              selected = "GO"),
                  actionButton("run_pathway", "Run Analysis", class = "btn-primary")
                ),
                
                # Results
                box(
                  title = "Pathway Analysis Results", 
                  width = 12, 
                  solidHeader = TRUE, 
                  status = "primary",
                  
                  tabsetPanel(
                    id = "results_tabs",
                    
                    tabPanel(
                      title = "Plot", 
                      value = "plot",
                      br(),
                      downloadButton("download_pathway_plot_pdf", "Download Pathway Plot as PDF"),
                      br(),br(),
                      plotOutput("pathway_plot", height = "700px")
                    ),
                    
                    tabPanel(
                      title = "Data Table", 
                      value = "table",
                      br(),
                      downloadButton("download_pathway_data_csv", "Download Pathway Data Table as CSV"),
                      br(),br(),
                      DT::dataTableOutput("pathway_results")
                    )
                  )
                )
              )
      )
    )
  )
)
