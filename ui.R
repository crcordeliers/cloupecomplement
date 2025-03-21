if (!require("BiocManager")) install.packages("BiocManager", quiet = TRUE)
if (!require("devtools")) install.packages("devtools", quiet = TRUE)
if (!require("pacman")) install.packages("pacman", quiet = TRUE)
if (!require("enrichR")) install_github("wjawaid/enrichR")
if (!require("presto")) install_github("immunogenomics/presto")
if (!require("ggheatmapper")) install_github("csgroen/ggheatmapper")

pacman::p_load(shiny, shinydashboard, ggplot2, shinyWidgets, dplyr, ggbeeswarm,
               Seurat, reshape2, ggpubr, ggheatmapper, viridis, clusterProfiler,
               org.Hs.eg.db, biomaRt, fgsea, msigdbr, tidyverse, readxl, devtools,
               enrichR, callr, shinyjs, gtools, WriteXLS)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
                     textInput("cellranger_out", "Enter CellRanger Output Folder Path", 
                               placeholder = "/path/to/cellranger/outs")
              ),
              column(6, 
                     fileInput("cluster_csv", "Choose Cluster CSV File", 
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
                         value = 100, min = 0, step = 10)
          )
        ),
        actionButton("load_data", "Load Data"),
        br(),
        verbatimTextOutput("data_info"),
        
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
        var selectize = $('#' + inputId).selectize()[0].selectize;
        
        if (selectize) {
          selectize.onPaste = function(e) {
            var pastedData = (e.originalEvent || e).clipboardData.getData('text');
            var genes = pastedData.split(/[\\s,]+/).filter(Boolean);
            selectize.addItems(genes);
            e.preventDefault();
          };
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
