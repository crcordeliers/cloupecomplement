if (!require("pacman")) install.packages("pacman", quiet = TRUE)
pacman::p_load(shiny, shinydashboard, ggplot2, shinyWidgets, dplyr, ggbeeswarm,
               Seurat, reshape2, ggpubr, pheatmap, viridis)

ui <- dashboardPage(
  dashboardHeader(title = "cLoupeComplement"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Loading", tabName = "data_loading", icon = icon("upload")),
      menuItem("Violin & Beeswarm Plots", tabName = "plots", icon = icon("chart-simple")),
      menuItem("Heatmap & Dotplots", tabName = "hmap_dotplot", icon = icon("chart-bar")),
      menuItem("Diffexp", tabName = "diffexp", icon = icon("table")),
      menuItem("Pathway Analysis", tabName = "pathway", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(
      # Data loading tab
      tabItem(
        tabName = "data_loading",
        h2("Data Loading"),
        textInput("cellranger_out", "Enter CellRanger Output Folder Path", placeholder = "/path/to/cellranger/outs"),
        fileInput("cluster_csv", "Choose Cluster CSV File", accept = ".csv", placeholder = "No file selected"),
        numericInput("gene_expression_cutoff", "Minimum % of cells expressing gene:", value = 1, min = 0, max = 100, step = 1),
        numericInput("spot_gene_cutoff", "Minimum number of genes expressed per spot:", value = 100, min = 0, step = 10),
        actionButton("load_data", "Load Data"),
        verbatimTextOutput("data_info")
      ),
      
      # Violin & Beeswarm Plots tab
      tabItem(tabName = "plots",
              h2("Violin & Beeswarm Plots"),
              fluidRow(
                column(6, downloadButton("download_png", "Download as PNG")),
                column(6, downloadButton("download_pdf", "Download as PDF"))
              ),
              selectizeInput("gene_select", "Select Gene of Interest:",
                             choices = NULL, multiple = FALSE),
              checkboxInput("show_comparisons", "Show pairwise comparisons", FALSE),
              conditionalPanel(
                condition = "input.show_comparisons == true",
                fluidRow(
                  column(4, selectizeInput("comparison_select", "Select clusters to compare:",
                                 choices = NULL, multiple = TRUE, options = list(maxItems = 2))),
                  column(4, checkboxInput("display_pval", "Show comparisons as pval", FALSE)),
                ),
                actionButton("add_comparison", "Add Comparison"),
                actionButton("remove_comparison", "Remove Last Comparison"),
                verbatimTextOutput("current_comparisons")
              ),
              plotOutput("violinPlot"),
              plotOutput("beeswarmPlot")
      ),
      
      # Heatmap & Dotplot tab
      tabItem(tabName = "hmap_dotplot",
              h2("Heatmap & Dotplots"),
              selectizeInput("gene_select_dotplot", "Select Genes of Interest:",
                             choices = NULL, multiple = TRUE),
              plotOutput("heatmapPlot"),
              plotOutput("dotPlot")
      ),
      
      # Diffexp tab
      tabItem(
        tabName = "diffexp",
        h2("Differential Expression Analysis"),
        selectInput("selected_cluster", "Select Cluster:", choices = NULL, selected = NULL),
        br(), br(),
        DT::dataTableOutput("diffexp_table"),
        downloadButton("download_diffexp", "Download Differential Expression Results")
      ),
      
      # Pathway Analysis tab (Work in progress)
      tabItem(tabName = "pathway",
              h2("Pathway Analysis - Work in Progress")
      )
    )
  )
)