ui <- dashboardPage(
  dashboardHeader(title = "cLoupeComplement"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Loading", tabName = "data_loading", icon = icon("upload")),
      menuItem("Violin & Beeswarm Plots", tabName = "plots", icon = icon("chart-simple")),
      menuItem("Heatmap & Dotplots", tabName = "hmap_dotplot", icon = icon("chart-bar")),
      menuItem("Diffexp", tabName = "diffexp", icon = icon("table")),
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
            ),
          )
        ),
        
        # preprocessing options
        fluidRow(
          box(
            width = 12,
            title = "Preprocessing",
            solidHeader = TRUE,
            status = "warning",
            
            selectInput("species", "Select Species", choices = c("Human", "Mouse"), selected = "Human"),
            numericInput("gene_expression_cutoff", "Minimum % of cells expressing gene:", 
                         value = 1, min = 0, max = 100, step = 1),
            numericInput("spot_gene_cutoff", "Minimum number of genes expressed per spot:", 
                         value = 100, min = 0, step = 10)
          )
        ),
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
      
      # Pathway Analysis tab
      tabItem(tabName = "pathway_analysis",
              fluidRow(
                box(title = "Pathway Analysis Settings", width = 4,
                    checkboxInput("use_custom_diffexp", "Use Custom Differential Expression Results", value = FALSE),
                    conditionalPanel(
                      condition = "input.use_custom_diffexp == true",
                      fileInput("diffexp_file", "Upload Differential Expression CSV", accept = ".csv")
                    ),
                    selectInput("pathway_method", "Select Method:", choices = c("clusterProfiler", "fgsea")),
                    actionButton("run_pathway", "Run Analysis")
                ),
                box(title = "Pathway Analysis Results", width = 8,
                    tableOutput("pathway_results"),
                    plotOutput("pathway_plot")
                )
              )
      )
    )
  )
)