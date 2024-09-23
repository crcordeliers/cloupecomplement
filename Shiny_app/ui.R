library(shiny)
library(shinydashboard)

# Define UI for the app
ui <- dashboardPage(
  dashboardHeader(title = "LoupeBrowserComplementaryAnalyses"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Loading", tabName = "data_loading", icon = icon("upload")),
      menuItem("Violin & Beeswarm Plots", tabName = "plots", icon = icon("chart-bar")),
      menuItem("Diffexp", tabName = "diffexp", icon = icon("table")),
      menuItem("Pathway Analysis", tabName = "pathway", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(
      # Data loading tab
      tabItem(tabName = "data_loading",
              h2("Data Loading"),
              textInput("cellranger_out", "Enter CellRanger Output Folder Path", 
                        placeholder = "/path/to/cellranger/outs"),
              fileInput("cluster_csv", "Choose Cluster CSV File", 
                        accept = ".csv", placeholder = "No file selected"),
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
              plotOutput("violinPlot"),
              plotOutput("beeswarmPlot")
      ),
      
      # Diffexp tab (Work in progress)
      tabItem(tabName = "diffexp",
              h2("Diffexp - Work in Progress")
      ),
      
      # Pathway Analysis tab (Work in progress)
      tabItem(tabName = "pathway",
              h2("Pathway Analysis - Work in Progress")
      )
    )
  )
)
