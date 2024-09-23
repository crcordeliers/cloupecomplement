if (!require("pacman")) install.packages("pacman", quiet = TRUE)
pacman::p_load(shiny, shinydashboard, ggplot2, shinyWidgets, dplyr, ggbeeswarm,
               Seurat, reshape2)

# Load UI and Server from external files
source("ui.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)
