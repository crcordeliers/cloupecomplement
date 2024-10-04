# Load UI and Server from external files
source("server.R")
source("ui.R")

# Run the app
shinyApp(ui = ui, server = server)
