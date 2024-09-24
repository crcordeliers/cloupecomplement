# Load UI and Server from external files
source("ui.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)
