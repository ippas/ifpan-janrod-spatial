# Load necessary packages
.libPaths("/usr/local/lib/R/site-library")

library(shiny)

# Source the separate scripts
source("global.R")
source("ui.R")
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
