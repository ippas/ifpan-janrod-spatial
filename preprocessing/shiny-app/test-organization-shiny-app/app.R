# Load necessary packages
library(shiny)

# Source the separate scripts
source("global.R")
source("ui.R")
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
