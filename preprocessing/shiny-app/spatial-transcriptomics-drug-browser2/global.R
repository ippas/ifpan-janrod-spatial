# read function
source("/home/rstudio/preprocessing/functions/functions-spatial-data.R")
source("/home/rstudio/preprocessing/functions/statistics-functions.R")
source("/home/rstudio/preprocessing/functions/visualization-functions.R")
source("/home/rstudio/preprocessing/functions/umi-per-spot.R")
source("/home/rstudio/preprocessing/functions/spatial_cluster_select.R")

require(shiny)
require(shinydashboard)


load("/home/rstudio/results/risperidone/risperidone-half.RData")
load("/home/rstudio/results/pz1190/pz1190-half.RData")
load("/home/rstudio/results/clozapine/clozapine-half.RData")

# Define the DPI
dpi <- 100
