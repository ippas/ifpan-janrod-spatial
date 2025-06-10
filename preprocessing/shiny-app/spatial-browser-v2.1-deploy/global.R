# read function
source("/home/rstudio/preprocessing/functions/functions-spatial-data.R")
source("/home/rstudio/preprocessing/functions/statistics-functions.R")
source("/home/rstudio/preprocessing/functions/visualization-functions.R")
source("/home/rstudio/preprocessing/functions/umi-per-spot.R")
source("/home/rstudio/preprocessing/functions/spatial_cluster_select.R")

require(shiny)
require(shinydashboard)


load("/home/rstudio/results/risperidone/risperidone-half.RData")
message("Loaded: risperidone-half.RData")

load("/home/rstudio/results/pz1190/pz1190-half.RData")
message("Loaded: pz1190-half.RData")

load("/home/rstudio/results/clozapine/clozapine-half.RData")
message("Loaded: clozapine-half.RData")

# load("/home/rstudio/results/ldopa/ldopa.RData")
# message("Loaded: ldopa.RData")

# load("/home/rstudio/results/risperidone-3q29/risWtSalWt.RData")
# message("Loaded: risWtSalWt.RData")

load("/home/rstudio/results/risperidone-3q29/risWtSalWt-risDelsalDel.RData")
message("Loaded: risWtSalWt-risDelsalDel.RData")

# Define the DPI
dpi <- 100


# helper for null-safe fallback
`%||%` <- function(a, b) if (!is.null(a)) a else b

samples_salDel <- c("S13839Nr3", "S13839Nr6", "S13839Nr11", 
                    "S13839Nr14", "S13839Nr18", "S13839Nr26")


samples_risDel <- c("S13839Nr2", "S13839Nr7", "S13839Nr8", 
                    "S13839Nr9", "S13839Nr15", "S13839Nr17", 
                    "S13839Nr19", "S13839Nr21", "S13839Nr23")

# ==== Dataset registry ====
dataset_registry <- list(
  "Risperidone" = list(
    spatial_data = "risperidone_st_data_half",
    summary_data = "risperidone_summary_statistics_half",
    case_samples = "samples_risperidone",
    control_samples = "samples_saline"
  ),
  "PZ-1190" = list(
    spatial_data = "pz1190_st_data_half",
    summary_data = "pz1190_summary_statistics_half",
    case_samples = "samples_pz1190",
    control_samples = "samples_saline"
  ),
  "risWtSalWt" = list(
    spatial_data = "ris3q29_st_data",
    summary_data = "risWtSalWt_summary_statistics",
    case_samples = "samples_risWt",
    control_samples = "samples_salWt"
  ),
  "risDelSalDel" = list(
    spatial_data = "ris3q29_st_data",
    summary_data = "risDelSalDel_summary_statistics",
    case_samples = "samples_risDel",
    control_samples = "samples_salDel"
  )
)
