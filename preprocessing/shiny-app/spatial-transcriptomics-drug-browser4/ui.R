# Define the user interface
ui <- fluidPage(
  titlePanel("Gene and Peak Selection"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("drug", "Select Drug", choices = c("Risperidone", "PZ-1190", "Clozapine", "risWtSalWt", "risDelSalDel")),
      
      selectInput("spatial_data", "Select spatial data", choices = NULL),
      
      HTML("<h5><b>Parameters for choosing data</b></h5>"),
      
      selectInput("summary_data", "Select summary data", choices = NULL),
      
      selectInput("experiment_samples", "Select experiment samples", choices = NULL),
      
      selectInput("metric",
                  "Select metric",
                  choices = c("mean", 
                              "median", "sum")),
      
      radioButtons("data_type_statistics", 
                   "Choose data type for statistics:", 
                   choices = c("raw_data", 
                               "seurat",
                               "range_normalize",
                               "quantile_metric",
                               "quantile_normalize")),
      
      radioButtons("resolution_statistics", 
                   "Choose resolution for statistics:", 
                   choices = c("0.05",
                               "0.1",
                               "0.15",
                               "0.2",
                               "0.4",
                               "0.8", 
                               "1")),
      
      HTML("<h5><b> Parameters for filtering statistics </b></h5>"), 
      
      # Custom CSS to restrict resizing to vertical only
      tags$head(tags$style(HTML("#gene_names { resize: vertical;   height: 355px;  }"))),
      
      # UI layout
      fluidRow(
        # Left column with parameters
        column(6,
               numericInput("control_mean_threshold", "Control Mean Threshold:", value = 0.2, step = 0.1),
               numericInput("experiment_mean_threshold", "Experiment Mean Threshold:", value = 0.2, step = 0.1),
               numericInput("log2ratio_threshold", "Log2Ratio Threshold:", value = 0.5, step = 0.1),
               numericInput("t_test_threshold", "T-Student Test Threshold:", value = 0.05, step = 0.1),
               numericInput("wilcoxon_test_threshold", "Wilcoxon Test Threshold:", value = 0.05, step = 0.1),
               numericInput("ks_test_threshold", "Kolmogorov-Smirnov Test Threshold:", value = 0.05, step = 0.1)
        ),
        
        # Right column with gene names input and submit button
        column(6,
               textAreaInput("gene_names", "Enter Gene Names:", rows = 5),
               actionButton("clear_genes", "Clear Gene Names", rows = 5)
               # actionButton("submit_genes", "Submit Gene Names")
        )
      ),
      
      HTML("<h5><b> Visualization parameters </b></h5>"), 
      
      uiOutput("data_type_visualization"),
      
      checkboxInput("zero_normalize", 
                    "Zero Normalize", 
                    value = FALSE),
      
      checkboxInput("tif_image", 
                    "Background image", 
                    value = TRUE),
      
      checkboxInput("show_legend", "Show Legend", value = FALSE),
      
      fluidRow(
        column(6, uiOutput("gene")),
        column(6, selectInput("peak", "Select Peak:", choices = NULL))
      ),
      
      fluidRow(
        column(6, selectInput("num_columns", "Number of Columns to Display:", choices = c(4, 2))),
        column(6, numericInput("spot_size", "Spot size:", 1, min = 0.1, max = 5, step = 0.1)),
      ),
      
      fluidRow(
        column(6, numericInput("width_plot", "Width:", 1100, min = 100, step = 10)),
        column(6, numericInput("height_plot", "Height:", 950, min = 100, step = 10))
      ),
      
      fluidRow(
        column(6, 
               numericInput("min_percentile", 
                            "Input min percentile to remove:", 
                            value = 0, min = 0, max = 1, step = 0.001)
        ),
        column(6, 
               numericInput("max_percentile", 
                            "Input max percentile to remove 2:", 
                            value = 1, min = 0, max = 1, step = 0.001)
        )
      ),
      
      fluidRow(
        column(6,
               h5(tags$b("Saline")),
               lapply(c("S6269Nr1", "S6269Nr3", "S7788Nr1", "S7788Nr2", "S7788Nr3", "S7788Nr11"),
                      function(x) checkboxInput(inputId = paste("saline_", gsub("Nr|S", "", x), sep = ""),
                                                label = x, value = TRUE))),
        column(6,
               h5(tags$b("Experimental")),
               lapply(c("S6230Nr3", "S6230Nr4", "S6269Nr2", "S6269Nr4", "S7788Nr15", "S7788Nr16"),
                      function(x) checkboxInput(inputId = paste("experimental_", gsub("Nr|S", "", x), sep = ""),
                                                label = x, value = TRUE))),
        column(6,
               h5(tags$b("salWt")),
               lapply(
                 # c("S5023Nr1", "S5023Nr4", "S5295Nr1", "S5295Nr6", "S5295Nr7"),
                 c("S13839Nr1", "S13839Nr5", "S13839Nr12", "S13839Nr20", "S13839Nr22"),
                      function(x) checkboxInput(inputId = paste("saline_", gsub("Nr|S", "", x), sep = ""),
                                                label = x, value = F))),
        column(6,
               h5(tags$b("risWt")),
               lapply(
                 # c("S5023Nr2", "S5023Nr3", "S5295Nr3", "S5295Nr8"),
                 c("S13839Nr4", "S13839Nr10", "S13839Nr13", "S13839Nr16", "S13839Nr24", "S13839Nr25"),
                      function(x) checkboxInput(inputId = paste("experimental_", gsub("Nr|S", "", x), sep = ""),
                                                label = x, value = F))),
        column(6,
               h5(tags$b("salDel")),
               lapply(
                 c("S13839Nr3",  "S13839Nr6",  "S13839Nr11", "S13839Nr14", "S13839Nr18", "S13839Nr26"),
                 function(x) checkboxInput(inputId = paste("saline_", gsub("Nr|S", "", x), sep = ""),
                                           label = x, value = F))),
        column(6,
               h5(tags$b("risDel")),
               lapply(
                 c("S13839Nr2", "S13839Nr7", "S13839Nr8", "S13839Nr9", "S13839Nr15", "S13839Nr17", "S13839Nr19", "S13839Nr21", "S13839Nr23"),
                 function(x) checkboxInput(inputId = paste("experimental_", gsub("Nr|S", "", x), sep = ""),
                                           label = x, value = F)))
        
      ),
      
      actionButton("update_samples", "Update Samples"),
      
      HTML("<h5><b> Select clusters </b></h5>"), 
      
      # Dynamically create columns and distribute checkboxes
      lapply(1:3, function(col) {
        column(4,
               lapply(0:32, function(i) {
                 if (i %% 3 == col - 1) {
                   checkboxInput(inputId = paste("cluster", i, sep = ""), 
                                 label = paste("Cluster", i), 
                                 value = TRUE)
                 }
               })
        )
      }),
      
      # sliderInput("cluster", "Select number cluster:", min = 0, max = 25, value = 0),
      sliderInput("resolution", "Select resolution:", min = 0.05, max = 2, step = 0.05, value = 0.05),
      actionButton("update_values", "Update Resolution and Cluster")
      
    ),
    
    mainPanel(
      
      # tableOutput("filter_statistics"),
      dataTableOutput("filter_statistics"),
      
      # # Add this to your UI layout
      fluidRow(
        column(4, 
               div(style="display: inline-flex; align-items: center;", 
                   tags$span("Enter Filename:", style="font-weight: bold; margin-right: 10px;"),
                   textInput("filename", label = NULL, value = paste("results_", Sys.Date(), ".tsv", sep = ""))
               )
        ),
        column(4, downloadButton('downloadData', 'Download TSV'))
      ),
      
      verbatimTextOutput("selected_info"),
      
      plotOutput("peakPlot", height = "100%", width = "100%"),  
      
      downloadButton('downloadSpatialFeaturePlotPNG', 'Download PNG'),
      
      downloadButton('downloadSpatialFeaturePlotSVG', 'Download SVG'),
      
      plotOutput("clusterPlot", height = "100%", width = "100%"),
      downloadButton('downloadSpatialClusterPNG', 'Download PNG'),
      downloadButton('downloadSpatialClusterSVG', 'Download SVG'),
      
      # plotOutput("interestClusterPlot", height = "100%", width = "100%"),
      # downloadButton('downloadInterestClusterPNG', 'Download PNG'),
      # downloadButton('downloadInterestClusterSVG', 'Download SVG'),
      
      # plotOutput("plot_umis", height = "100%", width = "100%")
    )
  )
)
