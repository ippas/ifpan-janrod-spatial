# # Load the libraries
# install.packages("shinydashboard")
# install.packages("shiny")


source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")
require(shiny)
require(shinydashboard)

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))

# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]


load("results/risperidone/risperidone.RData")
load("results/risperidone/pz1190.RData")

# Define the user interface
ui <- fluidPage(
  titlePanel("Gene and Peak Selection"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("spatial_data",
                  "Select spatial data",
                  choices = c("risperidone_st_data_half",
                              "pz1190_st_data_half")),
      
      # uiOutput("spatial_data"),
      
      HTML("<h5><b> Parameters for choosing data </b></h5>"), 
      
      selectInput("summary_data",
                  "Select summary data",
                  choices = c("risperidone_summary_statistics_half",
                              "pz1190_summary_statistics_half")),
      
      selectInput("experiment_samples",
                  "Select experiment samples: ",
                  choices = c("samples_risperidone", "samples_pz1190")),
      
      selectInput("metric",
                  "Select metric",
                  choices = c("mean", 
                              "median", "sum", "IQR", "diff_range", "var", "skewness", "kurtosis")),
      
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
                               "0.8")),
      
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
               actionButton("submit_genes", "Submit Gene Names")
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
        column(6, numericInput("height_plot", "Height:", 750, min = 100, step = 10))
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
      
      sliderInput("cluster", "Select number cluster:", min = 0, max = 25, value = 0),
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
                   textInput("filename", label = NULL, value = paste("results_", Sys.Date(), ".csv", sep = ""))
               )
        ),
        column(4, downloadButton('downloadData', 'Download CSV'))
      ),
      
      verbatimTextOutput("selected_info"),
      
      plotOutput("peakPlot", height = "100%", width = "100%"),  
      
      plotOutput("clusterPlot", height = "100%", width = "100%"),
      
      plotOutput("interestClusterPlot", height = "100%", width = "100%"),
      
      plotOutput("plot_umis", height = "100%", width = "100%")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  vals <- reactiveValues(cluster =  0, resolution = 0.05)
  
  observeEvent(input$update_values, {
    vals$cluster <- input$cluster
    vals$resolution <- input$resolution
  })
  
  
  spatial_annotate <- reactive({get(input$spatial_data)[[input$data_type_visualization]]$annotate})
  
  data_type_vector <-  reactive({get(input$spatial_data) %>%
      names() %>% 
      .[grepl(paste(c("raw_data", "range_normalize", "quantile_normalize", "seurat"), collapse = "|"), .)]})
  
  # 1.5
  observeEvent(input$spatial_data, {
    output$data_type_visualization <- renderUI({
      selectInput("data_type_visualization", 
                  "Select Data Type:", 
                  choices = data_type_vector())
    })
    
    output$gene <- renderUI({
      selectInput("gene", 
                  "Select Gene:", 
                  choices = spatial_annotate() %>% .$gene_name)
    })
  })
  
  # Render a DataTable in the Shiny app
  output$filter_statistics <- renderDataTable({
    
    # The 'filter_data_statistics' function is called with inputs from the Shiny app's UI
    # This function seems to process/ filter the summary data based on several conditions
    # The filtered data is stored in the 'filter_data' object
    filter_data_statistics(summary_data = get(input$summary_data), 
                           data_type = input$data_type_statistics, 
                           resolution = input$resolution_statistics,
                           metric = input$metric,
                           control_mean_threshold = input$control_mean_threshold,
                           experiment_mean_threshold = input$experiment_mean_threshold,
                           log2ratio_threshold = input$log2ratio_threshold,
                           t_test_threshold = input$t_test_threshold,
                           wilcoxon_test_threshold = input$wilcoxon_test_threshold,
                           ks_test_threshold = input$ks_test_threshold) -> filter_data
    
    # The 'filter_data' object is passed through a pipeline:
    # - Numeric columns are rounded to 4 decimal places with the 'mutate_if' function
    # - The column "condition" is removed from the data with the 'select' function
    datatable_data <- filter_data %>% 
      mutate_if(is.numeric, round, 4) %>% 
      select(-condition)
    

    # This goes within your server function
    output$downloadData <- downloadHandler(
      filename = function() {
        # Use the filename inputted by the user. 
        # If you still want to add the current date, you can do that too.
        input$filename
      },
      content = function(file) {
        # Assuming 'datatable_data' contains the data you want to export
        write.csv(datatable_data, file, row.names = FALSE)
      }
    )
    
    # The DataTable is created with the 'datatable' function
    # It is configured with options like paging, scrolling, export buttons, and so on
    # The 'formatStyle' function is used to prevent the content in the "peak" column from wrapping
    datatable(datatable_data,
              options = list(paging = TRUE,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             # autoWidth = TRUE,
                             server = FALSE,
                             pageLength = 20,
                             buttons = c('csv', 'excel'),
                             columnDefs = list(list(targets = '_all', className = 'dt-center'))
                             # list(targets = c(0, 8, 9), visible = FALSE))
              ),
              extensions = 'Buttons',
              filter = 'top',
              selection = 'single',
              rownames = FALSE
    ) %>% formatStyle("peak","white-space"="nowrap")

  })
  
  # This goes within your server function
  observeEvent(input$submit_genes, {
    # Split the text from the text area into a vector of genes
    # Splitting by newline characters
    gene_vector <- strsplit(input$gene_names, split = "\n")[[1]]
    
    # Remove leading/trailing whitespaces from the gene names
    gene_vector <- trimws(gene_vector)
    
    # Remove leading/trailing special characters from the gene names
    # This uses regex to remove non-alphanumeric characters at the start and end of each gene name
    gene_vector <- gsub("^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$", "", gene_vector)
    
    # Remove empty elements from the vector
    gene_vector <- gene_vector[nzchar(gene_vector)]
    
    # Convert elements to lowercase
    # gene_vector <- tolower(gene_vector)
    
    # Print the vector of genes to the console
    print(gene_vector)
    
    # Filter the data based on other criteria
    filter_data_statistics(summary_data = get(input$summary_data), 
                           data_type = input$data_type_statistics, 
                           resolution = input$resolution_statistics,
                           metric = input$metric,
                           control_mean_threshold = input$control_mean_threshold,
                           experiment_mean_threshold = input$experiment_mean_threshold,
                           log2ratio_threshold = input$log2ratio_threshold,
                           t_test_threshold = input$t_test_threshold,
                           wilcoxon_test_threshold = input$wilcoxon_test_threshold,
                           ks_test_threshold = input$ks_test_threshold) -> filter_data
    
    filter_data <- filter_data[filter_data$gene %in% gene_vector, ]
    
    print(filter_data)
    
  })
  
  ### code to perform visualization
  
  # Observe when the gene selection changes and update the peak dropdown accordingly
  observeEvent(input$gene, {
    df <- spatial_annotate()
    # Filter the data based on the selected gene
    filtered_df <- df[df$gene_name == input$gene,]
    
    # Update the choices in the peak dropdown
    updateSelectInput(session, "peak", choices = unique(filtered_df$peak_id))
  })
  
  # Show the selected gene and peak
  output$selected_info <- renderPrint({
    paste("Selected gene:", input$gene, "\nSelected peak:", input$peak)
  })
  
  # Add a plot for the selected peak
  output$peakPlot <- renderPlot({
    spatial_feature_plot(spatial_data = get(input$spatial_data),
                         type_data = {{input$data_type_visualization}},
                         peak_id = {{input$peak}},
                         # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                         samples = c(samples_saline, get(input$experiment_samples)),
                         min_percentile = input$min_percentile,
                         max_percentile = input$max_percentile,
                         normalization = {{input$zero_normalize}},
                         size = input$spot_size,
                         tif_image = input$tif_image) +
      plot_layout(ncol = as.numeric({{input$num_columns}}))
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) })
  
  output$clusterPlot <- renderPlot({
    # visualize interest cluster
    spatial_cluster(spatial_data = get(input$spatial_data),
                    resolution = vals$resolution,
                    samples = c(samples_saline, get(input$experiment_samples)),
                    palette = palette_allen, 
                    size= 1.2, 
                    ncol = 4)
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) + 200})
  
  output$interestClusterPlot <- renderPlot({
    spatial_interest_cluster(cluster = vals$cluster,
                             spatial_data = get(input$spatial_data),
                             resolution = vals$resolution,
                             samples = c(samples_saline, get(input$experiment_samples)),
                             size= 1.3,
                             ncol = 4)
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) + 200})
  
  output$plot_umis <- renderPlot({
    plot_umi_reads_cluster(spatial_data = get(input$spatial_data), 
                           data_type = "raw_data",
                           resolution = {{vals$resolution}},
                           cluster = vals$cluster) -> p1
    
    plot_umi_reads_cluster(spatial_data = get(input$spatial_data),
                           data_type = {{input$data_type_visualization}},
                           resolution = {{vals$resolution}},
                           cluster = vals$cluster) -> p2
    
    # Combine the plots using the patchwork package
    p1/p2
  }, height = 800)
  
}


  # Run the app
shinyApp(ui = ui, server = server)

