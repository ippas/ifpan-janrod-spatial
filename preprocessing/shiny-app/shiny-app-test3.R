# Load the libraries
library(shiny)

# Load your data
# Replace this with the actual path to your data
df <-risperidone_st_data_half$raw_data$annotate


# Define the user interface
ui <- fluidPage(
  titlePanel("Gene and Peak Selection"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("spatial_data",
                  "Select spatial data",
                  choices = c("risperidone_st_data_half")),
      
      HTML("<h5><b> Parameters for choosing data </b></h5>"), 
      
      selectInput("summary_data",
                  "Select summary data",
                  choices = c("risperidone_summary_statistics", 
                              "risperidone_summary_statistics_remove_ris")),
      
      radioButtons("data_type_statistics", 
                   "Choose data type for statistics:", 
                   choices = c("raw_data", 
                               "range_normalize",
                               "quantile_normalize")),
      
      radioButtons("resolution_statistics", 
                   "Choose resolution for statistics:", 
                   choices = c("0.1", 
                               "0.4",
                               "0.8")),
      
      HTML("<h5><b> Parameters for filtering statistics </b></h5>"), 
      
      # Numeric input field in the same line as text
      fluidRow(
        column(6, numericInput("control_mean_threshold", "Control mean:", value = 0, min = 0, step = 0.1)),
        column(6, numericInput("experiment_mean_threshold", "Experiment mean:", value = 0, min = 0, step = 0.1))
      ),
      
      fluidRow(
        column(6, numericInput("log2ratio_threshold", "log2ratio:", value = 0.5, min = 0, step = 0.1))
      ),
      
      fluidRow(
        column(6, numericInput("t_test_threshold", "t-Student test:", value = 0.05, min = 0, max = 1, step = 0.01))
      ),
      
      fluidRow(
        column(6, numericInput("wilcoxon_test_threshold", "Wilcoxon test:", value = 0.05, min = 0, max = 1, step = 0.01))
      ),
      
      fluidRow(
        column(6, numericInput("ks_test_threshold", "Kolmogorov-Smirnov test:", value = 0.05, min = 0, max = 1, step = 0.01))
      ),
      
      HTML("<h5><b> Visualization parameters </b></h5>"), 
      
      selectInput("gene", 
                  "Select Gene:", 
                  choices = unique(df$gene_name)),
      
      selectInput("peak", 
                  "Select Peak:", 
                  choices = NULL),
      
      radioButtons("type_data", 
                   "Choose data type:", 
                   choices = c("raw_data", "quantile_normalize")),
      
      selectInput("type_data_visualization",
                  "Select data type to visualization:",
                  choices = ),
      
      checkboxInput("zero_normalize", 
                    "Zero Normalize", 
                    value = FALSE),
      
      selectInput("num_columns",
                  "Number of Columns to Display:",
                  choices = c(4, 2)),
      
      numericInput("min_percentile", 
                   "Input min percentile to remove:", 
                   value = 0, min = 0, max = 1, step = 0.001),
      
      numericInput("max_percentile", 
                   "Input max percentile to remove 2:", 
                   value = 1, min = 0, max = 1, step = 0.001),
      
      sliderInput("cluster", "Select number cluster:", min = 0, max = 25, value = 0),
      sliderInput("resolution", "Select resolution:", min = 0.1, max = 2, step = 0.1, value = 0.1)
      
      
    ),
    
    mainPanel(
      
      # tableOutput("filter_statistics"),
      dataTableOutput("filter_statistics"),
      
      verbatimTextOutput("selected_info"),
      plotOutput("peakPlot",
                 height = 1100),  # Add this to display a plot
      
      plotOutput("clusterPlot",
                 height = 1100),
      
      plotOutput("interestClusterPlot",
                 height = 1100)
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  spatial_annotate <- reactive({get(input$spatial_data)[[input$type_data]]$annotate})
  #   if(input$spatial_data == "risperidone_st_data_half") {
  #     df <- risperidone_st_data_half$raw_data$annotate
  #   } else if (input$data_source == "other_data") {
  #     df <- other_data  # make sure to define 'other_data'
  #   }
  #   df
  # })
  
  # Render a DataTable in the Shiny app
  output$filter_statistics <- renderDataTable({
    
    # The 'filter_data_statistics' function is called with inputs from the Shiny app's UI
    # This function seems to process/ filter the summary data based on several conditions
    # The filtered data is stored in the 'filter_data' object
    filter_data_statistics(summary_data = get(input$summary_data), 
                           data_type = input$data_type_statistics, 
                           resolution = input$resolution_statistics,
                           metric = "mean",
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
    
    # The DataTable is created with the 'datatable' function
    # It is configured with options like paging, scrolling, export buttons, and so on
    # The 'formatStyle' function is used to prevent the content in the "peak" column from wrapping
    datatable(datatable_data,
              options = list(paging = TRUE,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             # autoWidth = TRUE,
                             server = FALSE,
                             dom = 'Bfrtip',
                             pageLength = 20,
                             buttons = c('csv', 'excel'),
                             columnDefs = list(list(targets = '_all', className = 'dt-center'))
                             # list(targets = c(0, 8, 9), visible = FALSE))
              ),
              extensions = 'Buttons',
              filter = 'top',
              selection = 'single',
              rownames = FALSE
    )  %>% formatStyle("peak","white-space"="nowrap")
  }) 
  
  
  # Observe when the gene selection changes and update the peak dropdown accordingly
  observeEvent(input$gene, {
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
    spatial_feature_plot(spatial_data = risperidone_st_data_half,
                         type_data = {{input$type_data}},
                         # type_data = "quantile_normalize",
                         peak_id = {{input$peak}},
                         # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                         samples = c(samples_saline, samples_risperidone),
                         min_percentile = input$min_percentile,
                         max_percentile = input$max_percentile,
                         normalization = {{input$zero_normalize}}) +
      plot_layout(ncol = as.numeric({{input$num_columns}}))
  })
  
  output$clusterPlot <- renderPlot({
    # visualize interest cluster
    spatial_cluster(spatial_data = risperidone_st_data_half,
                    resolution = input$resolution,
                    samples = c(samples_saline[-1], samples_risperidone[-3]),
                    palette = palette_allen, 
                    size= 1.2, 
                    ncol = 4)
  })
  
  output$interestClusterPlot <- renderPlot({
    # visualize interest cluster
    spatial_interest_cluster(cluster = input$cluster,
                             # seurat_object = integrated_analysis,
                             spatial_data = risperidone_st_data_half,
                             resolution = input$resolution,
                             samples = c(samples_saline, samples_risperidone),
                             size= 1.3,
                             ncol = 4)
  })
}

# Run the app
shinyApp(ui = ui, server = server)

