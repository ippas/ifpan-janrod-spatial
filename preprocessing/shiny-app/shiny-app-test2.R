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
      selectInput("gene", 
                  "Select Gene:", 
                  choices = unique(df$gene_name)),
      
      selectInput("peak", 
                  "Select Peak:", 
                  choices = NULL),
      
      radioButtons("type_data", 
                   "Choose data type:", 
                   choices = c("raw_data", "quantile_normalize")),
      
      checkboxInput("zero_normalize", 
                    "Zero Normalize", 
                    value = FALSE),
      
      selectInput("num_columns",
                  "Number of Columns to Display:",
                  choices = c(2, 4)),
      
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

