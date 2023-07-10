
ui <- fluidPage(
  # Application title
  titlePanel("Spatial Clustering App"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Inputs
    sidebarPanel(
      numericInput("resolution", "Resolution:", value = 0.2, min = 0),
      numericInput("size", "Size:", value = 1.0, min = 0),
      numericInput("ncol", "Number of Columns:", value = 4, min = 1),
      sliderInput("nfeatures", "Number of Features:", min = 500, max = 3000, value = 2500, step = 250),
      sliderInput("dims", "Dimensions:", min = 10, max = 60, value = 40, step = 5)
    ),
    
    # Output
    mainPanel(
      textOutput("numClusters"),
      plotOutput("clusterPlot", width = "1200px", height = "1000px")
    )
  )
)

server <- function(input, output) {
  
  # You may need to define samples_saline, samples_risperidone, 
  # palette_allen and other relevant variables here or read them from files
  
  clusters_data <- reactive({
    # Retrieve clusters based on user-selected nfeatures and dims.
    nfeatures_dims <- paste0("nfeatures", input$nfeatures, "_", "dims", input$dims)
    clusters <- cluster_results_ris[[nfeatures_dims]]$clusters
    
    # Return the clusters
    clusters
  })
  
  output$numClusters <- renderText({
    # Calculate the number of clusters
    num_clusters <- clusters_data()
    
    resolution_name <- paste0("cluster_resolution_", input$resolution)
    
    num_clusters[[resolution_name]] %>% unique() %>% length() -> num_clusters
    
    paste("Number of clusters:", num_clusters)
  })
  
  output$clusterPlot <- renderPlot({
    # Assign clusters to tmp (assuming that's how your function expects the data)
    # tmp <- list()
    tmp$clusters <- clusters_data()
    
    # Executing the spatial_cluster function
    spatial_cluster(spatial_data = tmp,
                    resolution = input$resolution,
                    samples = c(samples_saline, samples_risperidone),
                    palette = palette_allen,
                    size = input$size,
                    ncol = input$ncol)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
