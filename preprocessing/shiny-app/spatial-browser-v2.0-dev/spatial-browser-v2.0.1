ui <- fluidPage(
  titlePanel("Spatial Data Browser"),
  
  # Global Data Selection Bar (at the top)
  fluidRow(
    column(3, selectInput("dataset", "Select Dataset:", choices = c("dataset1", "dataset2"))),
    column(3, selectInput("summary", "Select Summary Data:", choices = c("summary1", "summary2"))),
    column(3, selectInput("samples", "Select Sample Set:", choices = c("samples1", "samples2")))
  ),
  tags$hr(),
  
  # Tab Layout
  tabsetPanel(
    id = "main_tabs",
    
    # --- Tab 1: Clustering Visualization ---
    tabPanel("🔍 Clustering Visualization",
             sidebarLayout(
               sidebarPanel(
                 h4("Settings"),
                 selectInput("clust_method", "Clustering Method:", choices = c("K-means", "Hierarchical")),
                 selectInput("dim_reduction", "Dimension Reduction:", choices = c("UMAP", "PCA"))
               ),
               mainPanel(
                 h4("Clustering Result"),
                 plotOutput("clust_plot", height = "400px")
               )
             )
    ),
    
    # --- Tab 2: Statistics Summary ---
    tabPanel("📊 Statistics Summary",
             sidebarLayout(
               sidebarPanel(
                 h4("Statistics Options"),
                 selectInput("stat_metric", "Metric:", choices = c("mean", "median", "log2Ratio")),
                 sliderInput("pval_thresh", "P-value Threshold:", min = 0, max = 1, value = 0.05)
               ),
               mainPanel(
                 h4("Summary Table"),
                 tableOutput("stats_table")
               )
             )
    ),
    
    # --- Tab 3: Gene Expression ---
    tabPanel("🧬 Gene Expression",
             sidebarLayout(
               sidebarPanel(
                 h4("Gene Selection"),
                 selectInput("gene", "Select Gene:", choices = c("Gene1", "Gene2", "Gene3")),
                 selectInput("display_mode", "Display Mode:", choices = c("Spatial", "Violin", "Boxplot"))
               ),
               mainPanel(
                 h4("Expression Plot"),
                 plotOutput("expr_plot", height = "400px")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  # Placeholder plot for clustering
  output$clust_plot <- renderPlot({
    plot(1:10, main = "Clustering Preview (mockup)")
  })
  
  # Placeholder statistics table
  output$stats_table <- renderTable({
    data.frame(
      Gene = c("GeneA", "GeneB"),
      log2Ratio = c(1.2, -0.8),
      p_value = c(0.01, 0.03)
    )
  })
  
  # Placeholder expression plot
  output$expr_plot <- renderPlot({
    hist(rnorm(100), main = paste("Mock Expression:", input$gene))
  })
}

shinyApp(ui, server)
