

# Jeśli uruchamiasz z wyższego katalogu – ustaw katalog roboczy:
setwd("preprocessing/shiny-app/spatial-browser-v2.0-dev")

ui <- fluidPage(
  
  # ===== HEADER WITH PLOTS =====
  fluidRow(
    column(2, plotOutput("logo_left", height = "80px")),
    column(8, div(style = "text-align: center; margin-top: 15px;",
                  tags$h2("Spatial Data Browser"))),
    column(2, plotOutput("logo_right", height = "80px"))
  ),
  
  tags$hr(),
  
  # ===== DATASET SELECTOR =====
  fluidRow(
    column(4, selectInput("dataset", "Select Dataset:", choices = c("dataset1", "dataset2")))
  ),
  
  tags$hr(),
  
  # ===== TABS =====
  tabsetPanel(
    id = "main_tabs",
    
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
    
    tabPanel("🔘 UMAP Visualization",
             sidebarLayout(
               sidebarPanel(
                 h4("UMAP Settings"),
                 selectInput("color_by", "Color By:", choices = c("Cluster", "Sample", "Gene Expression"))
               ),
               mainPanel(
                 h4("UMAP Plot"),
                 plotOutput("umap_plot", height = "400px")
               )
             )
    ),
    
    tabPanel("📈 Statistics Summary",
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
    
    tabPanel("🧠 Gene Expression (Spatial)",
             sidebarLayout(
               sidebarPanel(
                 h4("Gene Selection"),
                 selectInput("gene", "Select Gene:", choices = c("Gene1", "Gene2", "Gene3")),
                 selectInput("display_mode", "Display Mode:", choices = c("Spatial", "Violin", "Boxplot"))
               ),
               mainPanel(
                 h4("Spatial Expression Plot"),
                 plotOutput("expr_spatial_plot", height = "400px")
               )
             )
    ),
    
    tabPanel("🧠 Gene Expression (Spatial, Aggregate Cluster)",
             sidebarLayout(
               sidebarPanel(
                 h4("Aggregate Settings"),
                 selectInput("gene_agg", "Select Gene:", choices = c("Gene1", "Gene2", "Gene3")),
                 selectInput("cluster_agg", "Cluster Level:", choices = c("cluster_1", "cluster_2", "cluster_3"))
               ),
               mainPanel(
                 h4("Aggregated Spatial Expression"),
                 plotOutput("spatial_agg_plot", height = "400px")
               )
             )
    ),
    
    tabPanel("🔘 Gene Expression (UMAP)",
             sidebarLayout(
               sidebarPanel(
                 h4("Gene Selection"),
                 selectInput("gene_umap", "Select Gene:", choices = c("Gene1", "Gene2", "Gene3"))
               ),
               mainPanel(
                 h4("UMAP Colored by Gene Expression"),
                 plotOutput("umap_expr_plot", height = "400px")
               )
             )
    ),
    
    tabPanel("🧬 Gene Expression (Barplot)",
             sidebarLayout(
               sidebarPanel(
                 h4("Gene Selection"),
                 selectInput("gene_barplot", "Select Gene:", choices = c("Gene1", "Gene2", "Gene3"))
               ),
               mainPanel(
                 h4("Barplot of Expression"),
                 plotOutput("barplot_expr", height = "400px")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # ===== LOGO PLOTS =====
  output$logo_left <- renderPlot({
    img <- readPNG("spatial-browser-logo.png")
    grid::grid.raster(img)
  })
  
  output$logo_right <- renderPlot({
    img <- readPNG("ifpas-logo.png")
    grid::grid.raster(img)
  })
  
  # ===== MOCKUP CONTENT =====
  output$clust_plot <- renderPlot({
    plot(1:10, main = "Clustering Preview (mockup)")
  })
  
  output$umap_plot <- renderPlot({
    plot(runif(100), runif(100), col = "gray", pch = 16, main = "UMAP Plot")
  })
  
  output$stats_table <- renderTable({
    data.frame(
      Gene = c("GeneA", "GeneB"),
      log2Ratio = c(1.2, -0.8),
      p_value = c(0.01, 0.03)
    )
  })
  
  output$expr_spatial_plot <- renderPlot({
    hist(rnorm(100), main = paste("Spatial Expression:", input$gene))
  })
  
  output$spatial_agg_plot <- renderPlot({
    image(matrix(runif(100), 10), main = paste("Aggregate:", input$gene_agg, input$cluster_agg))
  })
  
  output$umap_expr_plot <- renderPlot({
    plot(runif(100), runif(100), col = rainbow(100), pch = 16,
         main = paste("UMAP Expression:", input$gene_umap))
  })
  
  output$barplot_expr <- renderPlot({
    barplot(rnorm(5), names.arg = LETTERS[1:5],
            main = paste("Expression:", input$gene_barplot))
  })
}

shinyApp(ui, server)
