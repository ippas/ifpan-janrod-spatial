library(shiny)
library(png)
library(grid)

# SET WD if needed
setwd("preprocessing/shiny-app/spatial-browser-v2.0-dev")

# ==== Dataset registry ====
dataset_registry <- list(
  "Risperidone" = list(
    spatial_data = "risperidone_st_data_half",
    summary_data = "risperidone_summary_statistics_half",
    case_samples = "samples_risperidone",
    control_samples = "samples_control_risperidone"
  ),
  "PZ-1190" = list(
    spatial_data = "pz1190_st_data_half",
    summary_data = "pz1190_summary_statistics_half",
    case_samples = "samples_pz1190",
    control_samples = "samples_control_pz1190"
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

# ==== UI ====
ui <- fluidPage(
  
  # ==== HEADER ====
  fluidRow(
    column(
      6,
      fluidRow(
        column(2, plotOutput("logo_left", height = "70px")),
        column(10, div(style = "margin-top: 15px;", tags$h2("Spatial Data Browser")))
      )
    ),
    column(
      6,
      div(style = "text-align: right;",
          span(""))  # empty for now, logo IFPAN removed
    )
  ),
  
  tags$hr(),
  
  # ==== DATASET SELECTOR ====
  fluidRow(
    column(4,
           selectInput("dataset", "Select Dataset", choices = names(dataset_registry)))
  ),
  
  tags$hr(),
  
  # ==== TABS ====
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

# ==== SERVER ====
server <- function(input, output, session) {
  
  # === REACTIVE ACCESS TO SELECTED DATASET CONFIG ===
  dataset_config <- reactive({
    req(input$dataset)
    dataset_registry[[input$dataset]]
  })
  
  # === LOGO ===
  output$logo_left <- renderPlot({
    img <- readPNG("spatial-browser-logo.png")
    grid::grid.raster(img)
  })
  
  # === MOCKUP RENDERINGS ===
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

# ==== APP ====
shinyApp(ui, server)
