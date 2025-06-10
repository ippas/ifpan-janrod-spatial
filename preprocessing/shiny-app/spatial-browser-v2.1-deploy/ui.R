# -----------------------------
# Section: UI Layout
# -----------------------------

ui <- fluidPage(
  
  # -----------------------------
  # Section: Header
  # -----------------------------
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
      div(style = "text-align: right;", span(""))
    )
  ),
  
  tags$hr(),
  
  # -----------------------------
  # Section: Dataset Selector
  # -----------------------------
  fluidRow(
    column(4,
           selectInput("dataset", "Select Dataset", choices = names(dataset_registry))
    )
  ),
  
  tags$hr(),
  
  # -----------------------------
  # Section: Main Tabs
  # -----------------------------
  tabsetPanel(
    id = "main_tabs",
    
    # -----------------------------
    # Tab: Clustering Visualization
    # -----------------------------
    tabPanel("🔍 Clustering Visualization",
             sidebarLayout(
               sidebarPanel(
                 h4("Visualization Options"),
                 fluidRow(
                   column(6, numericInput("spot_size", "Spot Size", value = 1, min = 0.1, max = 5, step = 0.1)),
                   column(6, selectInput("num_columns", "Number of Columns", choices = 1:7, selected = 4))
                 ),
                 fluidRow(
                   column(6, numericInput("width_plot", "Width (px):", value = 1100, min = 100, step = 50)),
                   column(6, numericInput("height_plot", "Height (px):", value = 950, min = 100, step = 50))
                 ),
                 fluidRow(
                   column(6, numericInput("dpi_svg", "SVG DPI:", value = 100, min = 50, max = 600, step = 10))
                 ),
                 checkboxInput("tif_image", "Show Background Image", value = TRUE),
                 tags$hr(),
                 h4("Clustering Resolution"),
                 sliderInput("cluster_resolution_slider", "Select resolution:", min = 0.05, max = 2, value = 0.4, step = 0.05),
                 tags$hr(),
                 h4("Cluster Selection"),
                 div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                     actionButton("select_all_clusters", "Select All"),
                     actionButton("deselect_all_clusters", "Deselect All")
                 ),
                 fluidRow(
                   lapply(0:2, function(col_index) {
                     column(4,
                            checkboxGroupInput(
                              inputId = paste0("cluster_selector_col", col_index + 1),
                              label = NULL,
                              choices = paste0("cluster_", (col_index * 10):(col_index * 10 + 9)),
                              selected = paste0("cluster_", (col_index * 10):(col_index * 10 + 9))
                            )
                     )
                   })
                 ),
                 br(),
                 actionButton("update_resolution", "Update Resolution and Cluster"),
                 tags$hr(),
                 h4("Samples in Dataset"),
                 uiOutput("sample_display_block"),
                 tags$hr(),
                 downloadButton("downloadSpatialClusterPNG", "Download PNG"),
                 downloadButton("downloadSpatialClusterSVG", "Download SVG")
               ),
               mainPanel(
                 h4("Cluster Visualization"),
                 plotOutput("clust_plot", height = "600px")
               )
             )
    ),
    
    # -----------------------------
    # Tab: UMAP Visualization
    # -----------------------------
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
    
    # -----------------------------
    # Tab: Statistics Summary
    # -----------------------------
    tabPanel("📈 Statistics Summary",
             sidebarLayout(
               
               # ---- SIDEBAR ----
               sidebarPanel(
                 h4("Parameters for filtering statistics"),
                 
                 # METRIC + DATA TYPE + RESOLUTION
                 fluidRow(
                   column(4,
                          selectInput("stat_summary_metric", "Select metric:",
                                      choices = c("mean", "median", "sum"),
                                      selected = "mean")
                   ),
                   column(4,
                          selectInput("stat_data_type", "Select data type:",
                                      choices = c("raw_data", "seurat", "range_normalize", "quantile_metric", "quantile_normalize"),
                                      selected = "raw_data")
                   ),
                   column(4,
                          selectInput("stat_resolution", "Select resolution:",
                                      choices = c("0.05", "0.1", "0.15", "0.2", "0.4", "0.8", "1"),
                                      selected = "0.4")
                   )
                 ),
                 
                 tags$hr(),
                 h4("Thresholds"),
                 
                 # CONTROL + EXPERIMENT + LOG2(FC)
                 fluidRow(
                   column(4,
                          numericInput("control_mean_thresh", "Control Mean:", value = 0.2, min = 0, step = 0.01)
                   ),
                   column(4,
                          numericInput("experiment_mean_thresh", "Experiment Mean:", value = 0.2, min = 0, step = 0.01)
                   ),
                   column(4,
                          numericInput("log2ratio_thresh", "log2(FC):", value = 0.5, step = 0.01)
                   )
                 ),
                 
                 tags$hr(),
                 h4("Statistical test"),
                 
                 # TEST TYPE + P-VALUE
                 fluidRow(
                   column(6,
                          selectInput("stat_test_type", "Test:",
                                      choices = c("T-Student", "Wilcoxon", "Kolmogorov-Smirnov"),
                                      selected = "T-Student")
                   ),
                   column(6,
                          numericInput("pvalue_thresh", "P-value:", value = 0.05, min = 0, max = 1, step = 0.001)
                   )
                 ),
                 
                 tags$hr(),
                 h4("Download Table"),
                 downloadButton("download_stats_tsv", "Download TSV"),
                 downloadButton("download_stats_csv", "Download CSV")
               ),
               
               # ---- MAIN PANEL ----
               mainPanel(
                 h4("Filtered Statistics Table"),
                 DT::dataTableOutput("stats_table")
               )
             )
    ),
    
    # -----------------------------
    # Tab: Gene Expression (Spatial)
    # -----------------------------
    tabPanel("🧠 Gene Expression (Spatial)",
             sidebarLayout(
               sidebarPanel(
                 h4("Visualization Options"),
                 fluidRow(
                   column(6, numericInput("spot_size_expr", "Spot Size", value = 1, min = 0.1, max = 5, step = 0.1)),
                   column(6, selectInput("num_columns_expr", "Number of Columns", choices = 1:7, selected = 4))
                 ),
                 fluidRow(
                   column(6, numericInput("width_plot_expr", "Width (px):", value = 1300, min = 100, step = 50)),
                   column(6, numericInput("height_plot_expr", "Height (px):", value = 950, min = 100, step = 50))
                 ),
                 fluidRow(
                   column(6, numericInput("dpi_svg_expr", "SVG DPI:", value = 100, min = 50, max = 600, step = 10))
                 ),
                 checkboxInput("tif_image_expr", "Show Background Image", value = TRUE),
                 checkboxInput("zero_normalize_expr", "Zero Normalize", value = FALSE),
                 checkboxInput("show_legend_expr", "Show Legend", value = TRUE),
                 tags$hr(),
                 h4("Gene & Peak Selection"),
                 uiOutput("data_type_visualization_expr"),
                 uiOutput("gene_selector_expr"),
                 uiOutput("peak_selector_expr"),
                 tags$hr(),
                 h4("Percentile Filtering"),
                 fluidRow(
                   column(6, numericInput("min_percentile_expr", "Input min percentile to remove:", value = 0, min = 0, max = 1, step = 0.01)),
                   column(6, numericInput("max_percentile_expr", "Input max percentile to remove:", value = 1, min = 0, max = 1, step = 0.01))
                 ),
                 tags$hr(),
                 h4("Clustering Resolution"),
                 sliderInput("cluster_resolution_slider_expr", "Select resolution:", min = 0.05, max = 2, value = 0.4, step = 0.05),
                 tags$hr(),
                 h4("Cluster Selection"),
                 div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                     actionButton("select_all_clusters_expr", "Select All"),
                     actionButton("deselect_all_clusters_expr", "Deselect All")
                 ),
                 fluidRow(
                   lapply(0:2, function(col_index) {
                     column(4,
                            checkboxGroupInput(
                              inputId = paste0("cluster_selector_expr_col", col_index + 1),
                              label = NULL,
                              choices = paste0("cluster_", (col_index * 10):(col_index * 10 + 9)),
                              selected = paste0("cluster_", (col_index * 10):(col_index * 10 + 9))
                            )
                     )
                   })
                 ),
                 br(),
                 actionButton("update_resolution_expr", "Update Resolution and Cluster"),
                 tags$hr(),
                 h4("Samples in Dataset"),
                 uiOutput("sample_display_block_expr"),
                 tags$hr(),
                 downloadButton("downloadSpatialClusterPNG_expr", "Download PNG"),
                 downloadButton("downloadSpatialClusterSVG_expr", "Download SVG")
               ),
               mainPanel(
                 h4("Spatial Gene Expression Plot"),
                 plotOutput("expr_spatial_plot_cluster", height = "600px")
               )
             )
    ),
    
    # -----------------------------
    # Tab: Gene Expression (Spatial, Aggregate Cluster)
    # -----------------------------
    tabPanel("🧠 Gene Expression (Spatial, Aggregate Cluster)",
             sidebarLayout(
               sidebarPanel(
                 h4("Visualization Options"),
                 
                 # Spot size and columns
                 fluidRow(
                   column(6, numericInput("spot_size_agg", "Spot Size", value = 1, min = 0.1, max = 5, step = 0.1)),
                   column(6, selectInput("num_columns_agg", "Number of Columns", choices = 1:7, selected = 4))
                 ),
                 
                 # Plot dimensions
                 fluidRow(
                   column(6, numericInput("width_plot_agg", "Width (px):", value = 1300, min = 100, step = 50)),
                   column(6, numericInput("height_plot_agg", "Height (px):", value = 950, min = 100, step = 50))
                 ),
                 
                 # DPI
                 fluidRow(
                   column(6, numericInput("dpi_svg_agg", "SVG DPI:", value = 100, min = 50, max = 600, step = 10))
                 ),
                 
                 # Toggles
                 checkboxInput("tif_image_agg", "Show Background Image", value = TRUE),
                 checkboxInput("zero_normalize_agg", "Zero Normalize", value = FALSE),
                 checkboxInput("show_legend_agg", "Show Legend", value = TRUE),
                 
                 tags$hr(),
                 h4("Gene & Peak Selection"),
                 uiOutput("data_type_visualization_agg"),
                 uiOutput("gene_selector_agg"),
                 uiOutput("peak_selector_agg"),
                 
                 tags$hr(),
                 h4("Clustering Resolution"),
                 selectInput("cluster_resolution_agg", "Select resolution:",
                             choices = c(
                               "cluster_resolution_0.05",
                               "cluster_resolution_0.1",
                               "cluster_resolution_0.15",
                               "cluster_resolution_0.2",
                               "cluster_resolution_0.4",
                               "cluster_resolution_0.8"
                             ),
                             selected = "cluster_resolution_0.4"),
                 
                 tags$hr(),
                 h4("Summary Metric"),
                 selectInput("summary_metric_agg", "Metric:",
                             choices = c("mean", "median", "sum"),
                             selected = "mean"),
                 
                 tags$hr(),
                 h4("Samples in Dataset"),
                 uiOutput("sample_display_block_agg"),
                 
                 tags$hr(),
                 downloadButton("downloadSpatialClusterPNG_agg", "Download PNG"),
                 downloadButton("downloadSpatialClusterSVG_agg", "Download SVG")
               ),
               mainPanel(
                 h4("Aggregated Spatial Expression"),
                 plotOutput("spatial_agg_plot", height = "600px")
               )
             )
    ),
    
    # -----------------------------
    # Tab: Gene Expression (UMAP)
    # -----------------------------
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
    
    # -----------------------------
    # Tab: Gene Expression (Barplot)
    # -----------------------------
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
