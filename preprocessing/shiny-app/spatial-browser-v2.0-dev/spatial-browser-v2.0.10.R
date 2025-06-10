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
                 
                 selectInput("stat_summary_metric", "Select metric:",
                             choices = c("mean", "median", "sum"),
                             selected = "mean"),
                 
                 tags$hr(),
                 h4("Choose data type for statistics:"),
                 radioButtons("stat_data_type", label = NULL,
                              choices = c("raw_data", "seurat", "range_normalize", "quantile_metric", "quantile_normalize"),
                              selected = "raw_data"),
                 
                 tags$hr(),
                 h4("Choose resolution for statistics:"),
                 fluidRow(
                   column(6,
                          checkboxGroupInput("stat_resolution_col1", label = NULL,
                                             choices = c("0.05", "0.1", "0.15"),
                                             selected = "0.4")
                   ),
                   column(6,
                          checkboxGroupInput("stat_resolution_col2", label = NULL,
                                             choices = c("0.2", "0.4", "0.8", "1"),
                                             selected = "0.4")
                   )
                 ),
                 
                 tags$hr(),
                 h4("Thresholds"),
                 numericInput("control_mean_thresh", "Control Mean Threshold:", value = 0.2, min = 0, step = 0.01),
                 numericInput("experiment_mean_thresh", "Experiment Mean Threshold:", value = 0.2, min = 0, step = 0.01),
                 numericInput("log2ratio_thresh", "Log2Ratio Threshold:", value = 0.5, step = 0.01),
                 
                 tags$hr(),
                 h4("Statistical Test"),
                 radioButtons("stat_test_type", "Choose test:",
                              choices = c("T-Student", "Wilcoxon", "Kolmogorov-Smirnov"),
                              selected = "T-Student"),
                 
                 numericInput("pvalue_thresh", "P-value Threshold:", value = 0.05, min = 0, max = 1, step = 0.001),
                 
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
    tabPanel("💜 Gene Expression (Spatial)",
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

# ==== SERVER ====
server <- function(input, output, session) {
  
  # -----------------------------
  # Section: Global / Shared Reactives
  # -----------------------------
  
  vals <- reactiveValues(
    resolution = 0.5,
    select_clusters = as.character(0:29),
    resolution_expr = 0.4,
    select_clusters_expr = as.character(0:29)
  )
  
  dataset_config <- reactive({
    req(input$dataset)
    dataset_registry[[input$dataset]]
  })
  
  # -----------------------------
  # Section: Clustering Visualization
  # -----------------------------
  
  observeEvent(input$update_resolution, {
    vals$resolution <- input$cluster_resolution_slider
    selected_clusters <- unlist(c(
      input$cluster_selector_col1 %||% character(0),
      input$cluster_selector_col2 %||% character(0),
      input$cluster_selector_col3 %||% character(0)
    ))
    vals$select_clusters <- gsub("^cluster_", "", selected_clusters)
  })
  
  observeEvent(input$select_all_clusters, {
    lapply(0:2, function(col_index) {
      updateCheckboxGroupInput(
        session,
        inputId = paste0("cluster_selector_col", col_index + 1),
        selected = paste0("cluster_", (col_index * 10):(col_index * 10 + 9))
      )
    })
  })
  
  observeEvent(input$deselect_all_clusters, {
    lapply(1:3, function(i) {
      updateCheckboxGroupInput(session, paste0("cluster_selector_col", i), selected = character(0))
    })
  })
  
  observeEvent(input$select_all_control, {
    updateCheckboxGroupInput(session, "control_selected",
                             selected = get(dataset_config()$control_samples))
  })
  
  observeEvent(input$deselect_all_control, {
    updateCheckboxGroupInput(session, "control_selected", selected = character(0))
  })
  
  observeEvent(input$select_all_case, {
    updateCheckboxGroupInput(session, "case_selected",
                             selected = get(dataset_config()$case_samples))
  })
  
  observeEvent(input$deselect_all_case, {
    updateCheckboxGroupInput(session, "case_selected", selected = character(0))
  })
  
  output$sample_display_block <- renderUI({
    req(dataset_config())
    control <- get(dataset_config()$control_samples)
    case <- get(dataset_config()$case_samples)
    
    fluidRow(
      column(6,
             h5(strong("Control Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_control", "Select All"),
                 actionButton("deselect_all_control", "Deselect All")
             ),
             checkboxGroupInput("control_selected", NULL, choices = control, selected = control)
      ),
      column(6,
             h5(strong("Experimental Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_case", "Select All"),
                 actionButton("deselect_all_case", "Deselect All")
             ),
             checkboxGroupInput("case_selected", NULL, choices = case, selected = case)
      )
    )
  })
  
  output$logo_left <- renderPlot({
    img <- png::readPNG("spatial-browser-logo.png")
    grid::grid.raster(img)
  })
  
  output$clust_plot <- renderPlot({
    req(input$dataset)
    selected_samples <- c(input$control_selected, input$case_selected)
    req(length(selected_samples) > 0, vals$select_clusters, vals$resolution)
    
    spatial_cluster_select(
      spatial_data = get(dataset_config()$spatial_data),
      resolution = vals$resolution,
      samples = selected_samples,
      palette = palette_allen,
      size = input$spot_size,
      select_clusters = vals$select_clusters,
      tif_image = input$tif_image,
      ncol = as.numeric(input$num_columns)
    ) + patchwork::plot_layout(ncol = as.numeric(input$num_columns))
  },
  width = function() as.integer(input$width_plot),
  height = function() as.integer(input$height_plot))
  
  output$downloadSpatialClusterPNG <- downloadHandler(
    filename = function() {
      paste0("cluster-", input$dataset, "-res", vals$resolution, ".png")
    },
    content = function(file) {
      png(file, width = as.integer(input$width_plot), height = as.integer(input$height_plot))
      print(output$clust_plot())
      dev.off()
    }
  )
  
  output$downloadSpatialClusterSVG <- downloadHandler(
    filename = function() {
      paste0("cluster-", input$dataset, "-res", vals$resolution, ".svg")
    },
    content = function(file) {
      svg(file, width = as.integer(input$width_plot) / 100, height = as.integer(input$height_plot) / 100)
      print(output$clust_plot())
      dev.off()
    }
  )
  
  # -----------------------------
  # Section: Gene Expression (Spatial)
  # -----------------------------
  
  observeEvent(input$update_resolution_expr, {
    vals$resolution_expr <- input$cluster_resolution_slider_expr
    
    selected_clusters <- unlist(c(
      input$cluster_selector_expr_col1 %||% character(0),
      input$cluster_selector_expr_col2 %||% character(0),
      input$cluster_selector_expr_col3 %||% character(0)
    ))
    vals$select_clusters_expr <- gsub("^cluster_", "", selected_clusters)
  })
  
  observeEvent(input$select_all_clusters_expr, {
    lapply(0:2, function(col_index) {
      updateCheckboxGroupInput(
        session,
        inputId = paste0("cluster_selector_expr_col", col_index + 1),
        selected = paste0("cluster_", (col_index * 10):(col_index * 10 + 9))
      )
    })
  })
  
  observeEvent(input$deselect_all_clusters_expr, {
    lapply(1:3, function(i) {
      updateCheckboxGroupInput(session, paste0("cluster_selector_expr_col", i), selected = character(0))
    })
  })
  
  observeEvent(input$select_all_control_expr, {
    updateCheckboxGroupInput(session, "control_selected_expr",
                             selected = get(dataset_config()$control_samples))
  })
  
  observeEvent(input$deselect_all_control_expr, {
    updateCheckboxGroupInput(session, "control_selected_expr", selected = character(0))
  })
  
  observeEvent(input$select_all_case_expr, {
    updateCheckboxGroupInput(session, "case_selected_expr",
                             selected = get(dataset_config()$case_samples))
  })
  
  observeEvent(input$deselect_all_case_expr, {
    updateCheckboxGroupInput(session, "case_selected_expr", selected = character(0))
  })
  
  output$sample_display_block_expr <- renderUI({
    req(dataset_config())
    control <- get(dataset_config()$control_samples)
    case <- get(dataset_config()$case_samples)
    
    fluidRow(
      column(6,
             h5(strong("Control Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_control_expr", "Select All"),
                 actionButton("deselect_all_control_expr", "Deselect All")
             ),
             checkboxGroupInput("control_selected_expr", NULL, choices = control, selected = control)
      ),
      column(6,
             h5(strong("Experimental Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_case_expr", "Select All"),
                 actionButton("deselect_all_case_expr", "Deselect All")
             ),
             checkboxGroupInput("case_selected_expr", NULL, choices = case, selected = case)
      )
    )
  })
  
  # Returns available annotation table based on spatial_data and selected data type
  spatial_annotate_expr <- reactive({
    req(dataset_config()$spatial_data, input$data_type_visualization_expr)
    get(dataset_config()$spatial_data)[[input$data_type_visualization_expr]]$annotate
  })
  
  # Returns available data types (e.g., raw, normalized, seurat)
  data_type_vector_expr <- reactive({
    req(dataset_config()$spatial_data)
    all_names <- names(get(dataset_config()$spatial_data))
    exclude <- c("samples", "sample_information", "bcs_information", "images_information",
                 "stability_results", "filtered_data", "colfilt_data", "clusters")
    setdiff(all_names, exclude)
  })
  
  # Update inputs when dataset changes
  observeEvent(input$dataset, {
    updateSelectInput(session, "gene_expr", choices = NULL)
    updateSelectInput(session, "peak_expr", choices = NULL)
    
    updateSelectInput(session, "data_type_visualization_expr",
                      choices = data_type_vector_expr(),
                      selected = data_type_vector_expr()[1])
    
    updateTextInput(session, "spatial_data", value = dataset_config()$spatial_data)
  })
  
  # UI: select gene based on data type
  observeEvent(input$data_type_visualization_expr, {
    req(input$data_type_visualization_expr)
    output$gene_expr <- renderUI({
      selectInput("gene_expr", "Select Gene:",
                  choices = spatial_annotate_expr()$gene_name)
    })
  })
  
  # UI: select peak based on selected gene
  observeEvent(input$gene_expr, {
    req(input$gene_expr)
    df <- spatial_annotate_expr()
    filtered_df <- df[df$gene_name == input$gene_expr, ]
    output$peak_expr <- renderUI({
      selectInput("peak_expr", "Select Peak:",
                  choices = unique(filtered_df$peak_id))
    })
  })
  
  # Returns annotation dataframe for selected spatial_data and data type
  annotate_df_expr <- reactive({
    req(dataset_config()$spatial_data, input$data_type_visualization_expr)
    get(dataset_config()$spatial_data)[[input$data_type_visualization_expr]]$annotate
  })
  
  # Returns available gene list
  available_genes_expr <- reactive({
    req(annotate_df_expr())
    unique(annotate_df_expr()$gene_name) %>% sort()
  })
  
  # UI: gene selection input
  output$gene_selector_expr <- renderUI({
    selectInput("gene_selection_expr", "Gene:", choices = available_genes_expr(), selected = "Sgk1")
  })
  
  # Returns peaks available for the selected gene
  available_peaks_expr <- reactive({
    req(annotate_df_expr(), input$gene_selection_expr)
    annotate_df_expr() %>%
      dplyr::filter(gene_name == input$gene_selection_expr) %>%
      dplyr::pull(peak_id) %>%
      unique() %>%
      sort()
  })
  
  # UI: peak selection input
  output$peak_selector_expr <- renderUI({
    selectInput("peak_selection_expr", "Peak:", choices = available_peaks_expr())
  })
  
  # UI: data type selection
  output$data_type_visualization_expr <- renderUI({
    selectInput("data_type_visualization_expr",
                "Select Data Type:",
                choices = data_type_vector_expr(),
                selected = "raw_data")
  })
  
  # Render spatial gene expression plot
  output$expr_spatial_plot_cluster <- renderPlot({
    spatial_feature_plot_cluster(
      spatial_data = get(dataset_config()$spatial_data),
      type_data = input$data_type_visualization_expr,
      peak_id = input$peak_selection_expr,
      clusters = vals$select_clusters_expr,
      samples = c(input$control_selected_expr, input$case_selected_expr),
      min_percentile = input$min_percentile_expr,
      max_percentile = input$max_percentile_expr,
      normalization = input$zero_normalize_expr,
      resolution = vals$resolution_expr,
      size = input$spot_size_expr,
      show_legend = input$show_legend_expr,
      tif_image = input$tif_image_expr
    ) + patchwork::plot_layout(ncol = as.numeric(input$num_columns_expr))
  },
  width = function() as.integer(input$width_plot_expr),
  height = function() as.integer(input$height_plot_expr))
  
  # Download PNG
  output$downloadSpatialClusterPNG_expr <- downloadHandler(
    filename = function() {
      paste0("gene_expr_", input$gene_selection_expr, "_", input$peak_selection_expr, ".png")
    },
    content = function(file) {
      png(file, width = input$width_plot_expr, height = input$height_plot_expr, res = input$dpi_svg_expr)
      print(
        spatial_feature_plot_cluster(
          spatial_data = get(dataset_config()$spatial_data),
          type_data = input$data_type_visualization_expr,
          peak_id = input$peak_selection_expr,
          clusters = vals$select_clusters_expr,
          samples = c(input$control_selected_expr, input$case_selected_expr),
          min_percentile = input$min_percentile_expr,
          max_percentile = input$max_percentile_expr,
          normalization = input$zero_normalize_expr,
          resolution = input$cluster_resolution_slider_expr,
          size = input$spot_size_expr,
          show_legend = input$show_legend_expr,
          tif_image = input$tif_image_expr
        ) +
          patchwork::plot_layout(ncol = as.numeric(input$num_columns_expr))
      )
      dev.off()
    }
  )
  
  # Download SVG
  output$downloadSpatialClusterSVG_expr <- downloadHandler(
    filename = function() {
      paste0("gene_expr_", input$gene_selection_expr, "_", input$peak_selection_expr, ".svg")
    },
    content = function(file) {
      svg(file,
          width = input$width_plot_expr / 100,
          height = input$height_plot_expr / 100)
      print(
        spatial_feature_plot_cluster(
          spatial_data = get(dataset_config()$spatial_data),
          type_data = input$data_type_visualization_expr,
          peak_id = input$peak_selection_expr,
          clusters = vals$select_clusters_expr,
          samples = c(input$control_selected_expr, input$case_selected_expr),
          min_percentile = input$min_percentile_expr,
          max_percentile = input$max_percentile_expr,
          normalization = input$zero_normalize_expr,
          resolution = input$cluster_resolution_slider_expr,
          size = input$spot_size_expr,
          show_legend = input$show_legend_expr,
          tif_image = input$tif_image_expr
        ) +
          patchwork::plot_layout(ncol = as.numeric(input$num_columns_expr))
      )
      dev.off()
    }
  )
  
  
  # -----------------------------
  # Section: Gene Expression (Spatial, Aggregate Cluster)
  # -----------------------------
  
  # UI: sample selection block (control vs experimental)
  output$sample_display_block_agg <- renderUI({
    req(dataset_config())
    control <- get(dataset_config()$control_samples)
    case <- get(dataset_config()$case_samples)
    
    fluidRow(
      column(6,
             h5(strong("Control Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_control_agg", "Select All"),
                 actionButton("deselect_all_control_agg", "Deselect All")
             ),
             checkboxGroupInput("control_selected_agg", NULL, choices = control, selected = control)
      ),
      column(6,
             h5(strong("Experimental Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_case_agg", "Select All"),
                 actionButton("deselect_all_case_agg", "Deselect All")
             ),
             checkboxGroupInput("case_selected_agg", NULL, choices = case, selected = case)
      )
    )
  })
  
  # Obsługa przycisków wyboru próbek
  observeEvent(input$select_all_control_agg, {
    updateCheckboxGroupInput(session, "control_selected_agg",
                             selected = get(dataset_config()$control_samples))
  })
  observeEvent(input$deselect_all_control_agg, {
    updateCheckboxGroupInput(session, "control_selected_agg", selected = character(0))
  })
  observeEvent(input$select_all_case_agg, {
    updateCheckboxGroupInput(session, "case_selected_agg",
                             selected = get(dataset_config()$case_samples))
  })
  observeEvent(input$deselect_all_case_agg, {
    updateCheckboxGroupInput(session, "case_selected_agg", selected = character(0))
  })
  
  # UI: gene select
  output$gene_selector_agg <- renderUI({
    req(annotate_df_expr())
    selectInput("gene_selection_agg", "Gene:",
                choices = available_genes_expr(), selected = "Sgk1")
  })
  
  # Returns peaks available for the selected gene (AGG)
  available_peaks_agg <- reactive({
    req(annotate_df_expr(), input$gene_selection_agg)
    annotate_df_expr() %>%
      dplyr::filter(gene_name == input$gene_selection_agg) %>%
      dplyr::pull(peak_id) %>%
      unique() %>%
      sort()
  })
  
  # UI: peak selection input (AGG)
  output$peak_selector_agg <- renderUI({
    selectInput("peak_selection_agg", "Peak:", choices = available_peaks_agg())
  })
  
  # UI: data type
  output$data_type_visualization_agg <- renderUI({
    selectInput("data_type_visualization_agg", "Select Data Type:",
                choices = data_type_vector_expr(), selected = "raw_data")
  })
  
  # Plot
  output$spatial_agg_plot <- renderPlot({
    req(input$data_type_visualization_agg, input$peak_selection_agg)
    
    spatial_feature_aggregated_cluster_plot(
      spatial_data = get(dataset_config()$spatial_data),
      summary_statistics = get(dataset_config()$summary_data),
      peak_id = input$peak_selection_agg,
      type_data = input$data_type_visualization_agg,
      samples = c(input$control_selected_agg, input$case_selected_agg),
      cluster_resolution = input$cluster_resolution_agg,
      summary_metric = input$summary_metric_agg,
      normalization = input$zero_normalize_agg,
      tif_image = input$tif_image_agg,
      show_legend = input$show_legend_agg,
      return_list = FALSE,
      size = input$spot_size_agg
    ) + patchwork::plot_layout(ncol = as.numeric(input$num_columns_agg))
  },
  width = function() as.integer(input$width_plot_agg),
  height = function() as.integer(input$height_plot_agg))
  
  # Download PNG
  output$downloadSpatialClusterPNG_agg <- downloadHandler(
    filename = function() {
      paste0("agg_expr_", input$gene_selection_agg, "_", input$peak_selection_agg, ".png")
    },
    content = function(file) {
      png(file, width = input$width_plot_agg, height = input$height_plot_agg, res = input$dpi_svg_agg)
      print(
        spatial_feature_aggregated_cluster_plot(
          spatial_data = get(dataset_config()$spatial_data),
          summary_statistics = get(dataset_config()$summary_data),
          peak_id = input$peak_selection_agg,
          type_data = input$data_type_visualization_agg,
          samples = c(input$control_selected_agg, input$case_selected_agg),
          cluster_resolution = input$cluster_resolution_agg,
          summary_metric = input$summary_metric_agg,
          normalization = input$zero_normalize_agg,
          tif_image = input$tif_image_agg,
          show_legend = input$show_legend_agg,
          return_list = FALSE,
          size = input$spot_size_agg
        ) + patchwork::plot_layout(ncol = as.numeric(input$num_columns_agg))
      )
      dev.off()
    }
  )
  
  # Download SVG
  output$downloadSpatialClusterSVG_agg <- downloadHandler(
    filename = function() {
      paste0("agg_expr_", input$gene_selection_agg, "_", input$peak_selection_agg, ".svg")
    },
    content = function(file) {
      svg(file,
          width = input$width_plot_agg / 100,
          height = input$height_plot_agg / 100)
      print(
        spatial_feature_aggregated_cluster_plot(
          spatial_data = get(dataset_config()$spatial_data),
          summary_statistics = get(dataset_config()$summary_data),
          peak_id = input$peak_selection_agg,
          type_data = input$data_type_visualization_agg,
          samples = c(input$control_selected_agg, input$case_selected_agg),
          cluster_resolution = input$cluster_resolution_agg,
          summary_metric = input$summary_metric_agg,
          normalization = input$zero_normalize_agg,
          tif_image = input$tif_image_agg,
          show_legend = input$show_legend_agg,
          return_list = FALSE,
          size = input$spot_size_agg
        ) + patchwork::plot_layout(ncol = as.numeric(input$num_columns_agg))
      )
      dev.off()
    }
  )
  
}


# helper for null-safe fallback
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==== APP ====
shinyApp(ui, server)
