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
                 h4("Visualization Options"),
                 
                 # Spot size + number of columns
                 fluidRow(
                   column(6,
                          numericInput("spot_size", "Spot Size", value = 1, min = 0.1, max = 5, step = 0.1)
                   ),
                   column(6,
                          selectInput("num_columns", "Number of Columns", choices = 1:7, selected = 4)
                   )
                 ),
                 
                 # Width + height in pixels
                 fluidRow(
                   column(6,
                          numericInput("width_plot", "Width (px):", value = 1100, min = 100, step = 50)
                   ),
                   column(6,
                          numericInput("height_plot", "Height (px):", value = 950, min = 100, step = 50)
                   )
                 ),
                 
                 # DPI for SVG export
                 fluidRow(
                   column(6,
                          numericInput("dpi_svg", "SVG DPI:", value = 100, min = 50, max = 600, step = 10)
                   )
                 ),
                 
                 checkboxInput("tif_image", "Show Background Image", value = TRUE),
                 
                 tags$hr(),
                 h4("Clustering Resolution"),
                 
                 sliderInput("cluster_resolution_slider",
                             label = "Select resolution:",
                             min = 0.05, max = 2, value = 0.4, step = 0.05),
                 
                 tags$hr(),
                 h4("Cluster Selection"),
                 
                 div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                     actionButton("select_all_clusters", "Select All"),
                     actionButton("deselect_all_clusters", "Deselect All")
                 ),
                 
                 # 3 columns for clusters
                 fluidRow(
                   lapply(0:2, function(col_index) {
                     column(
                       4,
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
               )
               
               ,
               mainPanel(
                 h4("Cluster Visualization"),
                 plotOutput("clust_plot", height = "600px")
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
    
    # === UI: Gene Expression (Spatial) Tab ===
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
                 uiOutput("data_type_visualization_expr"),   # dynamiczne selectInput z serwera
                 uiOutput("gene_selector_expr"),             # dynamiczne selectInput z serwera
                 uiOutput("peak_selector_expr"),             # dynamiczne selectInput z serwera
                 
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
                     column(
                       4,
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
               )
               ,
               
               mainPanel(
                 h4("Spatial Gene Expression Plot"),
                 plotOutput("expr_spatial_plot_cluster", height = "600px")
               )
             )
    )
    
    ,
    
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
  
  # ==== Reactive values ====
  vals <- reactiveValues(
    resolution = 0.5,
    select_clusters = as.character(0:29)  # domyślnie wszystkie klastry zaznaczone
  )
  
  # ==== Dataset config ====
  dataset_config <- reactive({
    req(input$dataset)
    dataset_registry[[input$dataset]]
  })
  
  # ==== Cluster selection logic ====
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
  
  # ==== Sample selection logic ====
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
  
  # ==== Sample display block ====
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
             checkboxGroupInput("control_selected", NULL,
                                choices = control,
                                selected = control)
      ),
      column(6,
             h5(strong("Experimental Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_case", "Select All"),
                 actionButton("deselect_all_case", "Deselect All")
             ),
             checkboxGroupInput("case_selected", NULL,
                                choices = case,
                                selected = case)
      )
    )
  })
  
  # ==== Logo ====
  output$logo_left <- renderPlot({
    img <- png::readPNG("spatial-browser-logo.png")
    grid::grid.raster(img)
  })
  
  # ==== Cluster plot ====
  output$clust_plot <- renderPlot({
    
    if (is.null(input$dataset) || input$dataset == "") {
      plot.new(); text(0.5, 0.5, "No dataset selected.", cex = 1.5, col = "red"); return()
    }
    
    selected_samples <- c(input$control_selected, input$case_selected)
    if (length(selected_samples) == 0) {
      plot.new(); text(0.5, 0.5, "No samples selected.", cex = 1.5, col = "red"); return()
    }
    
    if (is.null(vals$select_clusters) || length(vals$select_clusters) == 0) {
      plot.new(); text(0.5, 0.5, "No clusters selected.", cex = 1.5, col = "red"); return()
    }
    
    if (is.null(vals$resolution)) {
      plot.new(); text(0.5, 0.5, "Resolution not set.", cex = 1.5, col = "red"); return()
    }
    
    spatial_cluster_select(
      spatial_data = get(dataset_config()$spatial_data),
      resolution = vals$resolution,
      samples = selected_samples,
      palette = palette_allen,
      size = input$spot_size,
      select_clusters = vals$select_clusters,
      tif_image = input$tif_image,
      ncol = as.numeric(input$num_columns)
    ) +
      patchwork::plot_layout(ncol = as.numeric(input$num_columns))
    
  },
  # ✅ Dynamiczne wymiary
  width = function() as.integer(input$width_plot),
  height = function() as.integer(input$height_plot)
  )
  
  # ==== UMAP, summary, expr, etc. (placeholder examples) ====
  output$umap_plot <- renderPlot({
    plot(runif(100), runif(100), col = "gray", pch = 16, main = "UMAP Plot")
  })
  
  output$stats_table <- renderTable({
    data.frame(Gene = c("GeneA", "GeneB"), log2Ratio = c(1.2, -0.8), p_value = c(0.01, 0.03))
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
  
  output$downloadSpatialClusterPNG <- downloadHandler(
    filename = function() {
      paste0("cluster-", input$dataset, "-res", vals$resolution, ".png")
    },
    content = function(file) {
      png(file,
          width = as.integer(input$width_plot),
          height = as.integer(input$height_plot))
      print(
        spatial_cluster_select(
          spatial_data = get(dataset_config()$spatial_data),
          resolution = vals$resolution,
          samples = c(input$control_selected, input$case_selected),
          palette = palette_allen,
          size = input$spot_size,
          select_clusters = vals$select_clusters,
          tif_image = input$tif_image,
          ncol = as.numeric(input$num_columns)
        ) +
          patchwork::plot_layout(ncol = as.numeric(input$num_columns))
      )
      dev.off()
    }
  )
  
  output$downloadSpatialClusterSVG <- downloadHandler(
    filename = function() {
      paste0("cluster-", input$dataset, "-res", vals$resolution, ".svg")
    },
    content = function(file) {
      dpi <- 100
      svg(file,
          width = as.integer(input$width_plot) / dpi,
          height = as.integer(input$height_plot) / dpi)
      print(
        spatial_cluster_select(
          spatial_data = get(dataset_config()$spatial_data),
          resolution = vals$resolution,
          samples = c(input$control_selected, input$case_selected),
          palette = palette_allen,
          size = input$spot_size,
          select_clusters = vals$select_clusters,
          tif_image = input$tif_image,
          ncol = as.numeric(input$num_columns)
        ) +
          patchwork::plot_layout(ncol = as.numeric(input$num_columns))
      )
      dev.off()
    }
  )
  
  # ==== Dynamic UI for Gene Expression (Spatial) ====
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
             checkboxGroupInput("control_selected_expr", NULL,
                                choices = control,
                                selected = control)
      ),
      column(6,
             h5(strong("Experimental Group")),
             div(style = "display: flex; gap: 6px; margin-bottom: 8px;",
                 actionButton("select_all_case_expr", "Select All"),
                 actionButton("deselect_all_case_expr", "Deselect All")
             ),
             checkboxGroupInput("case_selected_expr", NULL,
                                choices = case,
                                selected = case)
      )
    )
  })
  
  # Zwraca dostępne adnotacje (zależne od spatial_data i typu danych)
  spatial_annotate_expr <- reactive({
    req(dataset_config()$spatial_data, input$data_type_visualization_expr)
    get(dataset_config()$spatial_data)[[input$data_type_visualization_expr]]$annotate
  })
  
  # Zwraca dostępne typy danych (raw, normalized, seurat itd.)
  data_type_vector_expr <- reactive({
    req(dataset_config()$spatial_data)
    get(dataset_config()$spatial_data) %>%
      names() %>%
      .[grepl(paste(c("raw_data", "range_normalize", "quantile_normalize", "seurat"), collapse = "|"), .)]
  })
  
  # Gdy zmieni się dataset, wyczyść pola i ustaw domyślny typ danych
  observeEvent(input$dataset, {
    # Przeładuj dynamiczne pola po zmianie datasetu
    updateSelectInput(session, "gene_expr", choices = NULL)
    updateSelectInput(session, "peak_expr", choices = NULL)
    
    # Zaktualizuj dostępne typy danych dla aktualnego spatial_data
    updateSelectInput(session, "data_type_visualization_expr",
                      choices = data_type_vector_expr(),
                      selected = data_type_vector_expr()[1])
    
    # Ustaw spatial_data jako ukryte input — potrzebne dla reactives bazujących na input$spatial_data
    updateTextInput(session, "spatial_data", value = dataset_config()$spatial_data)
  })
  # UI: wybór genu po zmianie typu danych
  observeEvent(input$data_type_visualization_expr, {
    req(input$data_type_visualization_expr)
    output$gene_expr <- renderUI({
      selectInput("gene_expr", "Select Gene:",
                  choices = spatial_annotate_expr()$gene_name)
    })
  })
  
  # UI: wybór peaku po wyborze genu
  observeEvent(input$gene_expr, {
    req(input$gene_expr)
    df <- spatial_annotate_expr()
    filtered_df <- df[df$gene_name == input$gene_expr, ]
    output$peak_expr <- renderUI({
      selectInput("peak_expr", "Select Peak:",
                  choices = unique(filtered_df$peak_id))
    })
  })
  
  data_type_vector_expr <- reactive({
    req(input$spatial_data)
    get(input$spatial_data) %>%
      names() %>%
      .[grepl(paste(c("raw_data", "range_normalize", "quantile_normalize", "seurat"), collapse = "|"), .)]
  })
  
  observeEvent(input$spatial_data, {
    output$data_type_visualization_expr <- renderUI({
      selectInput("data_type_visualization_expr",
                  "Select Data Type:",
                  choices = data_type_vector_expr())
    })
  })
  
  # observeEvent(input$data_type_visualization_expr, {
  #   req(input$data_type_visualization_expr)
  #   output$gene_expr <- renderUI({
  #     selectInput("gene_expr",
  #                 "Select Gene:",
  #                 choices = spatial_annotate_expr()$gene_name)
  #   })
  # })
  # 
  # observeEvent(input$gene_expr, {
  #   req(input$gene_expr, input$data_type_visualization_expr)
  #   df <- spatial_annotate_expr()
  #   filtered_df <- df[df$gene_name == input$gene_expr, ]
  #   output$peak_expr <- renderUI({
  #     selectInput("peak_expr",
  #                 "Select Peak:",
  #                 choices = unique(filtered_df$peak_id))
  #   })
  # })
  
  # Pobranie adnotacji dla wybranego spatial_data i type_data
  annotate_df_expr <- reactive({
    req(dataset_config()$spatial_data, input$data_type_visualization_expr)
    get(dataset_config()$spatial_data)[[input$data_type_visualization_expr]]$annotate
  })
  
  # Dynamiczna lista genów
  available_genes_expr <- reactive({
    req(annotate_df_expr())
    unique(annotate_df_expr()$gene_name) %>% sort()
  })
  
  # Dynamiczne UI do wyboru genu
  output$gene_selector_expr <- renderUI({
    selectInput("gene_selection_expr",
                "Gene:",
                choices = available_genes_expr(),
                selected = "Sgk1")
  })
  
  # Dynamiczna lista peaków zależna od wybranego genu
  available_peaks_expr <- reactive({
    req(annotate_df_expr(), input$gene_selection_expr)
    annotate_df_expr() %>%
      dplyr::filter(gene_name == input$gene_selection_expr) %>%
      dplyr::pull(peak_id) %>%
      unique() %>%
      sort()
  })
  
  # Dynamiczne UI do wyboru peaku
  output$peak_selector_expr <- renderUI({
    selectInput("peak_selection_expr",
                "Peak:",
                choices = available_peaks_expr())
  })
  
  data_type_vector_expr <- reactive({
    req(dataset_config()$spatial_data)
    all_names <- names(get(dataset_config()$spatial_data))
    exclude <- c(
      "samples", "sample_information", "bcs_information", "images_information",
      "stability_results", "filtered_data", "colfilt_data", "clusters"
    )
    setdiff(all_names, exclude)
  })
  
  output$data_type_visualization_expr <- renderUI({
    selectInput("data_type_visualization_expr",
                "Select Data Type:",
                choices = data_type_vector_expr(),
                selected = "raw_data")
  })
  
  output$expr_spatial_plot_cluster <- renderPlot({
    spatial_feature_plot_cluster(
      spatial_data = get(dataset_config()$spatial_data),
      type_data = input$data_type_visualization_expr,  # <-- DYNAMICZNY TYP DANYCH
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
}

# helper for null-safe fallback
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==== APP ====
shinyApp(ui, server)
