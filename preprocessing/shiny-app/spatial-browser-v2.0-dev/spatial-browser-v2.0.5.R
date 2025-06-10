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
                 
                 # Width + height
                 fluidRow(
                   column(6,
                          numericInput("width_plot", "Width:", value = 1100, min = 100, step = 50)
                   ),
                   column(6,
                          numericInput("height_plot", "Height:", value = 950, min = 100, step = 50)
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
}

# helper for null-safe fallback
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==== APP ====
shinyApp(ui, server)