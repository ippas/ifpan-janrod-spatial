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
  
  # -----------------------------
  # Section: Statistics Summary Table
  # -----------------------------
  
  filtered_stats <- reactive({
    req(dataset_config())
    
    filter_data_statistics(
      summary_data = get(dataset_config()$summary_data),
      data_type = input$stat_data_type,
      resolution = input$stat_resolution,
      metric = input$stat_summary_metric,
      control_mean_threshold = input$control_mean_thresh,
      experiment_mean_threshold = input$experiment_mean_thresh,
      log2ratio_threshold = input$log2ratio_thresh,
      t_test_threshold = if (input$stat_test_type == "T-Student") input$pvalue_thresh else 1,
      wilcoxon_test_threshold = if (input$stat_test_type == "Wilcoxon") input$pvalue_thresh else 1,
      ks_test_threshold = if (input$stat_test_type == "Kolmogorov-Smirnov") input$pvalue_thresh else 1
    ) %>%
      dplyr::mutate_if(is.numeric, ~ round(.x, 4)) %>%
      dplyr::select(-condition)  # jeśli istnieje taka kolumna
  })
  
  output$stats_table <- DT::renderDataTable({
    DT::datatable(
      filtered_stats(),
      options = list(
        paging = TRUE,
        scrollX = TRUE,
        scrollY = TRUE,
        pageLength = 20,
        columnDefs = list(list(targets = "_all", className = "dt-center"))
      ),
      extensions = "Buttons",
      filter = "top",
      selection = "single",
      rownames = FALSE
    ) %>%
      DT::formatStyle("peak", `white-space` = "nowrap")
  })
  
  output$download_stats_tsv <- downloadHandler(
    filename = function() {
      paste0("results_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      write.table(filtered_stats(), file, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  
  output$download_stats_csv <- downloadHandler(
    filename = function() {
      paste0("results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_stats(), file, row.names = FALSE)
    }
  )
  
  #  -----------------------------
  #   SERVER LOGIC FOR: Gene Expression (Barplot)
  #  -----------------------------
  
  # Cluster selection logic for barplot tab
  observeEvent(input$update_resolution_bar, {
    selected_clusters <- unlist(c(
      input$cluster_selector_bar_col1 %||% character(0),
      input$cluster_selector_bar_col2 %||% character(0),
      input$cluster_selector_bar_col3 %||% character(0)
    ))
    vals$select_clusters_bar <- gsub("^cluster_", "", selected_clusters)
  })
  
  observeEvent(input$select_all_clusters_bar, {
    lapply(0:2, function(col_index) {
      updateCheckboxGroupInput(
        session,
        inputId = paste0("cluster_selector_bar_col", col_index + 1),
        selected = paste0("cluster_", (col_index * 10):(col_index * 10 + 9))
      )
    })
  })
  
  observeEvent(input$deselect_all_clusters_bar, {
    lapply(1:3, function(i) {
      updateCheckboxGroupInput(session, paste0("cluster_selector_bar_col", i), selected = character(0))
    })
  })
  
  # UI for gene & peak selection (Barplot)
  output$data_type_visualization_bar <- renderUI({
    selectInput("data_type_visualization_bar",
                "Select Data Type:",
                choices = data_type_vector_expr(),
                selected = "raw_data")
  })
  
  output$gene_selector_bar <- renderUI({
    selectInput("gene_selection_bar", "Gene:", choices = available_genes_expr(), selected = "Sgk1")
  })
  
  available_peaks_bar <- reactive({
    req(annotate_df_expr(), input$gene_selection_bar)
    annotate_df_expr() %>%
      dplyr::filter(gene_name == input$gene_selection_bar) %>%
      dplyr::pull(peak_id) %>%
      unique() %>%
      sort()
  })
  
  output$peak_selector_bar <- renderUI({
    selectInput("peak_selection_bar", "Peak:", choices = available_peaks_bar())
  })
  
  # Render barplot expression
  output$barplot_expr <- renderPlot({
    req(input$stat_data_type_bar, input$stat_resolution_bar, input$stat_summary_metric_bar,
        input$gene_selection_bar, input$peak_selection_bar, input$stat_test_type_bar)
    
    resolution_label <- paste0("resolution_", input$stat_resolution_bar)
    
    barplot_feature_expression_clusters(
      summary_data = get(dataset_config()$summary_data),
      data_type = input$stat_data_type_bar,
      resolution = resolution_label,
      metric = input$stat_summary_metric_bar,
      selected_peak = input$peak_selection_bar,
      test = switch(input$stat_test_type_bar,
                    "T-Student" = "t_test",
                    "Wilcoxon" = "wilcoxon_test",
                    "Kolmogorov-Smirnov" = "ks_test"),
      same_y_scale = input$same_y_scale_bar,
      clusters = vals$select_clusters_bar,
      ncol = as.numeric(input$num_columns_bar)
    )
  },
  width = function() as.integer(input$width_plot_bar),
  height = function() as.integer(input$height_plot_bar))
  
  # Downloads for barplot
  output$downloadBarplotPNG <- downloadHandler(
    filename = function() {
      paste0("barplot_expr_", input$gene_selection_bar, "_", input$peak_selection_bar, ".png")
    },
    content = function(file) {
      png(file, width = input$width_plot_bar, height = input$height_plot_bar, res = input$dpi_svg_bar)
      print(
        barplot_feature_expression_clusters(
          summary_data = get(dataset_config()$summary_data),
          data_type = input$stat_data_type_bar,
          resolution = paste0("resolution_", input$stat_resolution_bar),
          metric = input$stat_summary_metric_bar,
          selected_peak = input$peak_selection_bar,
          test = switch(input$stat_test_type_bar,
                        "T-Student" = "t_test",
                        "Wilcoxon" = "wilcoxon_test",
                        "Kolmogorov-Smirnov" = "ks_test"),
          same_y_scale = input$same_y_scale_bar,
          clusters = vals$select_clusters_bar,
          ncol = as.numeric(input$num_columns_bar)
        )
      )
      dev.off()
    }
  )
  
  output$downloadBarplotSVG <- downloadHandler(
    filename = function() {
      paste0("barplot_expr_", input$gene_selection_bar, "_", input$peak_selection_bar, ".svg")
    },
    content = function(file) {
      svg(file,
          width = input$width_plot_bar / 100,
          height = input$height_plot_bar / 100)
      print(
        barplot_feature_expression_clusters(
          summary_data = get(dataset_config()$summary_data),
          data_type = input$stat_data_type_bar,
          resolution = paste0("resolution_", input$stat_resolution_bar),
          metric = input$stat_summary_metric_bar,
          selected_peak = input$peak_selection_bar,
          test = switch(input$stat_test_type_bar,
                        "T-Student" = "t_test",
                        "Wilcoxon" = "wilcoxon_test",
                        "Kolmogorov-Smirnov" = "ks_test"),
          same_y_scale = input$same_y_scale_bar,
          clusters = vals$select_clusters_bar,
          ncol = as.numeric(input$num_columns_bar)
        )
      )
      dev.off()
    }
  )
}
