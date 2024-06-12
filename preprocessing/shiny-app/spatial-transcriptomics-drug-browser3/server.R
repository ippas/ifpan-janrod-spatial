# Define the server logic
server <- function(input, output, session) {
  vals <- reactiveValues(cluster =  0, resolution = 0.05, select_clusters = c(0:32))
  
  observeEvent(input$update_values, {
    vals$cluster <- input$cluster
    vals$resolution <- input$resolution
    vals$select_clusters <- reactive_selected_clusters()
  })
  
  # Reactive expression to create a vector of selected clusters
  reactive_selected_clusters <- reactive({
    # Initialize an empty vector
    selected <- c()
    
    # Loop through each checkbox and check if it is selected
    for (i in 0:32) {
      if (input[[paste("cluster", i, sep = "")]]) {
        selected <- c(selected, i)
      }
    }
    
    # Return the vector of selected clusters
    selected
  })
  
  # # Reactive expression for selected saline samples
  # selected_saline_samples <- reactive({
  #   saline_ids <- c("S6269Nr1", "S6269Nr3", "S7788Nr1", "S7788Nr2", "S7788Nr3", "S7788Nr11")
  #   selected <- lapply(saline_ids, function(x) {
  #     if(input[[paste("saline_", gsub("Nr|S", "", x), sep = "")]]) x
  #   })
  #   unlist(selected)
  # })
  # 
  # # Reactive expression for selected experimental samples
  # selected_experimental_samples <- reactive({
  #   experimental_ids <- c("S6230Nr3", "S6230Nr4", "S6269Nr2", "S6269Nr4", "S7788Nr15", "S7788Nr16")
  #   selected <- lapply(experimental_ids, function(x) {
  #     if(input[[paste("experimental_", gsub("Nr|S", "", x), sep = "")]]) x
  #   })
  #   unlist(selected)
  # })
  
  selected_saline_samples <- eventReactive(input$update_samples, {
    saline_ids <- c( "S6269Nr1",  "S6269Nr3",  "S7788Nr1",  "S7788Nr2",  "S7788Nr3",  "S7788Nr11",
                     "S5023Nr1", "S5023Nr4", "S5295Nr1", "S5295Nr6", "S5295Nr7")
    selected <- lapply(saline_ids, function(x) {
      input_id <- paste("saline_", gsub("Nr|S", "", x), sep = "")
      if(isTruthy(input[[input_id]])) x
    })
    selected <- Filter(Negate(is.null), selected)  # Remove NULLs
    unlist(selected)
  }, ignoreNULL = FALSE)
  
  selected_experimental_samples <- eventReactive(input$update_samples, {
    experimental_ids <- c("S6230Nr3", "S6230Nr4", "S6269Nr2", "S6269Nr4", "S7788Nr15", "S7788Nr16",
                          "S5023Nr2", "S5023Nr3", "S5295Nr3", "S5295Nr8")
    selected <- lapply(experimental_ids, function(x) {
      input_id <- paste("experimental_", gsub("Nr|S", "", x), sep = "")
      if(isTruthy(input[[input_id]])) x
    })
    selected <- Filter(Negate(is.null), selected)  # Remove NULLs
    unlist(selected)
  }, ignoreNULL = FALSE)

  observeEvent(input$drug, {
    
    if (input$drug == "Risperidone") {
      
      updateSelectInput(session, "spatial_data",
                        choices = c("risperidone_st_data_half"),
                        selected = "risperidone_st_data_half")
      
      updateSelectInput(session, "summary_data",
                        choices = c("risperidone_summary_statistics_half"),
                        selected = "risperidone_summary_statistics_half")
      
      updateSelectInput(session, "experiment_samples",
                        choices = c("samples_risperidone"),
                        selected = "samples_risperidone")
      
    } else if (input$drug == "PZ-1190") {
      
      updateSelectInput(session, "spatial_data",
                        choices = c("pz1190_st_data_half"),
                        selected = "pz1190_st_data_half")
      
      updateSelectInput(session, "summary_data",
                        choices = c("pz1190_summary_statistics_half"),
                        selected = "pz1190_summary_statistics_half")
      
      updateSelectInput(session, "experiment_samples",
                        choices = c("samples_pz1190"),
                        selected = "samples_pz1190")
      
    }  else if (input$drug == "Ldopa") {
      
      updateSelectInput(session, "spatial_data",
                        choices = c("ldopa_st_data"),
                        selected = "ldopa_st_data")
      
      updateSelectInput(session, "summary_data",
                        choices = c("ldopa_summary_statistics"),
                        selected = "ldopa_summary_statistics")
      
      updateSelectInput(session, "experiment_samples",
                        choices = c("samples_ldopa"),
                        selected = "samples_ldopa")
      
    } else if (input$drug == "Clozapine") {
      
      updateSelectInput(session, "spatial_data",
                        choices = c("clozapine_st_data_half"),
                        selected = "clozapine_st_data_half")
      
      updateSelectInput(session, "summary_data",
                        choices = c("clozapine_summary_statistics_half"),
                        selected = "clozapine_summary_statistics_half")
      
      updateSelectInput(session, "experiment_samples",
                        choices = c("samples_clozapine"),
                        selected = "samples_clozapine")
      
    }
    
    # You can add more conditions for additional drugs here
    
  })
  
  drug_samples <- list(
    "Risperidone" = c("S6230Nr3", "S6230Nr4", "S6269Nr2", "S6269Nr4", "S7788Nr15", "S7788Nr16"),
    "PZ-1190" = c("S7788Nr4", "S7788Nr6", "S7788Nr7", "S7788Nr9", "S7788Nr10", "S7788Nr14"),
    "Clozapine" = c("S6230Nr1", "S6230Nr2", "S7788Nr5", "S7788Nr8", "S7788Nr12", "S7788Nr13"),
    "Ldopa" = c("S5023Nr2", "S5023Nr3", "S5295Nr3", "S5295Nr8")
  )
  
  
  spatial_annotate <- reactive({get(input$spatial_data)[[input$data_type_visualization]]$annotate})
  
  data_type_vector <-  reactive({get(input$spatial_data) %>%
      names() %>% 
      .[grepl(paste(c("raw_data", "range_normalize", "quantile_normalize", "seurat"), collapse = "|"), .)]})
  
  observeEvent(input$spatial_data, {
    output$data_type_visualization <- renderUI({
      selectInput("data_type_visualization", 
                  "Select Data Type:", 
                  choices = data_type_vector())
    })
    
    output$gene <- renderUI({
      selectInput("gene", 
                  "Select Gene:", 
                  choices = spatial_annotate() %>% .$gene_name)
    })
  })
  
  observeEvent(input$clear_genes, {
    # Update the text area input with an empty string, effectively clearing it.
    updateTextAreaInput(session, "gene_names", value = "")
  })
  
  # Reactive expression for gene vector with debounce
  gene_vector <- reactive({
    gene_vector <- strsplit(input$gene_names, split = "\n")[[1]]
    gene_vector <- trimws(gene_vector)
    gene_vector <- gsub("^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$", "", gene_vector)
    gene_vector <- gene_vector[nzchar(gene_vector)]
    gene_vector <- tolower(gene_vector)
    print(gene_vector)
    
    return(gene_vector)
  }) %>% debounce(500)  # 500ms debounce
  
  # Reactive expression for filtered data
  filtered_data <- reactive({
    # Initial filtering steps
    data <- filter_data_statistics(
      summary_data = get(input$summary_data), 
      data_type = input$data_type_statistics, 
      resolution = input$resolution_statistics,
      metric = input$metric,
      control_mean_threshold = input$control_mean_threshold,
      experiment_mean_threshold = input$experiment_mean_threshold,
      log2ratio_threshold = input$log2ratio_threshold,
      t_test_threshold = input$t_test_threshold,
      wilcoxon_test_threshold = input$wilcoxon_test_threshold,
      ks_test_threshold = input$ks_test_threshold
    )
    
    # print("##################################################################")
    # print(input$summary_data)
    # print(data)
    
    # Get gene_vector
    gene_vector <- gene_vector()
    
    # Filter by genes if gene_vector is not empty
    if (length(gene_vector) > 0) {
      data <- data %>%
        # filter(gene %in% gene_vector)
        mutate(gene_filt = tolower(gene)) %>%
        filter(gene_filt %in% gene_vector) %>%
        select(-gene_filt)
      # data[data$gene %in% gene_vector, ]
    }
    
    # Return the filtered data
    return(data)
  })
  
  # Use the filtered_data for rendering the DataTable
  output$filter_statistics <- renderDataTable({
    datatable_data <- filtered_data() %>%
      mutate_if(is.numeric, round, 4) %>%
      select(-condition)
    
    datatable(datatable_data,
              options = list(paging = TRUE,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             server = FALSE,
                             pageLength = 20,
                             # buttons = c('csv', 'excel'),
                             columnDefs = list(list(targets = '_all', className = 'dt-center'))
              ),
              extensions = 'Buttons',
              filter = 'top',
              selection = 'single',
              rownames = FALSE
    ) %>% formatStyle("peak", "white-space" = "nowrap")
  })
  
  # Use the filtered_data for downloading
  output$downloadData <- downloadHandler(
    filename = function() {
      input$filename
    },
    content = function(file) {
      datatable_data <- filtered_data() %>%
        mutate_if(is.numeric, round, 4) %>%
        select(-condition)
      write.table(datatable_data, file, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  
  ### code to perform visualization
  
  
  # Observe when the gene selection changes and update the peak dropdown accordingly
  observeEvent(input$gene, {
    df <- spatial_annotate()
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
    # spatial_feature_plot(spatial_data = get(input$spatial_data),
    #                      type_data = {{input$data_type_visualization}},
    #                      peak_id = {{input$peak}},
    #                      # samples =  c(samples_saline[-1], samples_risperidone[-3]),
    #                      # samples = c(samples_saline, get(input$experiment_samples)),
    #                      # samples = c(selected_saline_samples(), selected_experimental_samples()),
    #                      c(selected_saline_samples(), selected_experimental_samples()),
    #                      min_percentile = input$min_percentile,
    #                      max_percentile = input$max_percentile,
    #                      normalization = {{input$zero_normalize}},
    #                      size = input$spot_size,
    #                      show_legend = input$show_legend,
    #                      tif_image = input$tif_image) +
    #   plot_layout(ncol = as.numeric({{input$num_columns}}))
    
    spatial_feature_plot_cluster(spatial_data = get(input$spatial_data),
                         type_data = {{input$data_type_visualization}},
                         peak_id = {{input$peak}},
                         clusters = vals$select_clusters,
                         # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                         # samples = c(samples_saline, get(input$experiment_samples)),
                         # samples = c(selected_saline_samples(), selected_experimental_samples()),
                         c(selected_saline_samples(), selected_experimental_samples()),
                         min_percentile = input$min_percentile,
                         max_percentile = input$max_percentile,
                         normalization = {{input$zero_normalize}},
                         resolution = vals$resolution,
                         size = input$spot_size,
                         show_legend = input$show_legend,
                         tif_image = input$tif_image
                         ) +
      plot_layout(ncol = as.numeric({{input$num_columns}}))
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) })
  
  # Define the server-side function for the download
  output$downloadSpatialFeaturePlotPNG <- downloadHandler(
    # Set the filename of the download
    # This function is called when the user clicks the download button
    # We name it according to the gene_name and peak from the user's input
    filename = function() {
      # The paste function concatenates the gene name, peak, and "_plot.png"
      # The sep = '' argument means that there will be no separation between the concatenated elements
      paste0(input$gene, '-', input$peak, '-', input$data_type_visualization, '.png')
    },
    
    # Define what will be downloaded
    # This function is also called when the user clicks the download button
    content = function(file) {
      # Open a png device
      png(filename = file, width = as.integer(input$width_plot), height = as.integer(input$height_plot))
      
      # Generate the plot
      # print(spatial_feature_plot(spatial_data = get(input$spatial_data),
      #                            type_data = {{input$data_type_visualization}},
      #                            peak_id = {{input$peak}},
      #                            # samples = c(samples_saline, get(input$experiment_samples)),
      #                            c(selected_saline_samples(), selected_experimental_samples()),
      #                            min_percentile = input$min_percentile,
      #                            max_percentile = input$max_percentile,
      #                            normalization = {{input$zero_normalize}},
      #                            size = input$spot_size,
      #                            show_legend = input$show_legend,
      #                            tif_image = input$tif_image) +
      #         plot_layout(ncol = as.numeric({{input$num_columns}}))) #+
      #         # plot_annotation(title = paste0(input$gene, ": ", input$peak))) 
      
      print(    spatial_feature_plot_cluster(spatial_data = get(input$spatial_data),
                                             type_data = {{input$data_type_visualization}},
                                             peak_id = {{input$peak}},
                                             clusters = vals$select_clusters,
                                             # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                                             # samples = c(samples_saline, get(input$experiment_samples)),
                                             # samples = c(selected_saline_samples(), selected_experimental_samples()),
                                             c(selected_saline_samples(), selected_experimental_samples()),
                                             min_percentile = input$min_percentile,
                                             max_percentile = input$max_percentile,
                                             normalization = {{input$zero_normalize}},
                                             resolution = vals$resolution,
                                             size = input$spot_size,
                                             show_legend = input$show_legend,
                                             tif_image = input$tif_image
      ) +
        plot_layout(ncol = as.numeric({{input$num_columns}})))
      
      # Close the png device
      dev.off()
    }
  )
  
  output$downloadSpatialFeaturePlotSVG <- downloadHandler(
    # # Define the DPI
    # dpi <- 96
    
    filename = function() {
      paste0(input$gene, '-', input$peak, '.svg')
    },
    content = function(file) {
      dpi <- 100  # You can adjust the DPI as needed
      # Open an SVG device, converting pixels to inches
      svg(filename = file, width = as.integer(input$width_plot) / dpi, height = as.integer(input$height_plot) / dpi)
      
      # Generate the plot
      # print(spatial_feature_plot(spatial_data = get(input$spatial_data),
      #                            type_data = {{input$data_type_visualization}},
      #                            peak_id = {{input$peak}},
      #                            # samples = c(samples_saline, get(input$experiment_samples)),
      #                            c(selected_saline_samples(), selected_experimental_samples()),
      #                            min_percentile = input$min_percentile,
      #                            max_percentile = input$max_percentile,
      #                            normalization = {{input$zero_normalize}},
      #                            size = input$spot_size,
      #                            show_legend = input$show_legend,
      #                            tif_image = input$tif_image) +
      #         plot_layout(ncol = as.numeric({{input$num_columns}})))# +
      #         # plot_annotation(title = paste0(input$gene, ": ", input$peak)))
      
      print(    spatial_feature_plot_cluster(spatial_data = get(input$spatial_data),
                                             type_data = {{input$data_type_visualization}},
                                             peak_id = {{input$peak}},
                                             clusters = vals$select_clusters,
                                             # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                                             # samples = c(samples_saline, get(input$experiment_samples)),
                                             # samples = c(selected_saline_samples(), selected_experimental_samples()),
                                             c(selected_saline_samples(), selected_experimental_samples()),
                                             min_percentile = input$min_percentile,
                                             max_percentile = input$max_percentile,
                                             normalization = {{input$zero_normalize}},
                                             resolution = vals$resolution,
                                             size = input$spot_size,
                                             show_legend = input$show_legend,
                                             tif_image = input$tif_image
      ) +
        plot_layout(ncol = as.numeric({{input$num_columns}})))
      
      # Close the SVG device
      dev.off()
    }
  )
  
  output$clusterPlot <- renderPlot({
    # visualize interest cluster
    spatial_cluster_select(
      spatial_data = get(input$spatial_data),
      resolution = vals$resolution,
      # samples = c(samples_saline, get(input$experiment_samples)),
      c(
        selected_saline_samples(),
        selected_experimental_samples()
      ),
      palette = palette_allen,
      size = input$spot_size,
      select_clusters = vals$select_clusters,
      tif_image = input$tif_image,
      ncol = 1
    ) +
      plot_layout(ncol = as.numeric({
        {
          input$num_columns
        }
      }))
  },
  width = function() {
    as.integer(input$width_plot)
  },  # cast the input values to integer
  height = function() {
    as.integer(input$height_plot)
  })
  
  
  # Download handler for downloading the cluster plot as a PNG image.
  output$downloadSpatialClusterPNG <- downloadHandler(
    filename = function() {
      # Set the filename of the downloaded file.
      # It is a .png file named 'cluster_plot'.
      paste0(input$gene, '-', input$peak, '-', input$data_type_visualization, '.svg')
    },
    content = function(file) {
      # Open a PNG device. The dimensions are specified in pixels.
      png(filename = file, width = as.integer(input$width_plot), height = as.integer(input$height_plot))
      
      # Generate the plot.
      # The 'print' function is used to draw the plot to the PNG device.
      print(
        spatial_cluster_select(
          spatial_data = get(input$spatial_data),
          resolution = vals$resolution,
          # samples = c(samples_saline, get(input$experiment_samples)),
          c(selected_saline_samples(), selected_experimental_samples()),
          palette = palette_allen,
          size = input$spot_size,
          select_clusters = vals$select_clusters,
          tif_image = input$tif_image,
          ncol = 1
        ) + 
          plot_layout(ncol = as.numeric({
            {
              input$num_columns
            }
          }))
      )
      
      # Close the PNG device. This is important because the file isn't actually written until the device is closed.
      dev.off()
    }
  )
  
  # Download handler for downloading the cluster plot as an SVG image.
  output$downloadSpatialClusterSVG <- downloadHandler(
    filename = function() {
      # Set the filename of the downloaded file.
      # It is a .svg file named 'cluster_plot'.
      paste0('cluster-plot-resolution', vals$resolution, '.svg')
    },
    content = function(file) {
      dpi <- 100  # You can adjust the DPI as needed
      # Open an SVG device. The dimensions are specified in inches, so we convert the pixel dimensions to inches by dividing by the DPI.
      svg(filename = file, width = (as.integer(input$width_plot) / dpi), height = ((as.integer(input$height_plot)) / dpi))
      
      # Generate the plot.
      # The 'print' function is used to draw the plot to the SVG device.
      print(
        spatial_cluster_select(
          spatial_data = get(input$spatial_data),
          resolution = vals$resolution,
          # samples = c(samples_saline, get(input$experiment_samples)),
          c(selected_saline_samples(), selected_experimental_samples()),
          palette = palette_allen,
          size = input$spot_size,
          select_clusters = vals$select_clusters,
          tif_image = input$tif_image,
          ncol = 1
        ) + 
          plot_layout(ncol = as.numeric({
            {
              input$num_columns
            }
          }))
      )
      
      # Close the SVG device. This is important because the file isn't actually written until the device is closed.
      dev.off()
    }
  )
  
  
  output$plot_umis <- renderPlot({
    plot_umi_reads_cluster(spatial_data = get(input$spatial_data), 
                           data_type = "raw_data",
                           resolution = {{vals$resolution}},
                           cluster = vals$cluster) -> p1
    
    plot_umi_reads_cluster(spatial_data = get(input$spatial_data),
                           data_type = {{input$data_type_visualization}},
                           resolution = {{vals$resolution}},
                           cluster = vals$cluster) -> p2
    
    # Combine the plots using the patchwork package
    p1/p2
  }, height = 800)
  
}
