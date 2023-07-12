# # Load the libraries
# install.packages("shinydashboard")
# install.packages("shiny")


source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")
require(shiny)
require(shinydashboard)

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))

# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]


load("results/risperidone/risperidone.RData")
load("results/pz1190/pz1190.RData")

library(future)
plan(multisession) 
# Define the user interface
ui <- fluidPage(
  titlePanel("Gene and Peak Selection"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("drug", "Select Drug", choices = c("Risperidone", "PZ-1190", "Clozapine")),
      
      selectInput("spatial_data", "Select spatial data", choices = NULL),
      
      HTML("<h5><b>Parameters for choosing data</b></h5>"),
      
      selectInput("summary_data", "Select summary data", choices = NULL),
      
      selectInput("experiment_samples", "Select experiment samples", choices = NULL),
      
      selectInput("metric",
                  "Select metric",
                  choices = c("mean", 
                              "median", "sum")),
      
      radioButtons("data_type_statistics", 
                   "Choose data type for statistics:", 
                   choices = c("raw_data", 
                               "seurat",
                               "range_normalize",
                               "quantile_metric",
                               "quantile_normalize")),
      
      radioButtons("resolution_statistics", 
                   "Choose resolution for statistics:", 
                   choices = c("0.05",
                               "0.1",
                               "0.15",
                               "0.2",
                               "0.4",
                               "0.8")),
      
      HTML("<h5><b> Parameters for filtering statistics </b></h5>"), 
      
      # Custom CSS to restrict resizing to vertical only
      tags$head(tags$style(HTML("#gene_names { resize: vertical;   height: 355px;  }"))),
      
      # UI layout
      fluidRow(
        # Left column with parameters
        column(6,
               numericInput("control_mean_threshold", "Control Mean Threshold:", value = 0.2, step = 0.1),
               numericInput("experiment_mean_threshold", "Experiment Mean Threshold:", value = 0.2, step = 0.1),
               numericInput("log2ratio_threshold", "Log2Ratio Threshold:", value = 0.5, step = 0.1),
               numericInput("t_test_threshold", "T-Student Test Threshold:", value = 0.05, step = 0.1),
               numericInput("wilcoxon_test_threshold", "Wilcoxon Test Threshold:", value = 0.05, step = 0.1),
               numericInput("ks_test_threshold", "Kolmogorov-Smirnov Test Threshold:", value = 0.05, step = 0.1)
        ),
        
        # Right column with gene names input and submit button
        column(6,
               textAreaInput("gene_names", "Enter Gene Names:", rows = 5),
               # column(12, progressBar("calcProgress", value = 0, color = "success"))
               fluidRow(
                 column(6, actionButton("download_plots", "Download plots", width = "100%")),
                 column(6, actionButton("abort_download", "Abort download", width = "100%"))
               ),
               downloadButton("downloadFile", "Download Zip File"),
               actionButton("clear_genes", "Clear Gene Names", rows = 5)
        )
      ),
      
      HTML("<h5><b> Visualization parameters </b></h5>"), 
      
      uiOutput("data_type_visualization"),
      
      checkboxInput("zero_normalize", 
                    "Zero Normalize", 
                    value = FALSE),
      
      checkboxInput("tif_image", 
                    "Background image", 
                    value = TRUE),
      
      fluidRow(
        column(6, uiOutput("gene")),
        column(6, selectInput("peak", "Select Peak:", choices = NULL))
      ),
      
      fluidRow(
        column(6, selectInput("num_columns", "Number of Columns to Display:", choices = c(4, 2))),
        column(6, numericInput("spot_size", "Spot size:", 1, min = 0.1, max = 5, step = 0.1)),
      ),
      
      fluidRow(
        column(6, numericInput("width_plot", "Width:", 1100, min = 100, step = 10)),
        column(6, numericInput("height_plot", "Height:", 750, min = 100, step = 10))
      ),
      
      fluidRow(
        column(6, 
               numericInput("min_percentile", 
                            "Input min percentile to remove:", 
                            value = 0, min = 0, max = 1, step = 0.001)
        ),
        column(6, 
               numericInput("max_percentile", 
                            "Input max percentile to remove 2:", 
                            value = 1, min = 0, max = 1, step = 0.001)
        )
      ),
      
      sliderInput("cluster", "Select number cluster:", min = 0, max = 25, value = 0),
      sliderInput("resolution", "Select resolution:", min = 0.05, max = 2, step = 0.05, value = 0.05),
      actionButton("update_values", "Update Resolution and Cluster")
      
    ),
    
    mainPanel(
      
      # tableOutput("filter_statistics"),
      dataTableOutput("filter_statistics"),
      
      # # Add this to your UI layout
      fluidRow(
        column(4, 
               div(style="display: inline-flex; align-items: center;", 
                   tags$span("Enter Filename:", style="font-weight: bold; margin-right: 10px;"),
                   textInput("filename", label = NULL, value = paste("results_", Sys.Date(), ".tsv", sep = ""))
               )
        ),
        column(4, downloadButton('downloadData', 'Download TSV'))
      ),
      
      verbatimTextOutput("selected_info"),
      
      plotOutput("peakPlot", height = "100%", width = "100%"),  
      
      plotOutput("clusterPlot", height = "100%", width = "100%"),
      
      plotOutput("interestClusterPlot", height = "100%", width = "100%"),
      
      plotOutput("plot_umis", height = "100%", width = "100%")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  vals <- reactiveValues(cluster =  0, resolution = 0.05)
  
  observeEvent(input$update_values, {
    vals$cluster <- input$cluster
    vals$resolution <- input$resolution
  })
  
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
    spatial_feature_plot(spatial_data = get(input$spatial_data),
                         type_data = {{input$data_type_visualization}},
                         peak_id = {{input$peak}},
                         # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                         samples = c(samples_saline, get(input$experiment_samples)),
                         min_percentile = input$min_percentile,
                         max_percentile = input$max_percentile,
                         normalization = {{input$zero_normalize}},
                         size = input$spot_size,
                         tif_image = input$tif_image) +
      plot_layout(ncol = as.numeric({{input$num_columns}}))
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) })
  
  output$clusterPlot <- renderPlot({
    # visualize interest cluster
    spatial_cluster(spatial_data = get(input$spatial_data),
                    resolution = vals$resolution,
                    samples = c(samples_saline, get(input$experiment_samples)),
                    palette = palette_allen, 
                    size= 1.2, 
                    ncol = 4)
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) + 200})
  
  output$interestClusterPlot <- renderPlot({
    spatial_interest_cluster(cluster = vals$cluster,
                             spatial_data = get(input$spatial_data),
                             resolution = vals$resolution,
                             samples = c(samples_saline, get(input$experiment_samples)),
                             size= 1.3,
                             ncol = 4)
  },
  width = function() { as.integer(input$width_plot) },  # cast the input values to integer
  height = function() { as.integer(input$height_plot) + 200})
  
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
  
  
  # # Create a reactive expression that will trigger when the download button is clicked
  # observeEvent(input$download_plots, {
  #   # Generate the plots and zip them
  #   generate_plots_and_zip()
  # })
  # 
  # generate_plots_and_zip <- function() {
  #   
  #   # Create a directory to store the png files
  #   dir.create("plots")
  #   
  #   # Get a subset of the data that contains the genes of interest.
  #   # The selected dataframe is then stored in 'peaks_df'.
  #   clozapine_st_data_half$raw_data$annotate %>% 
  #     filter(gene_name %in% clozapine_genes_vector) %>% 
  #     select(c(peak_id, gene_name)) -> peaks_df
  #   
  #   shiny::withProgress(message = 'Generating plots...', value = 0, {
  #     
  #     for(i in seq_len(nrow(peaks_df[1:10,]))) {
  #       # Extract the peak and gene values from the current row
  #       peak <- peaks_df[i, "peak_id"]
  #       gene <- peaks_df[i, "gene_name"]
  #       
  #       # Generate the filename for the png file
  #       png_filename <- file.path("plots", paste0("plot_", gene, "-", peak, ".png"))
  #       
  #       # Initialize the png device with the generated filename
  #       png(filename = png_filename, width = 1100, height = 750)
  #       
  #       # Generate the spatial feature plot
  #       plot <- spatial_feature_plot(
  #         spatial_data = clozapine_st_data_half,
  #         type_data = "raw_data",
  #         peak_id = peak,
  #         samples = c(samples_saline, samples_clozapine),
  #         normalization = TRUE,
  #         tif_image = TRUE) +
  #         plot_layout(ncol = 4) + 
  #         plot_annotation(title = paste0(gene, ": ", peak))
  #       
  #       # Print the plot to the png file
  #       print(plot)
  #       
  #       # Close the png device
  #       dev.off()
  #       
  #       # Update the progress bar
  #       shiny::incProgress(i / nrow(peaks_df[1:10,]))
  #       
  #     }
  #   })
  #   
  #   # Close the progress bar after the loop finishes
  #   close(pb)
  #   
  #   # Zip all the generated png files into 'plots.zip'
  #   zip("plots.zip", files = "plots")
  #   
  #   # Delete the 'plots' directory and all its contents
  #   unlink("plots", recursive = TRUE)
  # }
  # 
  # 
  # # Enable download of the zip file
  # output$downloadFile <- downloadHandler(
  #   filename = "plots.zip",
  #   content = function(file) {
  #     file.copy("plots.zip", file)
  #   }
  # )
  
  # # Create a reactiveValues object to control the loop
  # control <- reactiveValues(abort = FALSE)
  # 
  # # Trigger abort when abort button is clicked
  # observeEvent(input$abort_download, {
  #   control$abort <- TRUE
  # })
  # 
  # observeEvent(input$download_plots, {
  #   # Reset abort flag at the beginning of process
  #   control$abort <- FALSE
  #   generate_plots_and_zip()
  # })
  # 
  # generate_plots_and_zip <- function() {
  #   
  #   # Create a directory to store the png files
  #   dir.create("plots")
  #   
  #   # Get a subset of the data that contains the genes of interest.
  #   # The selected dataframe is then stored in 'peaks_df'.
  #   clozapine_st_data_half$raw_data$annotate %>% 
  #     filter(gene_name %in% clozapine_genes_vector) %>% 
  #     select(c(peak_id, gene_name)) -> peaks_df
  #   
  #   shiny::withProgress(message = 'Generating plots...', value = 0, {
  #     for(i in seq_len(nrow(peaks_df[1:10,]))) {
  #       
  #       # Check if abort button was pressed
  #       if (control$abort) {
  #         stop("Aborted by user")  # or return(), depending on what you want to do
  #       }
  #       
  #       # Extract the peak and gene values from the current row
  #       peak <- peaks_df[i, "peak_id"]
  #       gene <- peaks_df[i, "gene_name"]
  #       
  #       # Generate the filename for the png file
  #       png_filename <- file.path("plots", paste0("plot_", gene, "-", peak, ".png"))
  #       
  #       # Initialize the png device with the generated filename
  #       png(filename = png_filename, width = 1100, height = 750)
  #       
  #       # Generate the spatial feature plot
  #       plot <- spatial_feature_plot(
  #         spatial_data = clozapine_st_data_half,
  #         type_data = "raw_data",
  #         peak_id = peak,
  #         samples = c(samples_saline, samples_clozapine),
  #         normalization = TRUE,
  #         tif_image = TRUE) +
  #         plot_layout(ncol = 4) + 
  #         plot_annotation(title = paste0(gene, ": ", peak))
  #       
  #       # Print the plot to the png file
  #       print(plot)
  #       
  #       # Close the png device
  #       dev.off()
  #       
  #       # Update the progress bar
  #       shiny::incProgress(i / nrow(peaks_df[1:10,]))
  #     }
  #   })
  #   
  #   # Close the progress bar after the loop finishes
  #   close(pb)
  #   
  #   # Zip all the generated png files into 'plots.zip'
  #   zip("plots.zip", files = "plots")
  #   
  #   # Delete the 'plots' directory and all its contents
  #   unlink("plots", recursive = TRUE)
  # }
  # 
  # 
  # # Enable download of the zip file
  # output$downloadFile <- downloadHandler(
  #   filename = "plots.zip",
  #   content = function(file) {
  #     file.copy("plots.zip", file)
  #   }
  # )
  # This will hold the future object
  # This will hold the future object
  progress <- reactiveVal(0)
  should_abort <- reactiveVal(FALSE)
  
  # Create a reactive expression that will trigger when the download button is clicked
  observeEvent(input$download_plots, {
    should_abort(FALSE)  # Reset abort flag
    progress(0)  # Reset progress
    generate_plots_and_zip()
  })
  
  generate_plots_and_zip <- function() {
    # Create a directory to store the png files
    dir.create("plots")
    
    # Get a subset of the data that contains the genes of interest.
    clozapine_st_data_half$raw_data$annotate %>% 
      filter(gene_name %in% clozapine_genes_vector) %>% 
      select(c(peak_id, gene_name)) -> peaks_df
    
    for(i in seq_len(nrow(peaks_df[1:10,]))) {
      if (should_abort() == TRUE) {
        print("Download aborted")
        unlink("plots", recursive = TRUE)  # Delete the 'plots' directory and all its contents
        return()
      }
      
      # Extract the peak and gene values from the current row
      peak <- peaks_df[i, "peak_id"]
      gene <- peaks_df[i, "gene_name"]
      
      # Generate the filename for the png file
      png_filename <- file.path("plots", paste0("plot_", gene, "-", peak, ".png"))
      
      # Initialize the png device with the generated filename
      png(filename = png_filename, width = 1100, height = 750)
      
      # Generate the spatial feature plot
      plot <- spatial_feature_plot(
        spatial_data = clozapine_st_data_half,
        type_data = "raw_data",
        peak_id = peak,
        samples = c(samples_saline, samples_clozapine),
        normalization = TRUE,
        tif_image = TRUE) +
        plot_layout(ncol = 4) + 
        plot_annotation(title = paste0(gene, ": ", peak))
      
      # Print the plot to the png file
      print(plot)
      
      # Close the png device
      dev.off()
      
      # Update the progress
      progress(i / nrow(peaks_df[1:10,]) * 100)
    }
    
    # Zip all the generated png files into 'plots.zip'
    zip("plots.zip", files = "plots")
    
    # Delete the 'plots' directory and all its contents
    unlink("plots", recursive = TRUE)
  }
  
  observeEvent(input$abort_download, {
    should_abort(TRUE)
  })
  
  # Enable download of the zip file
  output$downloadFile <- downloadHandler(
    filename = "plots.zip",
    content = function(file) {
      file.copy("plots.zip", file)
    }
  )
  
  # Watch the progress reactiveValues and update the progressBar
  observe({
    shinyWidgets::updateProgressBar(session, "plot_progress", value = progress())
  })
  
  
}


# Run the app
shinyApp(ui = ui, server = server)

