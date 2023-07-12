generate_plots_and_zip <- function() {
  # Create a directory to store the png files
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  
  # Get a subset of the data that contains the genes of interest.
  # The selected dataframe is then stored in 'peaks_df'.
  clozapine_st_data_half$raw_data$annotate %>% 
    filter(gene_name %in% clozapine_genes_vector) %>% 
    select(c(peak_id, gene_name)) -> peaks_df
  
  # Setting up a text progress bar to track the progress of the loop.
  pb <- txtProgressBar(min = 0, max = nrow(peaks_df[1:10,]), style = 3)
  
  # Iterate over each row in the 'peaks_df' dataframe.
  # Generate a plot for each peak_id in peaks_df and save it as a png file.
  for(i in seq_len(nrow(peaks_df[1:10,]))) {
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
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar after the loop finishes
  close(pb)
  
  # Record the time before zipping files
  start_time <- Sys.time()
  
  # Zip all the generated png files into 'plots.zip'
  zip("plots.zip", files = "plots")
  
  # Record the time after zipping files
  end_time <- Sys.time()
  
  # Print out the time taken to zip the files
  print(end_time - start_time)
  
  # Delete the 'plots' directory and all its contents
  unlink("plots", recursive = TRUE)
}
