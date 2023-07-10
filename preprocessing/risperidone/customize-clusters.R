
nfeatures <- 2000
dims <- 1:30


sample_names <- list.files(path = path_to_data)





# read images to spatial transcriptoms for risperidone
images_risperidone_half <- create_images_tibble(path_to_data, sample_names)

info_peaks_risperidone <- read_gene_annotation(file_path = "data/risperidone/gene-annotation/peaks-annotate-reduction.tsv")

# read barcode data for risperidone
barcode_risperidone_half <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_risperidone_half
)


nfeatures <- 2000
dims <- 1:30


risperidone_integrate_half <- 
  tmp <- integrate_data_seurat(path_to_data, nfeatures, dims)



tmp <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_risperidone,
                                                barcode_info = barcode_risperidone_half,
                                                images_info = images_risperidone_half,
                                                integrated_data = risperidone_integrate_half,
                                                peaks_info = info_peaks_risperidone)

tmp <- add_clusters_data(spatial_data = risperidone_st_data_half,
                                              integrated_data = risperidone_integrate_half,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

tmp <- evaluate_clustering_stability(spatial_data = risperidone_st_data_half,
                                                          seurat_object = risperidone_integrate_half,
                                                          resolution_start = 0.05,
                                                          resolution_end = 2,
                                                          resolution_step = 0.05)



process_spatial_data <- function(dims, nfeatures) {
  # Input: dims and nfeatures
  
  # Measure start time
  start_time <- Sys.time()
  
  # Code block
  integrated_data <- integrate_data_seurat(path_to_data, nfeatures, dims)
  
  tmp <- create_spatial_data(sample_names = sample_names,
                             metadata = metadata_risperidone,
                             barcode_info = barcode_risperidone_half,
                             images_info = images_risperidone_half,
                             integrated_data = integrated_data,
                             peaks_info = info_peaks_risperidone)
  
  tmp <- add_clusters_data(spatial_data = tmp,
                           integrated_data = integrated_data,
                           resolution_start = 0.05,
                           resolution_end = 2,
                           resolution_step = 0.05)
  
  tmp <- evaluate_clustering_stability(spatial_data = tmp,
                                       seurat_object = integrated_data,
                                       resolution_start = 0.05,
                                       resolution_end = 1,
                                       resolution_step = 0.05)
  
  # Measure end time
  end_time <- Sys.time()
  
  # Calculate time taken
  time_taken <- end_time - start_time
  cat("Time taken:", time_taken, "seconds\n")
  
  # Return the result
  return(tmp)
}

# Example usage:
nfeatures <- 2000
dims <- 30

nfeatures_vect <- c(500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000)
dims_vect <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)

cluster_results_ris <- list()

for(nfeatures in nfeatures_vect){
  for(dims in dims_vect){
    element_name <- paste0("nfeatures", nfeatures, "_dims", dims)
    
    result <- process_spatial_data(1:dims, nfeatures)
    
    cluster_results_ris[[element_name]][["clusters"]] <- result$clusters
    cluster_results_ris[[element_name]][["stability_results"]] <- result$stability_results
  }
}


cluster_results_ris %>%
  # Add names as a column
  purrr::imap(~ mutate(.x$stability_results, element_name = .y)) %>%
  # Bind all the data frames into a single data frame
  bind_rows() %>% filter(resolution == 0.2) %>% filter(num_clusters == 12) %>% filter(silhouette_score> 0.4)


risperidone_st_data_half -> tmp


tmp$clusters <- cluster_results_ris$nfeatures2500_dims40$clusters

spatial_cluster(spatial_data = tmp,
                resolution = 0.2,
                samples = c(samples_saline, samples_risperidone),
                palette = palette_allen,
                size= 1.0,
                ncol = 4)

spatial_cluster(spatial_data = risperidone_st_data_half,
                resolution = 0.2,
                samples = c(samples_saline, samples_risperidone),
                palette = palette_allen,
                size= 1.0,
                ncol = 4)
