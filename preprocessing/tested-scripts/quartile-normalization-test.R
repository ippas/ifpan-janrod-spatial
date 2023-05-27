BiocManager::install("preprocessCore")

library(preprocessCore)


add_quantile_normalization <- function(spatial_data, data_type = "raw_data", resolution = 0.8) {
  # Extract barcodes for both control and experiment groups
  # Join metadata and cluster information, and filter by cluster
  # A new column 'sample_barcode' is created by concatenating 'sample' and 'barcode'
  # Group label ('control' or 'experiment') is assigned based on sample
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
    pull(barcode)
  
  # Split barcodes into control and experiment groups based on group label
  barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
  
  # Extract control and experiment barcodes from the list
  control_barcodes <- barcodes_list[["control"]]
  experiment_barcodes <- barcodes_list[["experiment"]]
  
  # Prepare data for the statistical test depending on expression_unit
  # Preparing data for statistical tests
  if (expression_unit == "spot") {
    # Direct retrieval of expression data when unit is "spot"
    control_expression <- spatial_data[[data_type]]$data[, control_barcodes]
    experiment_expression <- spatial_data[[data_type]]$data[, experiment_barcodes]
}}


quantile_norm_cluster <- function(cluster){
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    pull(sample_barcode)
  
  cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] %>% as.matrix()
  
  norm_data <- normalize.quantiles(cluster_data)

  colnames(norm_data) <- colnames(cluster_data)
  rownames(norm_data) <- rownames(cluster_data)

  return(norm_data)
  
} 

start_time <- Sys.time()
tmp <- mclapply(unique_clusters, quantile_norm_cluster,  mc.cores = num_cores)
# Set names for the results
tmp <- setNames(tmp, paste0("cluster_", unique_clusters))
end_time <- Sys.time()
diff_time <- start_time - end_time

print(diff_time)

add_quantile_norm_data <- function(spatial_data, resolution = 0.8, num_cores = 24, data_type = "raw_data"){
  
  quantile_norm_cluster <- function(cluster){
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      pull(sample_barcode)
    
    cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] %>% as.matrix()
    
    norm_data <- normalize.quantiles(cluster_data)
    
    colnames(norm_data) <- colnames(cluster_data)
    rownames(norm_data) <- rownames(cluster_data)
    
    return(norm_data)
    
  } 
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  start_time <- Sys.time()
  results <- mclapply(unique_clusters, quantile_norm_cluster,  mc.cores = num_cores)
  # Set names for the results
  results <- setNames(results, paste0("cluster_", unique_clusters))
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  
  print(diff_time)
  
  do.call(cbind, results) %>% Matrix() -> quantile_normalize_data
  
  # Keep only the filtered metadata
  spatial_data$quantile_normalize$metadata <- spatial_data$raw_data$metadata
  
  # Keep only the filtered annotation data
  spatial_data$quantile_normalize$annotate <- spatial_data$raw_data$annotate
  
  # Keep only the filtered expression data
  spatial_data$quantile_normalize$data <- quantile_normalize_data
  
  return(spatial_data)
}


quantile_norm_cluster(cluster = 1) 

add_quantile_norm_data(spatial_data = risperidone_st_data_half,
                       resolution = 0.8,
                       num_cores = 24,
                       data_type = "raw_data") -> tmp2 



spatial_gene_plot(spatial_data = tmp,
                  type_data = "raw_data",
                  gene = "Egr3",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)



















# old data

cbind(results2$cluster_0$control_expression,
  results2$cluster_0$experiment_expression) -> data_matrix

# Save the rownames and colnames before normalization
rownames_original <- rownames(data_matrix)
colnames_original <- colnames(data_matrix)

# Perform quantile normalization
data_matrix_normalized <- normalize.quantiles(data_matrix)

# Restore the rownames and colnames
rownames(data_matrix_normalized) <- rownames_original
colnames(data_matrix_normalized) <- colnames_original
data_matrix_normalized %>% head
data_matrix %>% head
cluster <- 2

# Add filtered data to the spatial_data object
add_quartilenorm_data <- function(cluster, trim = 0.05) {
  
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    # Assign 'control' or 'experiment' group label based on sample
    # mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
    pull(sample_barcode) 
  
  samples <-spatial_data$samples

  # calculate mean for each sample
  mean_expression <- sapply(samples, function(sample_id) {
    samples_data <- spatial_data[[data_type]]$data[, all_barcodes[grepl(paste0(sample_id, "_"), all_barcodes)]]

    apply(samples_data, 1, mean, trim = trim)

  })
  
  ### prepare normalized data for mean
  # Save the rownames and colnames before normalization
  rownames_original <- rownames(mean_expression)
  colnames_original <- colnames(mean_expression)
  
  # Perform quantile normalization
  norm_mean_expression <- normalize.quantiles(mean_expression)
  
  # Restore the rownames and colnames
  rownames(norm_mean_expression) <- rownames_original
  colnames(norm_mean_expression) <- colnames_original
  
  # calculation of coefficients for multiplication of raw data
  norm_mean_by_mean <- norm_mean_expression/mean_expression  

  ### normalize raw data
  norm_raw_expression <- sapply(samples, function(sample_id) {
    samples_data <- spatial_data[[data_type]]$data[, all_barcodes[grepl(paste0(sample_id, "_"), all_barcodes)]]
    norm_vector <- norm_mean_by_mean[, sample_id]
    norm_raw_data <- samples_data * norm_vector
  })
  
  norm_cluster <- purrr::reduce(norm_raw_expression, cbind)
  
  return(norm_cluster)
}

data_type <- "raw_data"

start_time <- Sys.time()
# Use mclapply instead of lapply
raw_normalize <- mclapply(unique_clusters, add_quartilenorm_data, mc.cores = num_cores)
# Set names for the results
raw_normalize <- setNames(raw_normalize, paste0("cluster_", unique_clusters))
end_time <- Sys.time()
end_time - start_time

purrr::reduce(raw_normalize, cbind) -> tmp

tmp %>% class

as.matrix(tmp) %>% class

as.matrix(tmp) -> tmp

tmp[is.nan(tmp)] <- 0

spatial_data$raw_data$metadata
tmp_risperidone_st_data_half <- spatial_data

tmp_risperidone_st_data_half$norm_quartile$data <- tmp
tmp_risperidone_st_data_half$norm_quartile$annotate <- spatial_data$raw_data$annotate
tmp_risperidone_st_data_half$norm_quartile$metadata <- spatial_data$raw_data$metadata

spatial_gene_plot(spatial_data = tmp_risperidone_st_data_half,
                  type_data = "norm_quartile",
                  gene = "Socs5",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = F)

  

