#

zero_nonzero_matrix <- function(cluster) {
  # Extract barcodes for both control and experiment groups
  # Join metadata and cluster information, and filter by cluster
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    # Assign 'control' or 'experiment' group label based on sample
    mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
    select(group, sample_barcode)
  
  # Split barcodes into control and experiment groups based on group label
  barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
  
  # Extract control and experiment barcodes from the list
  control_barcodes <- barcodes_list[["control"]]
  experiment_barcodes <- barcodes_list[["experiment"]]
  
  control_expression <- spatial_data[[data_type]]$data[, control_barcodes]
  experiment_expression <- spatial_data[[data_type]]$data[, experiment_barcodes]
  
  # Define a function to count the number of zero and non-zero values in a vector
  zero_count <- function(vec) {
    return(sum(vec == 0))
  }
  
  non_zero_count <- function(vec) {
    return(sum(vec != 0))
  }
  
  # Calculate zero count for each sample in control group
  control_zero_count <- sapply(control_samples, function(sample_id) {
    control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
    apply(control_data, 1, zero_count)
  })
  # Calculate non-zero count for each sample in control group
  control_non_zero_count <- sapply(control_samples, function(sample_id) {
    control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
    apply(control_data, 1, non_zero_count)
  })
  
  # Calculate zero count for each sample in experiment group
  experiment_zero_count <- sapply(experiment_samples, function(sample_id) {
    experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
    apply(experiment_data, 1, zero_count)
  })
  
  # Calculate non-zero count for each sample in experiment group
  experiment_non_zero_count <- sapply(experiment_samples, function(sample_id) {
    experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
    apply(experiment_data, 1, non_zero_count)
  })
  
  control_zero_nonzero_ratio <- control_zero_count/control_non_zero_count
  
  experiment_zero_nonzero_ratio <- experiment_zero_count/experiment_non_zero_count
  
  # Return a list containing control and experiment barcodes, their expressions, and the p-values
  return(list(control_barcodes = control_barcodes, 
              experiment_barcodes = experiment_barcodes,
              control_expression = control_expression,
              experiment_expression = experiment_expression,
              peaks = rownames(control_expression),
              genes = spatial_data[[data_type]]$annotate$gene_name,
              control_zero_count = control_zero_count,
              control_non_zero_count = control_non_zero_count,
              experiment_zero_count = experiment_zero_count,
              experiment_non_zero_count = experiment_non_zero_count,
              control_zero_nonzero_ratio = control_zero_nonzero_ratio,
              experiment_zero_nonzero_ratio = experiment_zero_nonzero_ratio
  ))
}


data_type <- "raw_data"

start_time <- Sys.time()
# Use mclapply instead of lapply
results_zero_nonzero <- mclapply(unique_clusters, zero_nonzero_matrix)
# Set names for the results
results_zero_nonzero <- setNames(results_zero_nonzero, paste0("cluster_", unique_clusters))
end_time <- Sys.time()
end_time - start_time

results_zero_nonzero$cluster_0$control_zero_nonzero_ratio





######################################3
# 2. version
zero_nonzero_matrix <- function(cluster, min_spots = 30) {
  # Join metadata and cluster information, and filter by cluster
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
    select(group, sample_barcode)
  
  # Split barcodes into control and experiment groups
  barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
  control_barcodes <- barcodes_list[["control"]]
  experiment_barcodes <- barcodes_list[["experiment"]]
  
  # Extract control and experiment expression data
  control_expression <- spatial_data[[data_type]]$data[, control_barcodes]
  experiment_expression <- spatial_data[[data_type]]$data[, experiment_barcodes]
  
  # Functions to count zero and non-zero values
  zero_count <- function(vec) sum(vec == 0)
  non_zero_count <- function(vec) sum(vec != 0)
  
  # Helper function to apply the zero_count and non_zero_count functions to data
  apply_counts <- function(barcodes, samples, count_func) {
    sapply(samples, function(sample_id) {
      data <- spatial_data[[data_type]]$data[, barcodes[grepl(paste0(sample_id, "_"), barcodes)]]
      # Check if the number of spots per sample is less than min_spots
      if(ncol(data) < min_spots){
        return(rep(NA, nrow(data)))
      }else{
        return(apply(data, 1, count_func))
      }
      # Set rownames
      rownames(result) <- rownames(spatial_data[[data_type]]$data)
      return(result)
    })
  }
  
  # Calculate zero and non-zero counts
  control_zero_count <- apply_counts(control_barcodes, control_samples, zero_count)
  control_non_zero_count <- apply_counts(control_barcodes, control_samples, non_zero_count)
  experiment_zero_count <- apply_counts(experiment_barcodes, experiment_samples, zero_count)
  experiment_non_zero_count <- apply_counts(experiment_barcodes, experiment_samples, non_zero_count)
  
  # Calculate zero/non-zero ratios
  control_zero_nonzero_ratio <- control_zero_count / control_non_zero_count
  experiment_zero_nonzero_ratio <- experiment_zero_count / experiment_non_zero_count
  
  # Return results as a list
  list(
    control_barcodes = control_barcodes, 
    experiment_barcodes = experiment_barcodes,
    control_expression = control_expression,
    experiment_expression = experiment_expression,
    peaks = rownames(control_expression),
    genes = spatial_data[[data_type]]$annotate$gene_name,
    control_zero_count = control_zero_count,
    control_non_zero_count = control_non_zero_count,
    experiment_zero_count = experiment_zero_count,
    experiment_non_zero_count = experiment_non_zero_count,
    control_zero_nonzero_ratio = control_zero_nonzero_ratio,
    experiment_zero_nonzero_ratio = experiment_zero_nonzero_ratio
  )
}

start_time <- Sys.time()
# Use mclapply instead of lapply
results_zero_nonzero <- mclapply(unique_clusters, zero_nonzero_matrix, min_spots = 15)
# Set names for the results
results_zero_nonzero <- setNames(results_zero_nonzero, paste0("cluster_", unique_clusters))
end_time <- Sys.time()
end_time - start_time

results_zero_nonzero$cluster_21$control_zero_count
