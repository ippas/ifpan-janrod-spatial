# Calculate gene expression statistics between control and experiment groups
# Args:
#   spatial_data: A list containing the spatial transcriptomic data and associated metadata
#   stat_test: Statistical test to use for comparison (e.g. "t.test", "log2ratio", "wilcox.test")
#   type_data: The type of data to use (e.g. "normalized", "scaled")
#   resolution: The resolution at which to define clusters
#   per: The level of aggregation (either "spot" or "sample")
#   samples_vector1: Vector of sample names for group 1 (control)
#   samples_vector2: Vector of sample names for group 2 (experiment)
#   save_spatial_data: Whether to save the resulting statistics in the spatial_data object (default is FALSE)
# Calculate gene expression statistics between control and experiment groups
calculate_gene_expression_stats <- function(spatial_data, 
                                            stat_test,
                                            data_type,
                                            resolution,
                                            expression_unit,
                                            control_samples,
                                            experiment_samples,
                                            save_results = FALSE) {
  
  # Create a column name for the specified resolution
  resolution_column <- paste("cluster_resolution", resolution, sep = "_")
  
  # Get the unique cluster labels
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  # Initialize a list to store p-values for each cluster
  pvalue_list <- list()
  
  
  
  # Iterate through each unique cluster
  for (cluster in unique_clusters[0:21]) {
    print(cluster)
    
    # Get the barcodes for the control group in the current cluster
    control_barcodes <- spatial_data[[data_type]]$metadata %>% 
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% control_samples) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      pull(sample_barcode)
    
    # Get the barcodes for the experiment group in the current cluster
    experiment_barcodes <- spatial_data[[data_type]]$metadata %>% 
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% experiment_samples) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      pull(sample_barcode)
 
    # Prepare data for the statistical test
    if (expression_unit == "spot") {
      control_expression <- spatial_data[[data_type]]$data[, control_barcodes]
      experiment_expression <- spatial_data[[data_type]]$data[, experiment_barcodes]
    } else if (expression_unit == "sample") {
      control_expression <- sapply(control_samples, function(sample_id) {
        control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(sample_id, control_barcodes)]]
        apply(control_data, 1, mean, trim = 0.05)
      })
      
      experiment_expression <- sapply(experiment_samples, function(sample_id) {
        experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(sample_id, experiment_barcodes)]]
        apply(experiment_data, 1, mean, trim = 0.05)
      })
    }
    

    
    result_vector <- numeric(nrow(spatial_data[[data_type]]$data))
    
    # Perform the specified statistical test
    for (index in 1:nrow(spatial_data[[data_type]]$data)) {
      control_values <- control_expression[index, ]
      experiment_values <- experiment_expression[index, ]
      
      if (stat_test == "t.test") {
        result_vector[index] <- t.test(control_expression[index, ], experiment_expression[index, ], var.equal = TRUE)$p.value
      } else if (stat_test == "log2ratio") {
        result_vector[index] <- log2_ratio_result(control_values, experiment_values)
      } else if (stat_test == "wilcox.test") {
        result_vector[index] <- wilcox_test_result(control_values, experiment_values)
      }
    }
    
    pvalue_list[[paste("cluster", cluster, sep = "_")]] <- result_vector
  }
  
  result_matrix <- pvalue_list %>% do.call(rbind, .) %>% t
  rownames(result_matrix) <- rownames(spatial_data[[data_type]]$data)
  
  if (save_results) {
    spatial_data[[data_type]]$statistics[[resolution_column]][[expression_unit]][[stat_test]] <- result_matrix
    spatial_data
  }
  else if(save_results == FALSE){
    result_matrix
  }
  
}
    
    