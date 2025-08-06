compute_data_summary_v3 <- function(spatial_data,
                                 resolution = 0.8,
                                 trim = 0.05,
                                 num_cores = 4,
                                 control_samples,
                                 experiment_samples,
                                 data_type = "raw_data",
                                 min_number_spots = 20,
                                 metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis")) {
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort()
  
  compute_data_summary_cluster <- function(cluster) {
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_"),
             group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
      select(group, sample_barcode)
    
    barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
    control_barcodes <- barcodes_list[["control"]]
    experiment_barcodes <- barcodes_list[["experiment"]]
    
    control_expression_spot <- sapply(control_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)], drop = FALSE]
    })
    
    experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)], drop = FALSE]
    })
    
    compute_metrics <- function(expr_list, metrics) {
      result <- list()
      if ("expression_spot" %in% metrics) result$expression_spot <- expr_list
      if ("sum" %in% metrics) result$sum <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, sum(y))))
      if ("mean" %in% metrics) result$mean <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, mean(y, trim = trim, na.rm = TRUE))))
      if ("median" %in% metrics) result$median <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, median(y))))
      if ("IQR" %in% metrics) result$IQR <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, IQR(y))))
      if ("diff_range" %in% metrics) result$diff_range <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, diff(range(y)))))
      if ("var" %in% metrics) result$var <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, var(y))))
      if ("skewness" %in% metrics) result$skewness <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, skewness(y)))))
if ("kurtosis" %in% metrics) result$kurtosis <- sapply(expr_list, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, kurtosis(y)))))
return(result)
    }
    
    control_list <- compute_metrics(control_expression_spot, metrics)
    experiment_list <- compute_metrics(experiment_expression_spot, metrics)
    
    return(list(
      peak = spatial_data[[data_type]]$annotate$peak_id,
      gene = spatial_data[[data_type]]$annotate$gene_name,
      control = control_list,
      experiment = experiment_list
    ))
  }
  
  # Setup parallel backend
  plan(multisession, workers = num_cores)
  
  start_time <- Sys.time()
  
  results <- future_map(unique_clusters, ~try(compute_data_summary_cluster(.x), silent = TRUE))
  
  names(results) <- paste0("cluster_", unique_clusters)
  
  # Optionally filter out errors
  results <- results[!sapply(results, inherits, "try-error")]
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(results)
}