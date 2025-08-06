compute_global_expression_summary <- function(spatial_data,
                                              data_type = "raw_data",
                                              control_samples,
                                              experiment_samples,
                                              metric = c("mean", "sum"),
                                              trim = 0,
                                              verbose = FALSE) {
  expr_matrix <- spatial_data[[data_type]]$data
  metadata <- spatial_data[[data_type]]$metadata
  metadata$sample_barcode <- paste(metadata$sample, metadata$barcode, sep = "_")
  
  selected_samples <- c(control_samples, experiment_samples)
  selected_barcodes <- metadata %>%
    filter(sample %in% selected_samples) %>%
    pull(sample_barcode)
  
  expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% selected_barcodes, drop = FALSE]
  
  col_sample_map <- metadata %>%
    filter(sample_barcode %in% colnames(expr_matrix)) %>%
    select(sample_barcode, sample)
  
  sample_levels <- unique(col_sample_map$sample)
  
  if (verbose) {
    cat(">>> Computing metrics:", paste(metric, collapse = ", "), "\n")
  }
  
  results_list <- list()
  
  for (m in metric) {
    if (verbose) {
      cat(">>", m, "\n")
      pb <- txtProgressBar(min = 0, max = length(sample_levels), style = 3)
    }
    
    metric_df <- data.frame(
      peak = spatial_data[[data_type]]$annotate$peak_id,
      gene = spatial_data[[data_type]]$annotate$gene_name,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(sample_levels)) {
      s <- sample_levels[i]
      barcodes <- col_sample_map %>%
        filter(sample == s) %>%
        pull(sample_barcode)
      
      mat <- expr_matrix[, barcodes, drop = FALSE]
      
      values <- switch(m,
                       sum        = rowSums(mat, na.rm = TRUE),
                       mean       = apply(mat, 1, function(x) mean(x, trim = trim, na.rm = TRUE)),
                       median     = apply(mat, 1, function(x) median(x, na.rm = TRUE)),
                       var        = apply(mat, 1, function(x) var(x, na.rm = TRUE)),
                       IQR        = apply(mat, 1, function(x) IQR(x, na.rm = TRUE)),
                       diff_range = apply(mat, 1, function(x) diff(range(x, na.rm = TRUE))),
                       skewness   = apply(mat, 1, function(x) e1071::skewness(x, na.rm = TRUE)),
                       kurtosis   = apply(mat, 1, function(x) e1071::kurtosis(x, na.rm = TRUE)),
                       stop("Unsupported metric: ", m)
      )
      
      metric_df[[s]] <- values
      
      if (verbose) setTxtProgressBar(pb, i)
    }
    
    if (verbose) close(pb)
    
    results_list[[m]] <- metric_df
  }
  
  return(results_list)
}

compute_global_expression_multi <- function(spatial_data,
                                            data_types,
                                            control_samples,
                                            experiment_samples,
                                            metrics = c("mean", "sum"),
                                            trim = 0,
                                            verbose = TRUE) {
  #' Wrapper to compute global expression summaries across multiple data types
  #'
  #' @param spatial_data List with assays in spatial format (must contain named entries matching `data_types`)
  #' @param data_types Character vector with names of data_type slots to iterate over
  #' @param control_samples Character vector of control sample IDs
  #' @param experiment_samples Character vector of experimental sample IDs
  #' @param metrics Character vector of metrics to compute (e.g., "mean", "sum")
  #' @param trim Trim fraction for computing trimmed means (default 0)
  #' @param verbose Print progress
  #'
  #' @return Named list with structure:
  #'   result[[data_type]][[metric]] -> dataframe of peak/gene vs sample
  
  result_list <- list()
  
  for (dt in data_types) {
    if (verbose) {
      cat("\n=====================================\n")
      cat(">>> Processing data_type:", dt, "\n")
      cat("=====================================\n")
    }
    
    summary_res <- compute_global_expression_summary(
      spatial_data       = spatial_data,
      data_type          = dt,
      control_samples    = control_samples,
      experiment_samples = experiment_samples,
      metric             = metrics,
      trim               = trim,
      verbose            = verbose
    )
    
    result_list[[dt]] <- summary_res
  }
  
  return(result_list)
}
