summarize_and_test_v2 <- function(spatial_data,
                                  data_params_df,
                                  control_samples,
                                  experiment_samples,
                                  summary_metrics = c("mean", "median"),
                                  statistic_metrics = c("mean", "median"),
                                  mean_threshold = 0,
                                  trim = 0.05,
                                  min_number_spots = 20,   # <-- NOWY ARGUMENT
                                  num_cores = 20,
                                  verbose = TRUE) {
  
  result_summary_statistics <- list()
  n_steps <- nrow(data_params_df)
  
  pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  
  for (i in seq_len(n_steps)) {
    data_type  <- data_params_df[i, 1, drop = TRUE]
    resolution <- data_params_df[i, 2, drop = TRUE]
    data_type_name <- data_params_df[i, 3, drop = TRUE]
    quantile_normalization <- data_params_df[i, 4, drop = TRUE]
    
    if (verbose) {
      cat("\n=============================\n")
      cat(">>> ROW", i, "/", n_steps, ": ", data_type, "@", resolution, "\n")
      cat("=============================\n")
    }
    
    # ---- SUMMARY ----
    summary_data <- compute_data_summary_v2(
      spatial_data        = spatial_data,
      resolution          = resolution,
      trim                = trim,
      num_cores           = num_cores,
      control_samples     = control_samples,
      experiment_samples  = experiment_samples,
      data_type           = data_type,
      metrics             = summary_metrics,
      min_number_spots    = min_number_spots,  # <-- PRZEKAZANIE
      verbose             = verbose
    )
    
    result_summary_statistics[[data_type_name]][[paste0("resolution_", resolution)]] <- summary_data
    
    # ---- TESTS ----
    for (metric in statistic_metrics) {
      if (verbose) cat("\n>>> Testing metric:", metric, "\n")
      
      summary_data <- perform_statistical_tests_v2(
        spatial_data           = spatial_data,
        summary_data           = summary_data,
        metric                 = metric,
        resolution             = resolution,
        num_cores              = num_cores,
        mean_threshold         = mean_threshold,
        control_samples        = control_samples,
        experiment_samples     = experiment_samples,
        quantile_normalization = quantile_normalization,
        verbose                = verbose
      )
      
      error_mask <- sapply(summary_data, function(x) inherits(x, "cluster_error"))
      error_names <- names(summary_data)[error_mask]
      summary_data <- summary_data[!error_mask]
      
      if (verbose && length(error_names) > 0) {
        cat(">>> Removed", length(error_names), "cluster(s) due to errors in metric:", metric, "\n")
        cat(">>> Removed clusters:", paste(error_names, collapse = ", "), "\n")
      }
    }
    
    result_summary_statistics[[data_type_name]][[paste0("resolution_", resolution)]] <- summary_data
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(result_summary_statistics)
}
