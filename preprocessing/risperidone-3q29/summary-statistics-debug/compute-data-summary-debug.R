compute_data_summary_v2 <- function(spatial_data,
                                 resolution = 0.8,
                                 trim = 0.05,
                                 num_cores = 24,
                                 control_samples,
                                 experiment_samples,
                                 data_type = "raw_data",
                                 min_number_spots = 20,
                                 metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis"),
                                 verbose = FALSE) {
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort()
  
  if (verbose) {
    cat(">>> START compute_data_summary()\n")
    cat(">>> Resolution column:", resolution_column, "\n")
    cat(">>> Unique clusters:", paste(unique_clusters, collapse = ", "), "\n")
    cat(">>> Metrics:", paste(metrics, collapse = ", "), "\n\n")
  }
  
  # Wewnętrzna funkcja dla pojedynczego klastra
  compute_data_summary_cluster <- function(cluster) {
    if (verbose) cat(">>> Processing cluster:", cluster, "\n")
    
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
      select(group, sample_barcode)
    
    barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
    control_barcodes <- barcodes_list[["control"]]
    experiment_barcodes <- barcodes_list[["experiment"]]
    
    control_expression_spot <- sapply(control_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, grep(paste0(sample_id, "_"), control_barcodes), drop = FALSE]
    }, simplify = FALSE)
    
    experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, grep(paste0(sample_id, "_"), experiment_barcodes), drop = FALSE]
    }, simplify = FALSE)
    
    control_list <- list()
    experiment_list <- list()
    
    # Obliczenia dla zadanych metryk
    for (metric in metrics) {
      if (metric == "expression_spot") {
        control_list$expression_spot <- control_expression_spot
        experiment_list$expression_spot <- experiment_expression_spot
        next
      }
      
      control_list[[metric]] <- sapply(control_expression_spot, function(x) {
        apply(x, 1, function(y) {
          if (length(y) < min_number_spots) NA else switch(metric,
                                                           sum       = sum(y, na.rm = TRUE),
                                                           mean      = mean(y, trim = trim, na.rm = TRUE),
                                                           median    = median(y, na.rm = TRUE),
                                                           IQR       = IQR(y, na.rm = TRUE),
                                                           diff_range= diff(range(y, na.rm = TRUE)),
                                                           var       = var(y, na.rm = TRUE),
                                                           skewness  = skewness(y),
                                                           kurtosis  = kurtosis(y),
                                                           NA
          )
        })
      })
      
      experiment_list[[metric]] <- sapply(experiment_expression_spot, function(x) {
        apply(x, 1, function(y) {
          if (length(y) < min_number_spots) NA else switch(metric,
                                                           sum       = sum(y, na.rm = TRUE),
                                                           mean      = mean(y, trim = trim, na.rm = TRUE),
                                                           median    = median(y, na.rm = TRUE),
                                                           IQR       = IQR(y, na.rm = TRUE),
                                                           diff_range= diff(range(y, na.rm = TRUE)),
                                                           var       = var(y, na.rm = TRUE),
                                                           skewness  = skewness(y),
                                                           kurtosis  = kurtosis(y),
                                                           NA
          )
        })
      })
    }
    
    return(list(
      peak = spatial_data[[data_type]]$annotate$peak_id,
      gene = spatial_data[[data_type]]$annotate$gene_name,
      control = control_list,
      experiment = experiment_list
    ))
  }
  
  start_time <- Sys.time()
  
  results <- mclapply(unique_clusters, function(cl) {
    if (verbose) cat("\n>>> STARTING CLUSTER:", cl, "\n")
    tryCatch({
      compute_data_summary_cluster(cl)
    }, error = function(e) {
      if (verbose) message("!! ERROR in cluster ", cl, ": ", e$message)
      return(structure(list(error_message = e$message), class = "cluster_error"))
    })
  }, mc.cores = num_cores)
  
  names(results) <- paste0("cluster_", unique_clusters)
  
  # Zgłoś, które klastry miały błąd
  error_clusters <- sapply(results, function(x) inherits(x, "cluster_error"))
  if (any(error_clusters)) {
    message("\n!!! UWAGA: błędy w klastrach: ", paste(names(results)[error_clusters], collapse = ", "))
  }
  
  end_time <- Sys.time()
  if (verbose) cat("\n>>> Czas obliczeń:", round(end_time - start_time, 2), "sekund\n")
  
  return(results)
}


res <- compute_data_summary_v2(
  spatial_data        = ris3q29_st_data,
  resolution          = data_params_df[29, 2],
  trim                = 0.05,
  num_cores           = 20,
  control_samples     = samples_wt,
  experiment_samples  = samples_del,
  data_type           = data_params_df[5, 1],
  metrics             = c("mean", "median", "sum"),
  verbose             = TRUE
)

test_res <- perform_statistical_tests_v2(
  spatial_data         = ris3q29_st_data,
  summary_data         = res,
  metric               = c("mean"),
  resolution           = data_params_df[29, 2],
  num_cores            = 20,
  mean_threshold       = 0,
  control_samples      = samples_wt,
  experiment_samples   = samples_del,
  quantile_normalization = FALSE,
  verbose              = TRUE
)


result <- summarize_and_test_v2(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df[c(1,5), , drop = FALSE],
  control_samples      = samples_wt,
  experiment_samples   = samples_del,
  summary_metrics      = c("mean", "median", "sum"),
  statistic_metrics    = c("mean", "median"),
  mean_threshold       = 0,
  num_cores            = 20,
  verbose              = TRUE
)


result$range_normalize$resolution_0.2$cluster_15


# test_res$cluster_0$statistics$mean$



