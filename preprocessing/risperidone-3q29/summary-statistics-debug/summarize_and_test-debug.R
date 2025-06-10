summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05,
                   num_cores = 20,  # debug = lepiej 1 core
                   data_params_df = data_params_df[29, , drop = FALSE],
                   control_samples = samples_wt,
                   experiment_samples = samples_del,
                   mean_threshold = 0,
                   statistic_metrics = "mean",
                   metrics = "mean") -> wtDel_debug




summary_data <- compute_data_summary(
  spatial_data        = ris3q29_st_data,
  resolution          = 0.4,
  # trim                = 0.05,
  num_cores           = 15,
  control_samples     = samples_wt,
  experiment_samples  = samples_del,
  data_type           = "raw_data",
  min_number_spots    = 20,
  metrics             = c("mean", "median")
)


summary_data <- compute_data_summary(
  spatial_data        = ris3q29_st_data,
  resolution          = data_params_df[5, 2],
  trim                = 0.05,
  num_cores           = 20,
  control_samples     = samples_wt,
  experiment_samples  = samples_del,
  data_type           = data_params_df[5, 1],
  metrics             = "mean"
)



compute_data_summary_cluster <- function(cluster,
                                         spatial_data,
                                         resolution,
                                         control_samples,
                                         experiment_samples,
                                         data_type = "raw_data",
                                         trim = 0.05,
                                         min_number_spots = 20,
                                         metrics = c("mean"),
                                         verbose = FALSE) {
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
    select(group, sample_barcode)
  
  barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
  control_barcodes <- barcodes_list[["control"]]
  experiment_barcodes <- barcodes_list[["experiment"]]
  
  if (verbose) {
    cat("==== DEBUG: CLUSTER", cluster, "====\n")
    cat("Control barcodes:", length(control_barcodes), "\n")
    cat("Experiment barcodes:", length(experiment_barcodes), "\n\n")
  }
  
  control_expression_spot <- sapply(control_samples, function(sample_id) {
    matched <- grep(paste0(sample_id, "_"), control_barcodes, value = TRUE)
    if (verbose) cat("  Control sample:", sample_id, " ->", length(matched), "spots\n")
    spatial_data[[data_type]]$data[, matched, drop = FALSE]
  }, simplify = FALSE)
  
  experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
    matched <- grep(paste0(sample_id, "_"), experiment_barcodes, value = TRUE)
    if (verbose) cat("  Experiment sample:", sample_id, " ->", length(matched), "spots\n")
    spatial_data[[data_type]]$data[, matched, drop = FALSE]
  }, simplify = FALSE)
  
  control_list <- list()
  experiment_list <- list()
  
  if ("mean" %in% metrics) {
    control_list$mean <- sapply(control_expression_spot, function(x) {
      apply(x, 1, function(y) {
        if (length(y) < min_number_spots) NA else mean(y, trim = trim, na.rm = TRUE)
      })
    })
    
    experiment_list$mean <- sapply(experiment_expression_spot, function(x) {
      apply(x, 1, function(y) {
        if (length(y) < min_number_spots) NA else mean(y, trim = trim, na.rm = TRUE)
      })
    })
    
    if (verbose) {
      cat("\n--- Summary of `control_list$mean` (first 2 peaks):\n")
      print(head(control_list$mean, 2))
      cat("\n--- Summary of `experiment_list$mean` (first 2 peaks):\n")
      print(head(experiment_list$mean, 2))
    }
  }
  
  return(list(
    peak = spatial_data[[data_type]]$annotate$peak_id,
    gene = spatial_data[[data_type]]$annotate$gene_name,
    control = control_list,
    experiment = experiment_list
  ))
}


compute_data_summary_cluster(
  cluster             = 1,
  spatial_data        = ris3q29_st_data,
  resolution          = data_params_df[5, 2],
  control_samples     = samples_wt,
  experiment_samples  = samples_del,
  data_type           = data_params_df[5, 1],
  trim                = 0.05,
  metrics             = "mean",
  verbose             = TRUE
) -> res

res


results_list <- list()

for (cl in 0:20) {
  cat("\n==============================\n")
  cat(">>> START KLASTRA:", cl, "\n")
  cat("==============================\n")
  
  res <- compute_data_summary_cluster(
    cluster             = cl,
    spatial_data        = ris3q29_st_data,
    resolution          = data_params_df[5, 2],
    control_samples     = samples_wt,
    experiment_samples  = samples_del,
    data_type           = data_params_df[5, 1],
    trim                = 0.05,
    metrics             = "mean",
    verbose             = TRUE
  )
  
  results_list[[paste0("cluster_", cl)]] <- res
  
  cat("\n>>>>>>>>>>> KLASTER", cl, "ZAKOŃCZONY PRAWIDŁOWO <<<<<<<<<<<\n")
}
resolution_column <- paste0("cluster_resolution_0.4")
unique_clusters <- unique(ris3q29_st_data$clusters[[resolution_column]]) %>% sort()

