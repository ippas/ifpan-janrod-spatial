compute_data_summary_v5 <- function(spatial_data,
                                    resolution = 0.8,
                                    trim = 0.05,
                                    num_cores = 24,
                                    control_samples,
                                    experiment_samples,
                                    data_type = "raw_data",
                                    min_number_spots = 20,
                                    metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis")) {
  
  library(foreach)
  library(doParallel)
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort()
  spatial_data_local <- spatial_data  # lokalna kopia
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  message(">> Starting compute_data_summary_v4 using foreach and ", num_cores, " workers.")
  start_time <- Sys.time()
  
  results <- foreach(cluster = unique_clusters,
                     .export = c("skewness", "kurtosis"),
                     .packages = c("dplyr")) %dopar% {
                       message(">> Processing cluster ", cluster)
                       
                       tryCatch({
                         all_barcodes <- spatial_data_local[[data_type]]$metadata %>%
                           right_join(spatial_data_local$clusters, ., by = c("barcode", "sample")) %>%
                           filter((!!sym(resolution_column)) == cluster) %>%
                           mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
                           mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
                           select(group, sample_barcode)
                         
                         barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
                         control_barcodes <- barcodes_list[["control"]]
                         experiment_barcodes <- barcodes_list[["experiment"]]
                         
                         control_expression_spot <- sapply(control_samples, function(sample_id) {
                           spatial_data_local[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)], drop = FALSE]
                         })
                         
                         control_list <- list()
                         if ("expression_spot" %in% metrics) control_list$expression_spot <- control_expression_spot
                         if ("sum" %in% metrics) control_list$sum <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, sum(y))))
                         if ("mean" %in% metrics) control_list$mean <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, mean(y, trim = trim, na.rm = TRUE))))
                         if ("median" %in% metrics) control_list$median <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, median(y))))
                         if ("IQR" %in% metrics) control_list$IQR <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, IQR(y))))
                         if ("diff_range" %in% metrics) control_list$diff_range <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, diff(range(y)))))
                         if ("var" %in% metrics) control_list$var <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, var(y))))
                         if ("skewness" %in% metrics) control_list$skewness <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, skewness(y))))
                         if ("kurtosis" %in% metrics) control_list$kurtosis <- sapply(control_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, kurtosis(y))))
                         
                         experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
                           spatial_data_local[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)], drop = FALSE]
                         })
                         
                         experiment_list <- list()
                         if ("expression_spot" %in% metrics) experiment_list$expression_spot <- experiment_expression_spot
                         if ("sum" %in% metrics) experiment_list$sum <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, sum(y))))
                         if ("mean" %in% metrics) experiment_list$mean <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, mean(y, trim = trim, na.rm = TRUE))))
                         if ("median" %in% metrics) experiment_list$median <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, median(y))))
                         if ("IQR" %in% metrics) experiment_list$IQR <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, IQR(y))))
                         if ("diff_range" %in% metrics) experiment_list$diff_range <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, diff(range(y)))))
                         if ("var" %in% metrics) experiment_list$var <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, var(y))))
                         if ("skewness" %in% metrics) experiment_list$skewness <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, skewness(y))))
                         if ("kurtosis" %in% metrics) experiment_list$kurtosis <- sapply(experiment_expression_spot, function(x) apply(x, 1, function(y) ifelse(length(y) < min_number_spots, NA, kurtosis(y))))
                         
                         return(list(
                           peak = spatial_data_local[[data_type]]$annotate$peak_id,
                           gene = spatial_data_local[[data_type]]$annotate$gene_name,
                           control = control_list,
                           experiment = experiment_list
                         ))
                       }, error = function(e) {
                         message("!! ERROR in cluster ", cluster, ": ", e$message)
                         return(NULL)
                       })
                     }
  
  stopCluster(cl)
  end_time <- Sys.time()
  message(">> Total computation time: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  
  names(results) <- paste0("cluster_", unique_clusters)
  
  failed_clusters <- which(sapply(results, is.null))
  if (length(failed_clusters) > 0) {
    message("!! The following clusters failed: ", paste(unique_clusters[failed_clusters], collapse = ", "))
  } else {
    message(">> All clusters computed successfully.")
  }
  
  return(results)
}


perform_statistical_tests_v5 <- function(spatial_data,
                                         summary_data,
                                         metric = "mean",
                                         resolution = 0.8,
                                         num_cores = 8,
                                         mean_threshold = 0,
                                         control_samples,
                                         experiment_samples,
                                         quantile_normalization = FALSE,
                                         verbose = TRUE) {
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort()
  cluster_ids <- paste0("cluster_", unique_clusters)
  
  if (verbose) {
    message(">> Starting v5 foreach-parallel statistical tests for ", length(cluster_ids), " clusters...")
  }
  
  # Set up parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results <- foreach(cluster = cluster_ids,
                     .packages = c("preprocessCore"),
                     .export = c("summary_data", "spatial_data", "control_samples", "experiment_samples"),
                     .combine = "list",
                     .errorhandling = "pass") %dopar% {
                       
                       output <- tryCatch({
                         message("\n--- START CLUSTER: ", cluster)
                         
                         cluster_data <- summary_data[[cluster]]
                         if (is.null(cluster_data)) {
                           message("!!! cluster_data is NULL")
                           return(NULL)
                         }
                         
                         control_expression <- tryCatch(cluster_data$control[[metric]][, control_samples], error = function(e) {
                           message("!!! Failed to extract control_expression: ", e$message)
                           return(NULL)
                         })
                         
                         experiment_expression <- tryCatch(cluster_data$experiment[[metric]][, experiment_samples], error = function(e) {
                           message("!!! Failed to extract experiment_expression: ", e$message)
                           return(NULL)
                         })
                         
                         if (is.null(control_expression) || is.null(experiment_expression)) {
                           message("!!! Missing expression matrix — skipping ", cluster)
                           return(NULL)
                         }
                         
                         if (quantile_normalization) {
                           message(">>> Performing quantile normalization...")
                           expression_raw <- cbind(control_expression, experiment_expression)
                           expression_norm <- preprocessCore::normalize.quantiles(expression_raw)
                           
                           colnames(expression_norm) <- colnames(expression_raw)
                           rownames(expression_norm) <- rownames(expression_raw)
                           control_expression <- expression_norm[, control_samples]
                           experiment_expression <- expression_norm[, experiment_samples]
                           rm(expression_raw)
                         }
                         
                         n <- nrow(control_expression)
                         t_test <- numeric(n)
                         wilcoxon_test <- numeric(n)
                         ks_test <- numeric(n)
                         condition <- character(n)
                         
                         for (i in seq_len(n)) {
                           if (all(is.na(control_expression[i, ])) | all(is.na(experiment_expression[i, ]))) {
                             condition[i] <- "all_NA"
                             t_test[i] <- 1
                             wilcoxon_test[i] <- 1
                             ks_test[i] <- 1
                           } else if (mean(control_expression[i, ], na.rm = TRUE) < mean_threshold |
                                      mean(experiment_expression[i, ], na.rm = TRUE) < mean_threshold) {
                             condition[i] <- "low_mean"
                             t_test[i] <- 1
                             wilcoxon_test[i] <- 1
                             ks_test[i] <- 1
                           } else {
                             condition[i] <- "yes_mean"
                             t_test[i] <- tryCatch(t.test(control_expression[i, ], experiment_expression[i, ], var.equal = TRUE)$p.value,
                                                   error = function(e) { message("t-test error [", i, "]: ", e$message); return(1) })
                             wilcoxon_test[i] <- tryCatch(wilcox.test(control_expression[i, ], experiment_expression[i, ])$p.value,
                                                          error = function(e) { message("wilcoxon error [", i, "]: ", e$message); return(1) })
                             ks_test[i] <- tryCatch(ks.test(control_expression[i, ], experiment_expression[i, ])$p.value,
                                                    error = function(e) { message("ks-test error [", i, "]: ", e$message); return(1) })
                           }
                         }
                         
                         control_mean <- rowMeans(control_expression, na.rm = TRUE)
                         experiment_mean <- rowMeans(experiment_expression, na.rm = TRUE)
                         control_sd <- apply(control_expression, 1, sd, na.rm = TRUE)
                         experiment_sd <- apply(experiment_expression, 1, sd, na.rm = TRUE)
                         control_se <- control_sd / sqrt(ncol(control_expression))
                         experiment_se <- experiment_sd / sqrt(ncol(experiment_expression))
                         log2ratio <- log2(experiment_mean / control_mean)
                         
                         stats <- list(
                           t_test = t_test,
                           wilcoxon_test = wilcoxon_test,
                           ks_test = ks_test,
                           log2ratio = log2ratio,
                           condition = condition,
                           control_mean = control_mean,
                           experiment_mean = experiment_mean,
                           control_sd = control_sd,
                           experiment_sd = experiment_sd,
                           control_se = control_se,
                           experiment_se = experiment_se
                         )
                         
                         if (!is.null(cluster_data$control$skewness) && !is.null(cluster_data$experiment$kurtosis)) {
                           message(">>> Adding skewness/kurtosis info...")
                           stats$control_skewness <- apply(cluster_data$control$skewness, 1, mean, na.rm = TRUE)
                           stats$experiment_skewness <- apply(cluster_data$experiment$skewness, 1, mean, na.rm = TRUE)
                           stats$control_kurtosis <- apply(cluster_data$control$kurtosis, 1, mean, na.rm = TRUE)
                           stats$experiment_kurtosis <- apply(cluster_data$experiment$kurtosis, 1, mean, na.rm = TRUE)
                         }
                         
                         cluster_data$statistics[[metric]] <- stats
                         message("--- END CLUSTER: ", cluster)
                         return(cluster_data)
                       }, error = function(e) {
                         message("!!! Cluster ", cluster, " failed: ", e$message)
                         return(NULL)
                       })
                       
                       output
                     }
  
  stopCluster(cl)
  names(results) <- cluster_ids
  results <- results[!sapply(results, is.null)]
  return(results)
}


summarize_and_test_v5 <- function(spatial_data,
                                  trim = 0.05,
                                  num_cores = 24,
                                  data_params_df,
                                  control_samples,
                                  experiment_samples,
                                  mean_threshold = 0,
                                  metrics = c("mean", "median"),              # do podsumowań
                                  statistic_metrics = c("mean", "median")) {  # do testów
  
  result_summary_statistics <- list()
  
  for (i in 1:nrow(data_params_df)) {
    data_type <- data_params_df[i, 1]
    resolution <- data_params_df[i, 2]
    data_type_name <- data_params_df[i, 3]
    quantile_normalization <- data_params_df[i, 4]
    
    message(">> Running compute_data_summary_v3 for: ", data_type_name, " @ resolution = ", resolution)
    
    summary_data <- compute_data_summary_v5(
      spatial_data = spatial_data,
      resolution = resolution,
      trim = trim,
      num_cores = num_cores,
      control_samples = control_samples,
      experiment_samples = experiment_samples,
      data_type = data_type,
      metrics = metrics,
      min_number_spots = 20
    )
    
    for (metric in statistic_metrics) {
      message(">> Performing statistical tests for metric: ", metric)
      
      # summary_data <- perform_statistical_tests_v5(
      #   spatial_data = spatial_data,
      #   summary_data = summary_data,
      #   metric = metric,
      #   resolution = resolution,
      #   num_cores = num_cores,
      #   mean_threshold = mean_threshold,
      #   control_samples = control_samples,
      #   experiment_samples = experiment_samples,
      #   quantile_normalization = quantile_normalization
      # )
    }
    
    name_resolution <- paste0("resolution_", resolution)
    result_summary_statistics[[data_type_name]][[name_resolution]] <- summary_data
  }
  
  return(result_summary_statistics)
}

result_summary_statistics_v5 <- summarize_and_test_v5(
  spatial_data = ris3q29_st_data,
  trim = 0.05,
  num_cores = 20,
  data_params_df = data_params_df[3,],
  control_samples = samples_wt,
  experiment_samples = samples_del,
  mean_threshold = 0,
  metrics = c("sum", "mean", "median"),
  statistic_metrics = c("mean", "median")
)

result_summary_statistics_v5$raw_data$resolution_0.4$cluster_
