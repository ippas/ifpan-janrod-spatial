# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")


summarize_and_test <- function(spatial_data,
                               trim = 0.05,
                               num_cores = 24,
                               data_params_df,
                               control_samples,
                               experiment_samples,
                               mean_threshold = 0,
                               min_number_spots = 20,                    # <-- NOWE
                               metrics = c("mean", "median"),            # <-- do testów
                               statistic_metrics = c("mean", "median"))  # <-- do compute
{
  result_summary_statistics <- list()
  
  for(i in 1:nrow(data_params_df)){
    data_type <- data_params_df[i, 1]
    resolution <- data_params_df[i, 2]
    data_type_name <- data_params_df[i, 3]
    quantile_normalization <- data_params_df[i, 4]
    
    print(data_type)
    print(resolution)
    
    summary_data <- compute_data_summary(
      spatial_data = spatial_data,
      resolution = resolution,
      trim = trim,
      num_cores = num_cores,
      control_samples = control_samples,
      experiment_samples = experiment_samples,
      data_type = data_type,
      min_number_spots = min_number_spots,          # <-- PRZEKAZANIE
      metrics = metrics
    )
    
    for(metric in statistic_metrics){
      summary_data <- perform_statistical_tests(
        spatial_data = spatial_data,
        summary_data = summary_data,
        metric = metric,
        resolution = resolution,
        num_cores = num_cores,
        mean_threshold = mean_threshold,
        control_samples = control_samples,
        experiment_samples = experiment_samples,
        quantile_normalization = quantile_normalization
      )
    }
    
    name_resolution <- paste0("resolution_", resolution)
    result_summary_statistics[[data_type_name]][[name_resolution]] <- summary_data
  }
  
  return(result_summary_statistics)
}


compute_data_summary <- function(spatial_data,
                                 resolution = 0.8,
                                 trim = 0.05,
                                 num_cores = 24,
                                 control_samples,
                                 experiment_samples,
                                 data_type = "raw_data",
                                 min_number_spots = 20,
                                 metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis")) {   # Add metrics parameter
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <-
    unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  compute_data_summary_cluster <- function(cluster) {
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
      select(group, sample_barcode)
    
    # Split barcodes into control and experiment groups based on group label
    barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
    
    # Extract control and experiment barcodes from the list
    control_barcodes <- barcodes_list[["control"]]
    experiment_barcodes <- barcodes_list[["experiment"]]
    
    control_expression_spot <- sapply(control_samples, function(sample_id) {
      # Selecting the control data based on the sample ID
      control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
    })
    
    # Calculate only the metrics that are included in the metrics parameter for control group
    control_list <- list()
    if("expression_spot" %in% metrics) control_list$expression_spot <- control_expression_spot
    if("sum" %in% metrics) control_list$sum <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, sum(y))})})
    if("mean" %in% metrics) control_list$mean <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, mean(y, trim = trim, na.rm = TRUE))})})
    if("median" %in% metrics) control_list$median <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, median(y))})})
    if("IQR" %in% metrics) control_list$IQR <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, IQR(y))})})
    if("diff_range" %in% metrics) control_list$diff_range <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, diff(range(y)))})})
    if("var" %in% metrics) control_list$var <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, var(y))})})
    if("skewness" %in% metrics) control_list$skewness <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, skewness(y))})})
    if("kurtosis" %in% metrics) control_list$kurtosis <- sapply(control_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, kurtosis(y))})})
    
    
    # ... Repeat for other metrics ...
    experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
      # Selecting the experiment data based on the sample ID
      experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
    })
    
    
    # Calculate only the metrics that are included in the metrics parameter for experiment group
    experiment_list <- list()
    if("expression_spot" %in% metrics) experiment_list$expression_spot <- experiment_expression_spot
    if("sum" %in% metrics) experiment_list$sum <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, sum(y))})})
    if("mean" %in% metrics) experiment_list$mean <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, mean(y, trim = trim, na.rm = TRUE))})})
    if("median" %in% metrics) experiment_list$median <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, median(y))})})
    if("IQR" %in% metrics) experiment_list$IQR <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, IQR(y))})})
    if("diff_range" %in% metrics) experiment_list$diff_range <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, diff(range(y)))})})
    if("var" %in% metrics) experiment_list$var <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, var(y))})})
    if("skewness" %in% metrics) experiment_list$skewness <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, skewness(y))})})
    if("kurtosis" %in% metrics) experiment_list$kurtosis <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(y){ifelse(length(y) < min_number_spots, NA, kurtosis(y))})})
    
    
    return(list(peak = spatial_data[[data_type]]$annotate$peak_id,
                gene = spatial_data[[data_type]]$annotate$gene_name,
                control = control_list,
                experiment = experiment_list
    ))
  }
  
  start_time <- Sys.time()
  # Use mclapply instead of lapply
  results <- mclapply(unique_clusters, compute_data_summary_cluster, mc.cores = num_cores)
  # Set names for the results
  results <- setNames(results, paste0("cluster_", unique_clusters))
  end_time <-  Sys.time()
  diff_time <- end_time - start_time
  print(diff_time)
  
  is_try_error <- sapply(results, function(x) class(x)[1] == "try-error")
  results <- results[!is_try_error]
  
  return(results)
}


perform_statistical_tests <- function(spatial_data,
                                      summary_data,
                                      metric = "mean",
                                      resolution = 0.8,
                                      num_cores = 24,
                                      mean_threshold = 0,
                                      control_samples,
                                      experiment_samples,
                                      quantile_normalization = FALSE) {
  
  
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <-
    unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  
  perform_statistical_tests_cluster <- function(summary_data, cluster, metric = "mean", mean_threshold = 0){
    cluster_data <- summary_data[[cluster]]
    
    control_expression <- cluster_data$control[[metric]][, control_samples]
    experiment_expression <- cluster_data$experiment[[metric]][, experiment_samples]
    
    # Check if both control_expression and experiment_expression are completely NA
    if (all(is.na(control_expression)) && all(is.na(experiment_expression))) {
      # return a list of NA vectors for each element
      results <- list(
        t_test = rep(NA, nrow(control_expression)),
        wilcoxon_test = rep(NA, nrow(control_expression)),
        ks_test = rep(NA, nrow(control_expression)),
        log2ratio = rep(NA, nrow(control_expression)),
        condition = rep(NA, nrow(control_expression)),
        control_mean = rep(NA, nrow(control_expression)),
        experiment_mean = rep(NA, nrow(control_expression))
      )
      
      return(results)
      
    } 
    
    if(quantile_normalization == TRUE) {
      expression_raw <- cbind(control_expression, experiment_expression)
      
      expression_norm <- normalize.quantiles(expression_raw) 
      
      colnames(expression_norm) <- colnames(expression_raw)
      rownames(expression_norm) <- rownames(expression_raw)
      
      control_expression <-  expression_norm[, control_samples]
      experiment_expression <-  expression_norm[, experiment_samples]
      
      rm(expression_raw)
      
    }
    
    t_test <- numeric(nrow(control_expression))
    wilcoxon_test <- numeric(nrow(control_expression))
    ks_test <- numeric(nrow(control_expression))
    condition <- character(nrow(control_expression))
    
    # Iterate over each gene
    for (i in seq_len(length(cluster_data$peak))) {
      # Check if all values in the row are NA
      if (all(is.na(control_expression[i, ])) | all(is.na(experiment_expression[i, ]))) {
        condition[i] <- "all_NA"
        t_test[i] <- 1
        wilcoxon_test[i] <- 1
        ks_test[i] <- 1
      }
      # Check if the mean of either control or experiment group is below the threshold
      else if (mean(control_expression[i, ], na.rm = TRUE) < mean_threshold | mean(experiment_expression[i, ], na.rm = TRUE) < mean_threshold) {
        condition[i] <- "low_mean"
        t_test[i] <- 1
        wilcoxon_test[i] <- 1
        ks_test[i] <- 1
      } 
      else {
        condition[i] <- "yes_mean"
        
        # Perform a t-test if data are not constant
        t_test_result <- tryCatch({
          t.test(control_expression[i, ], experiment_expression[i, ], var.equal = TRUE)
        }, error = function(e) {
          # In case of any error, return a list with p-value as 1
          return(list(p.value = 1))
        })
        t_test[i] <- t_test_result$p.value
        
        # Perform a Wilcoxon test if mean is above the threshold
        wilcox_result <- tryCatch({
          wilcox.test(control_expression[i, ], experiment_expression[i, ])
        }, error = function(e) {
          # In case of any error, return a list with p-value as 1
          return(list(p.value = 1))
        })
        wilcoxon_test[i] <- wilcox_result$p.value
        
        # Perform a KS test if mean is above the threshold
        ks_result <- tryCatch({
          ks.test(control_expression[i, ], experiment_expression[i, ])
        }, error = function(e) {
          # In case of any error, return a list with p-value as 1
          return(list(p.value = 1))
        })
        ks_test[i] <- ks_result$p.value
      }
    }
    
    control_mean <- apply(control_expression, 1, mean, na.rm = T)
    experiment_mean <- apply(experiment_expression, 1, mean, na.rm = T)
    
    control_sd <- apply(control_expression, 1, sd, na.rm = T)
    experiment_sd <- apply(experiment_expression, 1, sd, na.rm = T)
    
    # # Calculate the Standard Error (SE) for each row in control_expression and experiment_expression
    control_se <- apply(control_expression, 1, function(row) {
      sd(row, na.rm = TRUE) / sqrt(sum(!is.na(row)))
    })
    
    experiment_se <- apply(experiment_expression, 1, function(row) {
      sd(row, na.rm = TRUE) / sqrt(sum(!is.na(row)))
    })
    
    
    log2ratio <- log2(rowMeans(experiment_expression, na.rm = TRUE) / rowMeans(control_expression, na.rm = TRUE)) %>%
      unname()
    
    if (!is.null(cluster_data$experiment$kurtosis) &
        !is.null(cluster_data$experiment$skewness)) {
      control_skewness <-  apply( cluster_data$control$skewness, 1, mean, na.rm = T)
      experiment_skewness <-  apply( cluster_data$experiment$skewness, 1, mean, na.rm = T)
      control_kurtosis <-  apply( cluster_data$control$kurtosis, 1, mean, na.rm = T)
      experiment_kurtosis <-  apply( cluster_data$experiment$kurtosis, 1, mean, na.rm = T)
      
      results_statatistics <- list(
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
        experiment_se = experiment_se,
        control_skewness = control_skewness,
        experiment_skewness = experiment_skewness,
        control_kurtosis = control_kurtosis,
        experiment_kurtosis = experiment_kurtosis
      )
    } else {
      results_statatistics <- list(
        t_test = t_test,
        wilcoxon_test = wilcoxon_test,
        ks_test = ks_test,
        log2ratio = log2ratio,
        condition = condition,
        control_mean = control_mean,
        experiment_mean = experiment_mean
      )
    }
    
    cluster_data$statistics[[metric]] <- results_statatistics
    
    return(cluster_data)
  }
  
  paste0("cluster_", unique_clusters) -> unique_clusters
  
  start_time <- Sys.time()
  # Use mclapply instead of lapply
  results <-
    mclapply(
      unique_clusters,
      perform_statistical_tests_cluster,
      summary_data = summary_data,
      mc.cores = num_cores,
      metric = metric,
      mean_threshold = mean_threshold
    )
  # Set names for the results
  results <- setNames(results, unique_clusters)
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  print(diff_time)
  
  is_try_error <- sapply(results, function(x) class(x)[1] == "try-error")
  results <- results[!is_try_error]
  
  return(results)
}


# executes seurat analysis for risperidone
path_to_data <- "data/risperidone/spaceranger_defaultReference_half_18.03.2026/"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
# sample_names <- list.files(path = path_to_data)

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))

# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]

sample_names <- c(samples_saline, samples_risperidone)

# # read peaks for risperidone
# info_peaks_risperidone <- read_gene_annotation(file_path = "data/risperidone/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for risperidone
risperidone_integrate_half_defaultReference <- integrate_data_seurat(path_to_data, sample_names, nfeatures, dims)

# read images to spatial transcriptoms for risperidone
images_risperidone_half_defaultReference <- create_images_tibble(path_to_data, sample_names)

# read barcode data for risperidone
barcode_risperidone_half_defaultReference <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_risperidone_half_defaultReference
)



risperidone_integrate_half_defaultReference@assays$RNA %>% .[1:10, 1:10]

data.frame(
  gene_name = rownames(risperidone_integrate_half_defaultReference@assays$RNA),
  peak_id   =   rownames(risperidone_integrate_half_defaultReference@assays$RNA)
) -> info_peaks_risperidone_defaultReference
  

# poprawić do 
risperidone_st_data_half_defaultReference <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_risperidone,
                                                barcode_info = barcode_risperidone_half_defaultReference,
                                                images_info = images_risperidone_half_defaultReference,
                                                integrated_data = risperidone_integrate_half_defaultReference,
                                                peaks_info = info_peaks_risperidone_defaultReference)

risperidone_st_data_half_defaultReference <- add_seurat_data(spatial_data = risperidone_st_data_half_defaultReference,
                                            integrated_data = risperidone_integrate_half_defaultReference)

risperidone_st_data_half_defaultReference <- add_filtered_data(spatial_data = risperidone_st_data_half_defaultReference,
                                              mean_expression_threshold = 0.5)

risperidone_st_data_half_defaultReference <- add_colfilt_data(spatial_data = risperidone_st_data_half_defaultReference,
                                             min_spot_threshold = 0,
                                             expression_threshold = 2)

risperidone_st_data_half_defaultReference <- add_range_normalize_data(spatial_data = risperidone_st_data_half_defaultReference,
                                                     range = 1500,
                                                     flatten = 1,
                                                     threshold = 500)

risperidone_st_data_half_defaultReference <- add_clusters_data(spatial_data = risperidone_st_data_half_defaultReference,
                                              integrated_data = risperidone_integrate_half_defaultReference,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

# risperidone_st_data_half <- evaluate_clustering_stability(spatial_data = risperidone_st_data_half,
#                                                           seurat_object = risperidone_integrate_half,
#                                                           resolution_start = 0.05,
#                                                           resolution_end = 2,
#                                                           resolution_step = 0.05)




# for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)){
#   risperidone_st_data_half_defaultReference <-
#     add_quantile_norm_data(
#       spatial_data = risperidone_st_data_half_defaultReference,
#       resolution = {{resolution}},
#       num_cores = 10,
#       data_type = "raw_data"
#     )
# }

save(risperidone_st_data_half_defaultReference,
     sample_names,
     samples_saline,
     samples_risperidone,
     risperidone_integrate_half_defaultReference,
     file = "results/risperidone/ris_DefaultRef_checkpointIntegrationData_17.03.2026.RData")


# load("results/risperidone/ris_DefaultRef_checkpointIntegrationData_17.03.2026.RData")

add_quantile_norm_data <- function(
  spatial_data,
  resolution = 0.8,
  data_type = "raw_data",
  verbose = TRUE
) {
  # This function normalizes gene expression data based on quantiles
  # within each cluster for spatial transcriptomics data.
  # The normalization is performed on a given resolution level.
  #
  # Args:
  #   spatial_data: A list containing spatial data
  #   resolution: The resolution level used for clusters
  #   data_type: The type of data to normalize
  #   verbose: Whether to print progress information
  #
  # Returns:
  #   Updated spatial_data object with quantile-normalized data
  
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' is required.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("Package 'rlang' is required.")
  }
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  name_data <- paste0("quantile_normalize_resolution_", resolution)
  
  if (!resolution_column %in% colnames(spatial_data$clusters)) {
    stop("Column '", resolution_column, "' not found in spatial_data$clusters.")
  }
  
  merged_metadata <- spatial_data[[data_type]]$metadata %>%
    dplyr::right_join(spatial_data$clusters, by = c("barcode", "sample")) %>%
    dplyr::mutate(sample_barcode = paste(sample, barcode, sep = "_"))
  
  unique_clusters <- sort(unique(merged_metadata[[resolution_column]]))
  n_clusters <- length(unique_clusters)
  
  if (verbose) {
    message("Starting quantile normalization")
    message("Resolution: ", resolution)
    message("Data type: ", data_type)
    message("Number of clusters: ", n_clusters)
  }
  
  results <- vector("list", length = n_clusters)
  names(results) <- paste0("cluster_", unique_clusters)
  
  total_start_time <- Sys.time()
  
  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    cluster_start_time <- Sys.time()
    
    if (verbose) {
      message("\n[", i, "/", n_clusters, "] Processing cluster: ", cluster)
    }
    
    all_barcodes <- merged_metadata %>%
      dplyr::filter((!!rlang::sym(resolution_column)) == cluster) %>%
      dplyr::pull(sample_barcode)
    
    if (length(all_barcodes) == 0) {
      warning("Cluster ", cluster, " has no barcodes. Skipping.")
      next
    }
    
    cluster_data <- spatial_data[[data_type]]$data[, all_barcodes, drop = FALSE] %>%
      as.matrix()
    
    if (verbose) {
      message("  Number of barcodes: ", length(all_barcodes))
      message("  Matrix dimensions: ", nrow(cluster_data), " x ", ncol(cluster_data))
    }
    
    norm_data <- preprocessCore::normalize.quantiles(cluster_data)
    
    colnames(norm_data) <- colnames(cluster_data)
    rownames(norm_data) <- rownames(cluster_data)
    
    results[[i]] <- norm_data
    
    rm(cluster_data, norm_data)
    invisible(gc())
    
    cluster_end_time <- Sys.time()
    cluster_time <- as.numeric(difftime(cluster_end_time, cluster_start_time, units = "secs"))
    total_time_so_far <- as.numeric(difftime(cluster_end_time, total_start_time, units = "mins"))
    
    if (verbose) {
      message("  Cluster done in: ", round(cluster_time, 2), " s")
      message("  Total elapsed: ", round(total_time_so_far, 2), " min")
    }
  }
  
  if (verbose) {
    message("\nCombining normalized clusters...")
  }
  
  quantile_normalize_data <- do.call(cbind, results) %>% Matrix::Matrix()
  
  spatial_data[[name_data]] <- list()
  spatial_data[[name_data]]$metadata <- spatial_data[[data_type]]$metadata
  spatial_data[[name_data]]$annotate <- spatial_data[[data_type]]$annotate
  spatial_data[[name_data]]$data <- quantile_normalize_data
  
  total_end_time <- Sys.time()
  total_time <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))
  
  if (verbose) {
    message("All done.")
    message("Total execution time: ", round(total_time, 2), " min")
  }
  
  return(spatial_data)
}

risperidone_st_data_half_defaultReference <-
  add_quantile_norm_data(
    spatial_data = risperidone_st_data_half_defaultReference,
    resolution = 0.4,
    data_type = "raw_data",
    # ram_limit_gb = 100,
    verbose = TRUE
  )

risperidone_st_data_half_defaultReference <-
  add_quantile_norm_data(
    spatial_data = risperidone_st_data_half_defaultReference,
    resolution = 0.3,
    data_type = "raw_data",
    # ram_limit_gb = 100,
    verbose = TRUE
  )

save(risperidone_st_data_half_defaultReference,
     sample_names,
     samples_saline,
     samples_risperidone,
     risperidone_integrate_half_defaultReference,
     file = "results/risperidone/ris_DefaultRef_checkpointSecondIntegrationData_17.03.2026.RData")
# 
# load("results/risperidone/ris_DefaultRef_checkpointSecondIntegrationData_17.03.2026.RData")



spatial_cluster(spatial_data = risperidone_st_data_half_defaultReference,
                resolution = 0.3,
                samples = c(samples_saline, samples_risperidone),
                palette = palette_allen,
                size= 1.0,
                ncol = 4)




data.frame(
  data_type = c(
    "raw_data", "raw_data",
    "quantile_normalize_resolution_0.3",
    "quantile_normalize_resolution_0.4"
  ),
  resolution = c(0.3, 0.4, 0.3, 0.4),
  data_type_name = c(
    "raw_data", "raw_data",
    "quantile_normalize",
    "quantile_normalize"
  )
) %>%
  mutate(
    quantile_normalization = F
  ) -> data_params_df

summarize_and_test(spatial_data = risperidone_st_data_half_defaultReference,
                   trim = 0.05, 
                   num_cores = 10,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   min_number_spots = 20,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> risperidone_summary_statistics_half_defaultReference



save(samples_saline,
     samples_risperidone,
     risperidone_integrate_half_defaultReference,
     risperidone_st_data_half_defaultReference,
     risperidone_summary_statistics_half_defaultReference,
     file = "results/risperidone/risperidone_defaultRef_17.03.2026.RData")


load("results/risperidone/risperidone_defaultRef_17.03.2026.RData")
load("results/risperidone/risperidone-half.RData")
# lapply(
#   names(risperidone_summary_statistics_half_defaultReference$quantile_normalize$resolution_0.4),
#   function(cluster_name) {
#     risperidone_summary_statistics_half_defaultReference$quantile_normalize$resolution_0.4[[cluster_name]]$statistics$mean %>%
#       as.data.frame() %>%
#       dplyr::filter(wilcoxon_test < 0.01) %>%
#       dplyr::filter(abs(log2ratio) > 0.8) %>%
#       dplyr::filter(experiment_mean > 0.2 | control_mean > 0.2)
#   }
# ) %>% unlist %>% unique %>% cat(sep = "\n")

risperidone_st_data_half_defaultReference$clusters %>% 
  select(sample, barcode, cluster_resolution_0.5) %>% 
  group_by(sample, cluster_resolution_0.5) %>% nest %>% 
  mutate(number_spot = map(data, ~nrow(.x))) %>% unnest(number_spot) %>% 
  filter(number_spot < 20)


# ##############################################################################
# supplementary analysi for the smae clusters as first risperidone.Rdata
# ##############################################################################
risperidone_st_data_half_defaultReference$clusters %>% head

risperidone_st_data_half$clusters %>% dim


risperidone_st_data_half_defaultReferenceClustersCustomRef <- risperidone_st_data_half_defaultReference

risperidone_st_data_half_defaultReferenceClustersCustomRef$clusters <- risperidone_st_data_half$clusters 

# spatial_cluster(spatial_data = risperidone_st_data_half_defaultReferenceClustersCustomRef,
#                 resolution = 0.4,
#                 samples = c(samples_saline, samples_risperidone),
#                 palette = palette_allen,
#                 size= 1.0,
#                 ncol = 4)
# spatial_interest_cluster(cluster = 6,
#                          # seurat_object = integrated_analysis,
#                          spatial_data = risperidone_st_data_half_defaultReferenceClustersCustomRef,
#                          resolution = 0.4,
#                          samples = c(samples_saline, samples_risperidone),
#                          size= 1,
#                          ncol = 4)


risperidone_st_data_half_defaultReferenceClustersCustomRef <-
  add_quantile_norm_data(
    spatial_data = risperidone_st_data_half_defaultReferenceClustersCustomRef,
    resolution = 0.4,
    data_type = "raw_data",
    verbose = TRUE
  )

risperidone_st_data_half_defaultReferenceClustersCustomRef <-
  add_quantile_norm_data(
    spatial_data = risperidone_st_data_half_defaultReferenceClustersCustomRef,
    resolution = 0.3,
    data_type = "raw_data",
    verbose = TRUE
  )

summarize_and_test(spatial_data = risperidone_st_data_half_defaultReferenceClustersCustomRef,
                   trim = 0.05, 
                   num_cores = 10,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   min_number_spots = 20,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> risperidone_summary_statistics_half_defaultReferenceClustersCustomRef


save(samples_saline,
     samples_risperidone,
     risperidone_integrate_half_defaultReference,
     risperidone_st_data_half_defaultReferenceClustersCustomRef,
     risperidone_summary_statistics_half_defaultReferenceClustersCustomRef,
     file = "results/risperidone/risperidone_defaultReferenceClustersCustomRef_18.03.2026.RData")


risperidone_summary_statistics_half_defaultReferenceClustersCustomRef$quantile_normalize$resolution_0.4$cluster_6$statistics$mean %>% 
  as.data.frame() %>% 
  filter(wilcoxon_test < 0.01) %>% 
  filter(log2ratio > 0.8) %>% 
  filter(control_mean > 0.2 | experiment_mean > 0.2)



lapply(
  names(risperidone_summary_statistics_half_defaultReferenceClustersCustomRef$quantile_normalize$resolution_0.4),
  function(cluster_name) {
    risperidone_summary_statistics_half_defaultReferenceClustersCustomRef$quantile_normalize$resolution_0.4[[cluster_name]]$statistics$mean %>%
      as.data.frame() %>%
      dplyr::filter(wilcoxon_test < 0.01) %>%
      dplyr::filter(abs(log2ratio) > 0.8) %>%
      dplyr::filter(experiment_mean > 0.2 | control_mean > 0.2)
  } 
)


lapply(
  names(risperidone_summary_statistics_half_defaultReferenceClustersCustomRef$quantile_normalize$resolution_0.4),
  function(cluster_name) {
    risperidone_summary_statistics_half_defaultReferenceClustersCustomRef$quantile_normalize$resolution_0.4[[cluster_name]]$statistics$mean %>%
      as.data.frame() %>% .["Gm14308", , drop = FALSE]
  } 
)


