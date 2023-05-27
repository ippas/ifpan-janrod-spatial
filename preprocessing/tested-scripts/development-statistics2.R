spatial_data <- risperidone_st_data_half
data_type <- "quantile_normalize"
trim <- 0.05

investigate_summary_cluster <- function(data, metric = "mean", gene){
  data.frame(gene = data$gene,
             peak = data$peak) %>% 
    filter(gene == {{gene}}) %>% 
    pull(peak) -> peak
  print(peak)
  
  cat("control\n")
  print(data$control[[metric]][peak, ])
  cat("experiment\n")
  print(data$experiment[[metric]][peak, ])
  t.test(data$control[[metric]][peak, ], data$experiment[[metric]][peak, ]) -> t_test
  print(t_test)
}

get_results_df <- function(summary_data, metric = "mean", cluster = "cluster_0") {
  statistics <- summary_data[[cluster]]$statistics[[metric]] %>% as.data.frame()
  
  peak <- summary_data[[cluster]]$peak
  gene <- summary_data[[cluster]]$gene
  
  cbind(peak, gene, statistics)
}


compute_data_summary <-
  function(spatial_data,
           resolution = 0.8,
           trim = 0.05,
           num_cores = 24,
           control_samples,
           experiment_samples,
           data_type = "quantile_normalize",
           metrics = c("sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis")) {   # Add metrics parameter
    
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
      control_list$expression_spot <- control_expression_spot  
      if("sum" %in% metrics) control_list$sum <- sapply(control_expression_spot, function(x){apply(x, 1, sum)})
      if("mean" %in% metrics) control_list$mean <- sapply(control_expression_spot, function(x){apply(x, 1, mean, trim = trim, na.rm = TRUE)})
      if("median" %in% metrics) control_list$median <- sapply(control_expression_spot, function(x){apply(x, 1, median)})
      if("IQR" %in% metrics) control_list$IQR <- sapply(control_expression_spot, function(x){apply(x, 1, IQR)})
      if("diff_range" %in% metrics) control_list$diff_range <- sapply(control_expression_spot, function(x){apply(x, 1, function(x){range(x) %>% diff()})})
      if("var" %in% metrics) control_list$var <- sapply(control_expression_spot, function(x){apply(x, 1, var)})
      if("skewness" %in% metrics) control_list$skewness <- sapply(control_expression_spot, function(x){apply(x, 1, skewness)})
      if("kurtosis" %in% metrics) control_list$kurtosis <- sapply(control_expression_spot, function(x){apply(x, 1, kurtosis)})
      
      
      # ... Repeat for other metrics ...
      experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
        # Selecting the experiment data based on the sample ID
        experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
      })
      
      
      # Calculate only the metrics that are included in the metrics parameter for experiment group
      experiment_list <- list()
      experiment_list$expression_spot <- experiment_expression_spot
      if("sum" %in% metrics) experiment_list$sum <- sapply(experiment_expression_spot, function(x){apply(x, 1, sum)})
      if("mean" %in% metrics) experiment_list$mean <- sapply(experiment_expression_spot, function(x){apply(x, 1, mean, trim = trim, na.rm = TRUE)})
      if("median" %in% metrics) experiment_list$median <- sapply(experiment_expression_spot, function(x){apply(x, 1, median)})
      if("IQR" %in% metrics) experiment_list$IQR <- sapply(experiment_expression_spot, function(x){apply(x, 1, IQR)})
      if("diff_range" %in% metrics) experiment_list$diff_range <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(x){range(x) %>% diff()})})
      if("var" %in% metrics) experiment_list$var <- sapply(experiment_expression_spot, function(x){apply(x, 1, var)})
      if("skewness" %in% metrics) experiment_list$skewness <- sapply(experiment_expression_spot, function(x){apply(x, 1, skewness)})
      if("kurtosis" %in% metrics) experiment_list$kurtosis <- sapply(experiment_expression_spot, function(x){apply(x, 1, kurtosis)})
      
      
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
    
    return(results)
  }


perform_statistical_tests <- function(spatial_data, summary_data, metric = "mean", resolution = 0.8, num_cores = 1, mean_threshold = 0,
                                      control_samples, experiment_samples) {
  
  
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <-
    unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  
  perform_statistical_tests_cluster <- function(summary_data, cluster, metric = "mean", mean_threshold = 0){
    cluster_data <- summary_data[[cluster]]
    
    control_expression <- cluster_data$control[[metric]][, control_samples]
    experiment_expression <- cluster_data$experiment[[metric]][, experiment_samples]
    
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
    
    log2ratio <- log2(rowMeans(experiment_expression, na.rm = TRUE) / rowMeans(control_expression, na.rm = TRUE)) %>%
      unname()
    
    results_statatistics <- list(
      t_test = t_test,
      wilcoxon_test = wilcoxon_test,
      ks_test = ks_test,
      log2ratio = log2ratio,
      condition = condition
    )
    
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
  
  # summary_data$statistics[[metric]] <- results
  
  return(results)
}





