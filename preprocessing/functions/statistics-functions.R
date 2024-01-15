# spatial_data <- risperidone_st_data_half
# data_type <- "quantile_normalize"
# trim <- 0.05

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

#####
compute_data_summary <-
  function(spatial_data,
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

perform_statistical_tests <-
  function(spatial_data,
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
      control_se <- apply(control_expression, 1, function(row) {sd(row) / sqrt(length(row))})
      experiment_se <- apply(experiment_expression, 1, function(row) {sd(row) / sqrt(length(row))})
      
      
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

#####
perform_statistical_tests2 <-
  function(spatial_data,
           summary_data,
           metric = "mean",
           resolution = 0.8,
           num_cores = 1,
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
    
    control_mean <- apply(control_expression, 1, mean)
    experiment_mean <- apply(experiment_expression, 1, mean)
    
    log2ratio <- log2(rowMeans(experiment_expression, na.rm = TRUE) / rowMeans(control_expression, na.rm = TRUE)) %>%
      unname()
    
    results_statatistics <- list(
      t_test = t_test,
      wilcoxon_test = wilcoxon_test,
      ks_test = ks_test,
      log2ratio = log2ratio,
      condition = condition,
      control_mean = control_mean,
      experiment_mean = experiment_mean
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



summarize_and_test <-function(spatial_data,
                              trim = 0.05,
                              num_cores = 24,
                              data_params_df,
                              control_samples,
                              experiment_samples,
                              mean_threshold = 0,
                              metrics = c("mean", "median")){
  
  result_summary_statistics <- list()
  
  for(i in 1:nrow(data_params_df)){
    data_type <- data_params_df[i,1]
    resolution <- data_params_df[i,2]
    data_type_name <- data_params_df[i,3]
    quantile_normalization <- data_params_df[i, 4]
    
    print(data_type)
    print(resolution)
    
    
    compute_data_summary(spatial_data = spatial_data,
                         resolution = resolution,
                         trim = trim,
                         num_cores = num_cores,
                         control_samples = control_samples,
                         experiment_samples = experiment_samples,
                         data_type = data_type,
                         metrics = metrics) -> summary_data 
    
    for(metric in  c("mean", "median")){
      perform_statistical_tests(spatial_data = spatial_data,
                                summary_data = summary_data, 
                                metric = metric, 
                                resolution = resolution,
                                num_cores = num_cores,
                                mean_threshold = mean_threshold,
                                control_samples = control_samples,
                                experiment_samples = experiment_samples,
                                quantile_normalization = quantile_normalization) -> summary_data
    }
    name_resolution <- paste0("resolution_", resolution)
    
    result_summary_statistics[[data_type_name]][[name_resolution]] <- summary_data
    
  }
  
  return(result_summary_statistics)
}


# Function to filter cluster statistics
filter_cluster_statistics <- function(summary_data,
                                      metric = "mean",
                                      cluster = "cluster_0",
                                      control_mean_threshold = 0,
                                      experiment_mean_threshold = 0,
                                      log2ratio_threshold = 0,
                                      t_test_threshold = 1,
                                      wilcoxon_test_threshold = 1,
                                      ks_test_threshold = 1,
                                      include_samples = FALSE) {
  
  # Function to filter cluster statistics with an option to include sample data
  # This function retrieves data from a specified cluster, 
  # filters based on provided thresholds for various statistical measures,
  # and returns a data frame of the filtered statistics.
  #
  # @param summary_data A list containing the summary data.
  # @param metric A character string specifying the metric to use. Default is "mean".
  # @param cluster A character string specifying the cluster to use. Default is "cluster_0".
  # @param control_mean_threshold A numeric value specifying the minimum acceptable mean for the control group. Default is 0.
  # @param experiment_mean_threshold A numeric value specifying the minimum acceptable mean for the experiment group. Default is 0.
  # @param log2ratio_threshold A numeric value specifying the minimum absolute log2ratio. Default is 0.
  # @param t_test_threshold A numeric value specifying the maximum acceptable t-test value. Default is 1.
  # @param wilcoxon_test_threshold A numeric value specifying the maximum acceptable Wilcoxon test value. Default is 1.
  # @param ks_test_threshold A numeric value specifying the maximum acceptable KS test value. Default is 1.
  # @param include_samples A logical value indicating whether to include sample data in the output. Default is FALSE.
  # @return A data frame of the filtered cluster statistics.
  
  # Retrieve specific cluster statistics and convert to data frame
  statistics <- summary_data[[cluster]]$statistics[[metric]] %>% as.data.frame()
  
  # Retrieve peak and gene data for specific cluster
  peak <- summary_data[[cluster]]$peak
  gene <- summary_data[[cluster]]$gene
  
  # If include_samples is TRUE, retrieve control and experiment samples for specific cluster
  if (include_samples) {
    control_samples <- summary_data[[cluster]]$control[[metric]]
    experiment_samples <- summary_data[[cluster]]$experiment[[metric]] 
    cluster_statistics_df <- cbind(peak, gene, control_samples, experiment_samples, cluster = cluster, statistics)
  } else {
    cluster_statistics_df <- cbind(peak, gene, cluster = cluster, statistics)
  }
  
  # Filter the data frame based on the provided thresholds
  cluster_statistics_df %>% 
    filter(control_mean >= control_mean_threshold | 
             experiment_mean >= experiment_mean_threshold) %>%
    filter(abs(log2ratio) >= log2ratio_threshold) %>%
    filter(t_test <= t_test_threshold) %>%
    filter(wilcoxon_test <= wilcoxon_test_threshold) %>%
    filter(ks_test <= ks_test_threshold) -> filter_cluster_statistics_df
  
  # Return the filtered data frame
  return(filter_cluster_statistics_df)
}

# Function to filter data statistics
filter_data_statistics <- function(summary_data,
                                   data_type, 
                                   resolution,
                                   ...) {
  
  # Function to filter data statistics
  # This function retrieves data of a specified type and resolution from the summary data, 
  # applies the `filter_cluster_statistics` function to each cluster,
  # and returns a combined data frame of the filtered statistics for all clusters.
  #
  # @param summary_data A list containing the summary data.
  # @param data_type A character string specifying the type of data to use (e.g., "RNA", "protein").
  # @param resolution A numeric or character value specifying the resolution to use.
  # @param ... Additional arguments to pass to the `filter_cluster_statistics` function.
  # @return A data frame of the filtered data statistics for all clusters.
  
  # Construct the name of the resolution field
  resolution_name <- paste0("resolution_", resolution)
  
  # Retrieve the names of all clusters for the specified data type and resolution
  cluster_vector <- summary_data[[data_type]][[resolution_name]] %>% names()
  
  # Apply the `filter_cluster_statistics` function to each cluster,
  # and combine the results into a single data frame
  lapply(cluster_vector, function(x){
    filter_cluster_statistics(summary_data = summary_data[[data_type]][[resolution_name]],
                              cluster = x,
                              ...)
  }) %>%
    do.call(rbind, .) -> filter_data_statistics_df
  
  # Return the filtered data frame
  return(filter_data_statistics_df)
} 


# Function to count unique genes and peaks per cluster
summary_filter_statistics <- function(data) {
  
  # Function to count unique genes and peaks per cluster
  # This function counts the number of unique peaks and genes in each cluster
  # and prints these numbers along with the total number of unique peaks and genes.
  #
  # @param data A data frame that includes at least the columns: peak, gene, and cluster.
  # @return A data frame with the number of unique genes and peaks per cluster.
  
  # Count the total number of unique peaks
  number_peaks <- data$peak %>% unique() %>% length()
  
  # Count the total number of unique genes
  number_genes <- data$gene %>% unique() %>% length()
  
  # For each cluster, count the number of unique genes and peaks
  data %>% 
    select(peak, gene, cluster) %>%   # Select relevant columns
    group_by(cluster) %>%             # Group by cluster
    nest() %>%                        # Nest the data (create a list-column of data frames)
    mutate(number_genes = purrr::map(data, ~length(unique(.x$gene)))) %>%   # Count unique genes per cluster
    mutate(number_peaks = purrr::map(data, ~length(unique(.x$peak)))) %>%   # Count unique peaks per cluster
    select(-data) %>%                 # Remove the nested data column
    unnest(cols = c(number_genes, number_peaks)) %>%  # Unnest the counts
    as.data.frame() -> genes_peak_per_cluster          # Convert tibble to data frame
  
  # Print the total number of unique peaks and genes
  cat(paste0("All significant peaks: ", number_peaks), "\n")
  cat(paste0("All significant genes: ", number_genes), "\n\n")
  
  # Print the number of unique genes and peaks per cluster
  # print(genes_peak_per_cluster)
  
  # Return the data frame
  return(genes_peak_per_cluster)
}

