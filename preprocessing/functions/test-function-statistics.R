# compute_data_summary <-
#   function(spatial_data,
#            resolution = 0.8,
#            trim = 0.05,
#            num_cores = 24,
#            control_samples,
#            experiment_samples,
#            data_type = "quantile_normalize",
#            min_number_spots = 20,
#            metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis")) {   # Add metrics parameter

compute_data_summary(spatial_data = risperidone_st_data_half,
                     data_type = "raw_data",
                     control_samples = samples_saline,
                     experiment_samples = samples_risperidone,
                     metrics = c("mean", "median", "skewness", "kurtosis")) -> tmp_summary


summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df[1,],
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis"))  -> tmp

perform_statistical_tests2(
  spatial_data = risperidone_st_data_half,
  summary_data = tmp_summary,
  metric = "mean",
  control_samples = samples_saline,
  experiment_samples = samples_risperidone
) -> tmp_stat2

tmp$raw_data$resolution_0.05$cluster_0$statistics$mean$t_test

tmp_summary$cluster_0

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
    
    return(results)
  }
