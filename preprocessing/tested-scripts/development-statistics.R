spatial_data <- risperidone_st_data_half
data_type = "raw_data"
resolution_column <- "cluster_resolution_0.8"
num_cores <- 24

unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort


control_samples <- samples_saline
experiment_samples <- samples_risperidone
# Define a function to process each cluster




###############################################################################3
# 7. version -> add condition with minimum 
###############################################################################3

process_cluster <- function(cluster, expression_unit, trim, min_spots = 30, mean_threshold = 0.8, quartile_normalization = FALSE) {
  # Extract barcodes for both control and experiment groups
  # Join metadata and cluster information, and filter by cluster
  # A new column 'sample_barcode' is created by concatenating 'sample' and 'barcode'
  # Group label ('control' or 'experiment') is assigned based on sample
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
  
  # Prepare data for the statistical test depending on expression_unit
  # Preparing data for statistical tests
  if (expression_unit == "spot") {
    # Direct retrieval of expression data when unit is "spot"
    control_expression <- spatial_data[[data_type]]$data[, control_barcodes]
    experiment_expression <- spatial_data[[data_type]]$data[, experiment_barcodes]
  } else if (expression_unit == "sample") {
    # Calculating mean expression values when unit is "sample"
    control_expression <- sapply(control_samples, function(sample_id) {
      # Selecting the control data based on the sample ID
      control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
      if (ncol(control_data) >= min_spots) { 
        # Calculate the mean when the number of columns is greater than or equal to the minimum spots 
        apply(control_data, 1, mean, trim = trim) 
      } else { 
        # Return a vector of NAs with names from the control data rows when the number of columns is less than minimum spots
        setNames(rep(NA, nrow(control_data)), rownames(control_data)) 
      }
    })
    
    experiment_expression <- sapply(experiment_samples, function(sample_id) {
      # Selecting the experiment data based on the sample ID
      experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
      if (ncol(experiment_data) >= min_spots) { 
        # Calculate the mean when the number of columns is greater than or equal to the minimum spots 
        apply(experiment_data, 1, mean, trim = trim) 
      } else { 
        # Return a vector of NAs with names from the experiment data rows when the number of columns is less than minimum spots
        setNames(rep(NA, nrow(experiment_data)), rownames(experiment_data)) 
      }
    })
  }
  
  control_raw_samples <- sapply(control_samples, function(sample_id) {
    # Selecting the control data based on the sample ID
    control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
  })
  
    experiment_raw_samples <- sapply(experiment_samples, function(sample_id) {
      # Selecting the experiment data based on the sample ID
      experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
    })
    
    control_median <- sapply(control_raw_samples, function(x){apply(x, 1, median)})
    control_IQR <- sapply(control_raw_samples, function(x){apply(x, 1, IQR)})
    control_range <- sapply(control_raw_samples, function(x){apply(x, 1, range)})
    control_var <- sapply(control_raw_samples, function(x){apply(x, 1, var)})
    control_skewness <- sapply(control_raw_samples, function(x){apply(x, 1, skewness)})
    control_kurtosis <- sapply(control_raw_samples, function(x){apply(x, 1, kurtosis)})
    
    control_metrics <- list(median = control_median,
                            IQR = control_IQR,
                            range = control_range,
                            var = control_var,
                            skewness = control_skewness,
                            kurtosis = control_kurtosis)
    
    experiment_median <- sapply(experiment_raw_samples, function(x){apply(x, 1, median)})
    experiment_IQR <- sapply(experiment_raw_samples, function(x){apply(x, 1, IQR)})
    experiment_range <- sapply(experiment_raw_samples, function(x){apply(x, 1, range)})
    experiment_var <- sapply(experiment_raw_samples, function(x){apply(x, 1, var)})
    experiment_skewness <- sapply(experiment_raw_samples, function(x){apply(x, 1, skewness)})
    experiment_kurtosis <- sapply(experiment_raw_samples, function(x){apply(x, 1, kurtosis)})
    
    experiment_metrics <- list(median = experiment_median,
                               IQR = experiment_IQR,
                               range = experiment_range,
                               var = experiment_var,
                               skewness = experiment_skewness,
                               kurtosis = experiment_kurtosis)
  
  if(quartile_normalization == TRUE){
    cbind(control_expression, experiment_expression) -> data_matrix
    rownames_original <- rownames(data_matrix)
    colnames_original <- colnames(data_matrix)
    
    # Perform quantile normalization
    data_matrix_normalized <- normalize.quantiles(data_matrix)
    
    # Restore the rownames and colnames
    rownames(data_matrix_normalized) <- rownames_original
    colnames(data_matrix_normalized) <- colnames_original
    
    control_expression <- data_matrix_normalized[, control_samples]
    experiment_expression <- data_matrix_normalized[, experiment_samples]
    
  }
  
  
  # Initialize a vector to store p-values for each gene
  t_test <- numeric(nrow(control_expression))
  wilcoxon_test <- numeric(nrow(control_expression))
  ks_test <- numeric(nrow(control_expression))
  condition <- character(nrow(control_expression))
  
  # Iterate over each gene
  for (i in seq_len(nrow(control_expression))) {
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
      condition[i] <- "yes_var_mean"
      
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

  # Calculate the log2 ratio of the means of the two groups
  log2ratio <- log2(rowMeans(experiment_expression, na.rm = TRUE) / rowMeans(control_expression, na.rm = TRUE)) %>%
    unname()
  
  # Return a list containing control and experiment barcodes, their expressions, and the p-values
  return(list(control_barcodes = control_barcodes, 
              experiment_barcodes = experiment_barcodes,
              control_expression = control_expression,
              experiment_expression = experiment_expression,
              control_metrics = control_metrics,
              experiment_metrics = experiment_metrics,
              peaks = rownames(control_expression),
              gene = spatial_data[[data_type]]$annotate$gene_name,
              t_test = t_test,
              wilcoxon_test = wilcoxon_test,
              ks_test = ks_test,
              log2ratio = log2ratio,
              condition = condition))
}


run_process_cluster_test <-
  function(spatial_data,
           resolution = 0.8,
           expression_unit = "sample",
           trim = 0.05,
           num_cores = 1,
           mean_threshold = 0.2,
           data_type = "quantile_normalize",
           quartile_normalization = FALSE) {
    resolution_column <- paste0("cluster_resolution_", resolution)
    unique_clusters <-
      unique(spatial_data$clusters[[resolution_column]]) %>% sort
    
    
    start_time <- Sys.time()
    results <-
      mclapply(
        unique_clusters,
        process_cluster,
        expression_unit = expression_unit,
        trim = trim,
        mc.cores = num_cores,
        mean_threshold = mean_threshold,
        quartile_normalization = quartile_normalization
      )
    results <- setNames(results, paste0("cluster_", unique_clusters))
    diff_time <- end_time <- Sys.time()
    
    print(diff_time)
    
    return(results)
  }

results <- run_process_cluster_test(
  spatial_data = risperidone_st_data_half,
  resolution = 0.8,
  expression_unit = "sample",
  trim = 0.05,
  num_cores = 24,
  mean_threshold = 0.2,
  data_type = "quantile_normalize") -> results_stat

  




start_time <- Sys.time()
# Use mclapply instead of lapply
results_range <- mclapply(unique_clusters, process_cluster, expression_unit = "sample", trim = 0.05, mc.cores = num_cores, mean_threshold = 0.2, quartile_normalization = FALSE)
# Set names for the results
results_range <- setNames(results_range, paste0("cluster_", unique_clusters))
end_time <- Sys.time()
end_time - start_time



data_type = "raw_data"
data_type = "range_normalize"
data_type = "quantile_normalize"


lapply(as.character(paste0("cluster_", unique_clusters[0:19])), function(x){
  results_range[[x]] %>% .[c(5, 6, 7, 8, 9, 10, 11)] %>% 
    data.frame() %>%
    # filter(t_test != 1) %>%
    mutate(FDR_t_test = p.adjust(t_test, method = "fdr")) %>%
    mutate(FDR_wilcox = p.adjust(wilcoxon_test, method = "fdr")) %>%
    # filter(wilcoxon_test < 0.05) %>%
    # filter(ks_test < 0.05) %>%
    filter(t_test < 0.05) %>%
    filter(log2ratio > 0.5) %>%
    filter(gene %in% michkor_anova_drug)
    # filter(gene  %in% c("Txnip", "Ddit4", "Arrdc2", "Polr3e", "Sgk1"))
    # arrange(t_test) 
})

risperidone_st_data_half$

data_type <- "range_normalize"

