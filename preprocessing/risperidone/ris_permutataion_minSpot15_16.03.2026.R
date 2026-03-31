message("1. Starting permutation analysis...")

library(dplyr)
library(tidyr)
library(magrittr)
library(parallel)
library(preprocessCore)
# library(moments)
library(tibble)
library(Matrix)

load("results/risperidone/risperidone-half.RData")


message("2. Loaded data ...")


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

permutation_spatial_fdr <-
  function(spatial_data,
           control_treatment,
           experiment_treatment,
           number_permutation = 10, 
           seed = NULL,
           log2ratio_threshold = 0.8,
           resolution = 0.4,
           mean_threshold = 0.2,
           test_type = "wilcoxon_test",
           num_cores = 24,
           data_type = "raw_data",
           min_number_spot = 20,
           quantile_normalization = TRUE) {
    
    # Setting the seed for reproducibility
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    # Accessing sample information
    sample_information <- spatial_data$sample_information
    
    # Summarizing treatment data
    treatment_info <- sample_information$treatment %>% 
      sort() %>% 
      table() %>% 
      as.data.frame() %>% 
      set_colnames(c("treatment", "freq"))
    
    n_control <- treatment_info %>% 
      filter(treatment == control_treatment) %>% 
      .$freq
    
    # Calculating frequencies
    n_experiment <- treatment_info %>% 
      filter(treatment == experiment_treatment) %>% 
      .$freq
    
    # Selecting a random sample
    # Generating random sample vectors
    random_sample_vectors <- lapply(1:number_permutation, function(x) sample(spatial_data$samples))
    
    significant_permutation_list <- list()
    
    for (number in 1:length(random_sample_vectors)) {
      print(number)
      
      control_samples <- random_sample_vectors[[number]][1:n_control]
      experiment_samples <- random_sample_vectors[[number]][(n_control + 1):length(random_sample_vectors[[number]])]
      
      compute_data_summary(
        spatial_data = spatial_data,
        resolution = resolution,
        trim = 0.05,
        num_cores = num_cores,
        control_samples = control_samples,
        experiment_samples = experiment_samples,
        data_type = data_type,
        min_number_spot = min_number_spot,
        metrics = c("mean")
      ) -> summary_data
      
      perform_statistical_tests(
        spatial_data = spatial_data,
        summary_data = summary_data, 
        metric = "mean", 
        resolution = resolution,
        num_cores = num_cores,
        mean_threshold = mean_threshold,
        control_samples = control_samples,
        experiment_samples = experiment_samples,
        quantile_normalization = quantile_normalization
      ) -> summary_data
      
      lapply(summary_data, function(data) {
        data$statistics$mean %>% 
          as.data.frame() %>% 
          filter(.data[[test_type]] < 0.01) %>% 
          filter(log2ratio > log2ratio_threshold) %>% 
          filter(control_mean > mean_threshold | experiment_mean > mean_threshold) %>% 
          rownames() %>%  length()
      }) %>% unlist() -> significant_permutation_vector
      
      print(significant_permutation_vector)
      significant_permutation_list[[number]] <- significant_permutation_vector
    }
    
    return(significant_permutation_list)
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

message("3. Loaded functions ...")

# trzeba zmienić parametry wewnątrz funkcji, albo dodać parametry do funkcji: 
# zmiana testu na t.test, log2ratio, control_mean, experimenta_mean, resolution, 

# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 0.5,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 1000,
#                         min_number_spot = 15,
#                         seed = 7) -> permuation_fdr_risperidone_res0.5
# 
# # print(permuation_fdr_risperidone_res0.5

# 
# 
# save(permuation_fdr_risperidone_res0.5, file = "results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.5_risperidone.RData")

# load("results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.5_risperidone.RData")
# 
# number_origin_signif_results <- c(
#   "cluster_0" = 6,
#   "cluster_1" = 1,
#   "cluster_2" = 3,
#   "cluster_3" = 5,
#   "cluster_4" = 2,
#   "cluster_5" = 16,
#   "cluster_6" = 25,
#   "cluster_7" = 14,
#   "cluster_8" = 16,
#   "cluster_9" = 6,
#   "cluster_10" = 9,
#   "cluster_11" = 89,
#   "cluster_12" = 42,
#   "cluster_13" = 17,
#   "cluster_14" = 18,
#   "cluster_15" = 0,
#   "cluster_16" = 0
# ) 
# 
# number_origin_signif_results_df <- tibble(
#   cluster = names(number_origin_signif_results),
#   signif_origin = as.numeric(number_origin_signif_results)
# )
# 
# permuation_fdr_risperidone_res0.5 %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "n_permutation") %>%
#   gather(value = "signif_random", key = "cluster", -n_permutation) %>%
#   # left_join(., number_origin_signif_results, by = "cluster") %>%
#   left_join(number_origin_signif_results_df, by = "cluster") %>% 
#   mutate(asses_random_origin = signif_random >= signif_origin) %>%   # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
#   # filter(cluster == "cluster_2")
#   group_by(cluster) %>%
#   summarise(sum_true = sum(asses_random_origin))  %>%
#   as.data.frame() %>% 
#   mutate(permutation_fdr = sum_true/1000) 




# 
# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 0.8,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 2,
#                         min_number_spot = 15,
#                         seed = 7) -> permuation_fdr_risperidone_res0.8
# 
# 
# save(permuation_fdr_risperidone_res0.8, file = "results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.8_risperidone.RData")
# 

# message("Analysis for log2ratio = 1")
# 
# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 1,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 1000,
#                         min_number_spot = 15,
#                         seed = 7) -> permuation_fdr_risperidone_res1
# 
# 
# save(permuation_fdr_risperidone_res1, file = "results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio1_risperidone.RData")

# 
# message("Analysis for log2ratio = 0.8")
# 
# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 0.8,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 1000,
#                         min_number_spot = 15,
#                         seed = 7) -> permuation_fdr_risperidone_res1
# 
# 
# save(permuation_fdr_risperidone_res1, file = "results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.8_risperidone.RData")
# 
# message("4. Done ...")
# 
# 
# number_origin_signif_results <- c(
#   "cluster_0" = 1,
#   "cluster_1" = 0,
#   "cluster_2" = 1,
#   "cluster_3" = 0,
#   "cluster_4" = 1,
#   "cluster_5" = 2,
#   "cluster_6" = 9,
#   "cluster_7" = 2,
#   "cluster_8" = 2,
#   "cluster_9" = 0,
#   "cluster_10" = 0,
#   "cluster_11" = 8,
#   "cluster_12" = 7,
#   "cluster_13" = 3,
#   "cluster_14" = 3,
#   "cluster_15" = 12,
#   "cluster_16" = 10
# )
# 
# number_origin_signif_results_df <- tibble(
#   cluster = names(number_origin_signif_results),
#   signif_origin = as.numeric(number_origin_signif_results)
# )
# 
# load("results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio1_risperidone.RData")
# 
# permuation_fdr_risperidone_res1 %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "n_permutation") %>%
#   gather(value = "signif_random", key = "cluster", -n_permutation) %>%
#   # left_join(., number_origin_signif_results, by = "cluster") %>%
#   left_join(number_origin_signif_results_df, by = "cluster") %>%
#   mutate(asses_random_origin = signif_random >= signif_origin) %>%   # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
#   # filter(cluster == "cluster_2")
#   group_by(cluster) %>%
#   summarise(sum_true = sum(asses_random_origin))  %>%
#   as.data.frame() %>%
#   mutate(permutation_fdr = sum_true/1000)


# permuation_fdr_risperidone_res0.5 %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "n_permutation") %>%
#   gather(value = "signif_random", key = "cluster", -n_permutation) %>%
#   left_join(., number_origin_signif_results, by = "cluster") %>%
#   mutate(asses_random_origin = signif_random >= signif_origin) %>%  # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
#   # filter(cluster == "cluster_2")
#   group_by(cluster) %>%
#   summarise(sum_true = sum(asses_random_origin))  %>%
#   as.data.frame() %>%
#   mutate(permutation_fdr = sum_true/1000) %>%
#   select(c(cluster, permutation_fdr)) %>%
#   write.table(file = "results/risperidone/permutation_results_16.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.5_risperidone.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


################################################################################
# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 0.5,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 1000,
#                         min_number_spot = 20,
#                         seed = 7) -> permuation_fdr_risperidone_res0.5
# # 
# # # print(permuation_fdr_risperidone_res0.5
# 
# save(permuation_fdr_risperidone_res0.5, file = "results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.5_risperidone.RData")
# 
# load("results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.5_risperidone.RData")
# 
# number_origin_signif_results <- c(
#   "cluster_0" = 6,
#   "cluster_1" = 1,
#   "cluster_2" = 3,
#   "cluster_3" = 5,
#   "cluster_4" = 2,
#   "cluster_5" = 16,
#   "cluster_6" = 25,
#   "cluster_7" = 14,
#   "cluster_8" = 16,
#   "cluster_9" = 6,
#   "cluster_10" = 9,
#   "cluster_11" = 89,
#   "cluster_12" = 42,
#   "cluster_13" = 17,
#   "cluster_14" = 17,
#   "cluster_15" = 0,
#   "cluster_16" = 0
# )
# 
# number_origin_signif_results_df <- tibble(
#   cluster = names(number_origin_signif_results),
#   signif_origin = as.numeric(number_origin_signif_results)
# )
# 
# permuation_fdr_risperidone_res0.5 %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "n_permutation") %>%
#   gather(value = "signif_random", key = "cluster", -n_permutation) %>%
#   # left_join(., number_origin_signif_results, by = "cluster") %>%
#   left_join(number_origin_signif_results_df, by = "cluster") %>%
#   mutate(asses_random_origin = signif_random >= signif_origin) %>%   # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
#   # filter(cluster == "cluster_2")
#   group_by(cluster) %>%
#   summarise(sum_true = sum(asses_random_origin))  %>%
#   as.data.frame() %>%
#   mutate(permutation_fdr = sum_true/1000)
# 




permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
                        control_treatment = "saline",
                        experiment_treatment = "risperidone",
                        data_type = "quantile_normalize_resolution_0.4",
                        mean_threshold = 0.2,
                        log2ratio_threshold = 0.8,
                        resolution = 0.4,
                        test_type = "wilcoxon_test",
                        num_cores = 10,
                        quantile_normalization = FALSE,
                        number_permutation = 1000,
                        min_number_spot = 20,
                        seed = 7) -> permuation_fdr_risperidone_res0.8


save(permuation_fdr_risperidone_res0.8, file = "results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.8_risperidone.RData")
# 
# 
load("results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.8_risperidone.RData")

number_origin_signif_results <- c(
  "cluster_0" = 2,
  "cluster_1" = 0,
  "cluster_2" = 2,
  "cluster_3" = 0,
  "cluster_4" = 1,
  "cluster_5" = 4,
  "cluster_6" = 13,
  "cluster_7" = 5,
  "cluster_8" = 6,
  "cluster_9" = 0,
  "cluster_10" = 4,
  "cluster_11" = 25,
  "cluster_12" = 22,
  "cluster_13" = 9,
  "cluster_14" = 6,
  "cluster_15" = 0,
  "cluster_16" = 0
)

number_origin_signif_results_df <- tibble(
  cluster = names(number_origin_signif_results),
  signif_origin = as.numeric(number_origin_signif_results)
)

permuation_fdr_risperidone_res0.8 %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  rownames_to_column(var = "n_permutation") %>%
  gather(value = "signif_random", key = "cluster", -n_permutation) %>%
  # left_join(., number_origin_signif_results, by = "cluster") %>%
  left_join(number_origin_signif_results_df, by = "cluster") %>%
  mutate(asses_random_origin = signif_random >= signif_origin) %>%   # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
  # filter(cluster == "cluster_2")
  group_by(cluster) %>%
  summarise(sum_true = sum(asses_random_origin))  %>%
  as.data.frame() %>%
  mutate(permutation_fdr = sum_true/1000)




# message("Analysis for log2ratio = 1")
# 
# permutation_spatial_fdr(spatial_data = risperidone_st_data_half,
#                         control_treatment = "saline",
#                         experiment_treatment = "risperidone",
#                         data_type = "quantile_normalize_resolution_0.4",
#                         mean_threshold = 0.2,
#                         log2ratio_threshold = 1,
#                         resolution = 0.4,
#                         test_type = "wilcoxon_test",
#                         num_cores = 10,
#                         quantile_normalization = FALSE,
#                         number_permutation = 1000,
#                         min_number_spot = 20,
#                         seed = 7) -> permuation_fdr_risperidone_res1
# 
# 
# save(permuation_fdr_risperidone_res1, file = "results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio0.8_risperidone.RData")
# 
# message("4. Done ...")


# number_origin_signif_results <- c(
#   "cluster_0" = 1,
#   "cluster_1" = 0,
#   "cluster_2" = 1,
#   "cluster_3" = 0,
#   "cluster_4" = 1,
#   "cluster_5" = 2,
#   "cluster_6" = 9,
#   "cluster_7" = 2,
#   "cluster_8" = 2,
#   "cluster_9" = 0,
#   "cluster_10" = 0,
#   "cluster_11" = 8,
#   "cluster_12" = 7,
#   "cluster_13" = 3,
#   "cluster_14" = 3,
#   "cluster_15" = 0,
#   "cluster_16" = 0
# )
# 
# number_origin_signif_results_df <- tibble(
#   cluster = names(number_origin_signif_results),
#   signif_origin = as.numeric(number_origin_signif_results)
# )
# load("results/risperidone/permutation_results_17.03.2026/permutation_results_res0.4_mean0.2p0.01log2ratio1_risperidone.RData")
# 
# permuation_fdr_risperidone_res1 %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "n_permutation") %>%
#   gather(value = "signif_random", key = "cluster", -n_permutation) %>%
#   # left_join(., number_origin_signif_results, by = "cluster") %>%
#   left_join(number_origin_signif_results_df, by = "cluster") %>%
#   mutate(asses_random_origin = signif_random >= signif_origin) %>%   # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
#   # filter(cluster == "cluster_2")
#   group_by(cluster) %>%
#   summarise(sum_true = sum(asses_random_origin))  %>%
#   as.data.frame() %>%
#   mutate(permutation_fdr = sum_true/1000)

