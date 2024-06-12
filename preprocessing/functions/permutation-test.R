function(spatial_data,
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


permutation_test <-
  function(spatial_data,
           control_treatment,
           experiment_treatment,
           number_permutation = 10, 
           seed = NULL) {
  
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
  # Generating 10 vectors of random samples
  random_sample_vectors <- lapply(1:number_permutation, function(x) sample(spatial_data$samples))
  
  # print(random_sample_vectors)
  
  random_sample_vectors[[1]][1:n_control] -> control_samples
  random_sample_vectors[[1]][(n_control+1):length(random_sample_vectors[[1]])] -> experiment_samples
  
  # print(random_sample_vectors)
  
  significant_permutation_list <- list()
  
  for(number in 1:length(random_sample_vectors)){
    print(number)
    
    random_sample_vectors[[number]][1:n_control] -> control_samples
    random_sample_vectors[[number]][(n_control+1):length(random_sample_vectors[[number]])] -> experiment_samples
    
    # print(control_samples)
    
    compute_data_summary(spatial_data = spatial_data,
                         resolution = 0.4,
                         trim = 0.05,
                         num_cores = 24,
                         control_samples = control_samples,
                         experiment_samples = experiment_samples,
                         data_type = "raw_data",
                         metrics = c("mean")) -> summary_data
    
    perform_statistical_tests(spatial_data = spatial_data,
                              summary_data = summary_data, 
                              metric = "mean", 
                              resolution = 0.4,
                              num_cores = 24,
                              mean_threshold = 0.2,
                              control_samples = control_samples,
                              experiment_samples = experiment_samples,
                              quantile_normalization = T) -> summary_data
    
    lapply(summary_data, function(data){
      data$statistics$mean %>% 
        as.data.frame() %>% 
        filter(wilcoxon_test < 0.01) %>% 
        filter(log2ratio > 0.8) %>% 
        filter(control_mean > 0.2 | experiment_mean > 0.2) %>% 
        rownames() %>%  length()
    }) %>% unlist() -> significant_permutation_vector
    
    print(significant_permutation_vector)
    significant_permutation_list[[number]] <- significant_permutation_vector
  }
  
  return(significant_permutation_list)
  
  # # Return a list with all the results
  # list(sample_information = sample_information,
  #      random_sample_vectors = random_sample_vectors,
  #      summary_data = summary_data)
}

permutation_test(spatial_data = risperidone_st_data_half,
                control_treatment = "saline",
                experiment_treatment = "risperidone",
                number_permutation = 10,
                seed = 10) -> tmp

number_signif_results_per_cluster <- c("cluster_0"=2, "cluster_2"=2, "cluster_4"=1, "cluster_5"=4, "cluster_6"=13, "cluster_7"=5, "cluster_8"=6, "cluster_10"=4, "cluster_11"=25, "cluster_12"=22, "cluster_13"=9, "cluster_14"=6)

permutation_test(spatial_data = risperidone_st_data_half,
                 control_treatment = "saline",
                 experiment_treatment = "risperidone",
                 number_permutation = 1000,
                 seed = 7) -> permutation_results_list_0.01

# save(permutation_results_list, "results/tmp.RData")
save(permutation_results_list_0.01, file = "results/risperidone/permutation_resutls_0.01.RData")
