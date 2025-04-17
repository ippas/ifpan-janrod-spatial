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



# trzeba zmienić parametry wewnątrz funkcji, albo dodać parametry do funkcji: 
# zmiana testu na t.test, log2ratio, control_mean, experimenta_mean, resolution, 

permutation_spatial_fdr(spatial_data = ldopa_st_data,
                        control_treatment = "saline",
                        experiment_treatment = "ldopa",
                        data_type = "quantile_normalize_resolution_0.8",
                        mean_threshold = 0.6,
                        log2ratio_threshold = 0.5,
                        resolution = 0.8,
                        test_type = "t_test",
                        num_cores = 30,
                        quantile_normalization = FALSE,
                        number_permutation = 1000,
                        seed = 7) -> permuation_fdr_ldopa


save(permuation_fdr_ldopa, file = "results/ldopa/permuation_fdr_ldopa.RData")

names(permuation_fdr_ldopa) <- paste0("permutation", seq_along(permuation_fdr_ldopa))

read.delim(file = "results/ldopa/ldopa-signif-resutls.tsv", 
           header = TRUE,
           # row.names = F,
           # quote = F,
           sep = "\t") %>% 
  select(c("peak", "gene", "cluster")) %>% 
  group_by(cluster) %>% 
  summarise(n = n()) %>% as.data.frame() %>% 
  set_colnames(c("cluster", "signif_origin"))-> number_origin_signif_results



permuation_fdr_ldopa %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "n_permutation") %>% 
  gather(value = "signif_random", key = "cluster", -n_permutation) %>% 
  left_join(., number_origin_signif_results, by = "cluster") %>% 
  mutate(asses_random_origin = signif_random >= signif_origin) %>%  # Poprawka: TRUE/FALSE zamiast znakowego "TRUE"/"FALSE"
  # filter(cluster == "cluster_2")
  group_by(cluster) %>% 
  summarise(sum_true = sum(asses_random_origin))  %>% 
  as.data.frame() %>% 
  mutate(permutation_fdr = sum_true/1000) %>% 
  select(c(cluster, permutation_fdr)) %>% 
  write.table(file = "results/ldopa/ldopa-cluster-permFDR.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
  