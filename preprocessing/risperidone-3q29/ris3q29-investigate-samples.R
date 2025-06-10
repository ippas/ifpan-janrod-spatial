samples_salDel <- metadata_ris3q29 %>% 
  filter(treatment == "saline" & mouse_genotype == "wtdel") %>%
  # filter(sample_ID != "S13839Nr26") %>% 
  .[, 1]

samples_risDel <- metadata_ris3q29 %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wtdel") %>%
  # filter(sample_ID != "S13839Nr23") %>% 
  .[, 1]


# 
summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05,
                   num_cores = 24,
                   data_params_df = {data_params_df %>% filter(data_type == "quantile_normalize_resolution_0.4")},
                   control_samples = samples_salDel,
                   experiment_samples = samples_risDel,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> tmp


data_params_df %>% 
  filter(data_type == "quantile_normalize_resolution_0.4")


tmp$quantile_normalize$resolution_0.4$cluster_0$control$mean %>% dim
tmp$quantile_normalize$resolution_0.4$cluster_0$experiment 

tmp$quantile_normalize$resolution_0.4 %>% 
  lapply(., function(x){
    cbind(x$control$mean, x$experiment$mean) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "peak_gene") %>% 
      filter(peak_gene %in% c("peak-422102-Cck", "peak-179117-Rps18", "peak-262621-Txnip")) %>% 
      select(c(peak_gene, S13839Nr26, S13839Nr23))
  }) %>% .[c("cluster_16", "cluster_13", "cluster_15")]



tmp$quantile_normalize$resolution_0.4 %>%
  .[c("cluster_10", "cluster_13", "cluster_15", "cluster_16", "cluster_18")] %>%
  { map2_dfr(.x = ., .y = names(.), ~ 
               cbind(.x$control$mean, .x$experiment$mean) %>%
               as.data.frame() %>%
               rownames_to_column("peak_gene") %>%
               mutate(cluster = .y)
  )
  } %>%
  filter(
    (peak_gene == "peak-422102-Cck"    & cluster == "cluster_16") |
      (peak_gene == "peak-179117-Rps18"  & cluster == "cluster_13") |
      (peak_gene == "peak-262621-Txnip"  & cluster == "cluster_15") |
      (peak_gene == "peak-60449-Pttg1"   & cluster == "cluster_13") |
      (peak_gene == "peak-297366-Errfi1" & cluster == "cluster_10") |
      (peak_gene == "peak-327389-Glcci1" & cluster == "cluster_18")
  ) %>% 
  select(c(peak_gene, cluster, everything())) 
  
  
