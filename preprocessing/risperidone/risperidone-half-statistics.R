
resolution <- 0.8

investigate_summary_cluster(data = tmp, metric = "mean", gene = "Fosb")


data.frame(data_type = c(
             rep("raw_data", 3),
             rep("range_normalize", 3),
             c(
               "quantile_normalize_resolution_0.1",
               "quantile_normalize_resolution_0.4",
               "quantile_normalize_resolution_0.8"
             )
           ),
           resolution = rep(c(0.1, 0.4, 0.8), 3),
           data_type_name = c(
             rep("raw_data", 3),
             rep("range_normalize", 3),
             rep("quantile_normalize", 3)
           )) -> data_params_df


data.frame(data_type = c(
  rep("raw_data", 6),
  rep("raw_data", 6),
  rep("range_normalize", 6),
  rep("seurat", 6),
  c(
    "quantile_normalize_resolution_0.05",
    "quantile_normalize_resolution_0.1",
    "quantile_normalize_resolution_0.15",
    "quantile_normalize_resolution_0.2",
    "quantile_normalize_resolution_0.4",
    "quantile_normalize_resolution_0.8"
  )
),
resolution = rep(c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8), 5),
data_type_name = c(
  rep("raw_data", 6),
  rep("quantile_metric", 6),
  rep("range_normalize", 6),
  rep("seurat", 6),
  rep("quantile_normalize", 6)
)) %>%
  mutate(quantile_normalization = ifelse(data_type_name == "quantile_metric", TRUE, FALSE)) -> data_params_df
  

summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   mean_threshold = 0,
                   metrics = c("mean", "median")) -> risperidone_summary_statistics



summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone[-3],
                   mean_threshold = 0,
                   metrics = c("mean", "median")) -> risperidone_summary_statistics_remove_ris

##################3

compute_data_summary(spatial_data = risperidone_st_data_half,
                     resolution = 0.8,
                     trim = 0.05,
                     num_cores = 24,
                     control_samples = samples_saline,
                     experiment_samples = samples_risperidone,
                     data_type = "raw_data",
                     metrics = c("expression_spot", "mean")) -> quantile_summary_data_0.8

quantile_summary_data_0.8$cluster_20$experiment$mean

# 
 
perform_statistical_tests(spatial_data = risperidone_st_data_half,
                          summary_data = quantile_summary_data,
                          metric = "mean",
                          resolution = resolution,
                          num_cores = 24,
                          mean_threshold = 0.2,
                          quantile_normalization = T) -> quantile_summary_data_0.8


tmp <-"risperidone_summary_statistics"
tmp <- "risperidone_summary_statistics_remove_ris"

filter_data_statistics(summary_data = get(tmp), 
                       data_type = "range_normalize", 
                       resolution = 0.1,
                       metric = "mean",
                       control_mean_threshold = 0.5,
                       experiment_mean_threshold = 0.5,
                       log2ratio_threshold = 0.5,
                       t_test_threshold = 0.05)




summary_filter_statistics(data = risperidone_filter_statistics)


risperidone_summary_statistics$raw_data$resolution_0.1$cluster_0$control$mean[2,]




spatial_interest_cluster(cluster = 5,
                         # seurat_object = integrated_analysis,
                         spatial_data = risperidone_st_data_half,
                         resolution = 0.4,
                         samples = c(samples_saline, samples_risperidone),
                         size= 1,
                         ncol = 4)

spatial_gene_plot(spatial_data = risperidone_st_data_half,
                  type_data = "quantile_normalize",
                  gene = "Homer1",
                  samples =  c(samples_saline[-1], samples_risperidone[-3]),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 5,
                  normalization = T) 