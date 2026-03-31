load("results/risperidone/risperidone-half.RData")

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

summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   mean_threshold = 0,
                   min_number_spots = 15,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> risperidone_summary_statistics_half_minNumberSpot15

# risperidone_summary_statistics_half_minNumberSpot15$raw_data$resolution_0.15$

save(samples_saline,
     samples_risperidone,
     risperidone_integrate_half,
     risperidone_st_data_half,
     risperidone_summary_statistics_half_minNumberSpot15, 
     file = "results/risperidone/risperidone-half_minNumberSpot15.RData")



summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   mean_threshold = 0,
                   min_number_spots = 20,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> risperidone_summary_statistics_half_minNumberSpot20

# risperidone_summary_statistics_half_minNumberSpot15$raw_data$resolution_0.15$

# save(samples_saline,
#      samples_risperidone,
#      risperidone_integrate_half,
#      risperidone_st_data_half,
#      risperidone_summary_statistics_half_minNumberSpot15, 
#      file = "results/risperidone/risperidone-half_minNumberSpot15.RData")

save(samples_saline,
     samples_risperidone,
     risperidone_integrate_half,
     risperidone_st_data_half,
     risperidone_summary_statistics_half_minNumberSpot20,
     file = "results/risperidone/risperidone-half_minNumberSpot_17.03.2026.RData")


# summarize_and_test_v2(
#   spatial_data         = ris3q29_st_data,
#   trim                 = 0.05,
#   data_params_df       = data_params_df,
#   control_samples      = samples_salDel,
#   experiment_samples   = samples_risDel,
#   summary_metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
#   statistic_metrics    = c("mean", "median", "sum"),
#   mean_threshold       = 0,
#   num_cores            = 16,
#   min_number_spots     = 15,
#   verbose              = TRUE
# ) -> risDelSalDel_summary_statistics


compute_data_summary_v2
