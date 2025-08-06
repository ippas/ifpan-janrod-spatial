# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/statistics-functions-v2.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")

source("preprocessing/functions/statistics-functions/statistics-functions-v3.R")


# load data
print("loading integrated data")
load(file = "results/risperidone-3q29/ris3q29_integrate.RData")

print("loading ris3q29_st_data")
load(file = "results/risperidone-3q29/ris3q29_st_data.RData")


metadata_ris3q29 <- read_metadata(file_path = "data/risperidone-3q29/metadata-3q29-ris.tsv", 
                                  treatments = c("saline", "risperidone")) 



# prepare datafrmae with data parameters
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


data_params_df %>% 
  filter(resolution %in% c(0.1, 0.2, 0.4)) %>% 
  filter(!(data_type_name %in% c("seurat", "quantile_metric"))) -> data_params_df


# visualization test
samples_salWt <- metadata_ris3q29 %>% 
  filter(treatment == "saline" & mouse_genotype == "wtwt") %>%
  .[, 1]

samples_risWt <- metadata_ris3q29 %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wtwt") %>%
  .[, 1]

samples_salDel <- metadata_ris3q29 %>% 
  filter(treatment == "saline" & mouse_genotype == "wtdel") %>%
  filter(sample_ID != "S13839Nr3") %>% 
  .[, 1]

samples_risDel <- metadata_ris3q29 %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wtdel") %>%
  .[, 1]

samples_wt <- metadata_ris3q29 %>% 
  filter(mouse_genotype == "wtwt") %>%
  .[, 1]

samples_del <- metadata_ris3q29 %>% 
  filter(mouse_genotype == "wtdel") %>%
  filter(sample_ID != "S13839Nr3") %>% 
  .[, 1]


summarize_and_test_v3(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df,
  control_samples      = samples_salWt,
  experiment_samples   = samples_salDel,
  metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
  statistic_metrics    = c("mean"),
  mean_threshold       = 0,
  num_cores            = 10
) -> salWtSalDel_summary_statistics

save(salWtSalDel_summary_statistics, file = "results/risperidone-3q29/salWtSalDel_summary_statistics.RData")
rm(salWtSalDel_summary_statistics)

summarize_and_test_v3(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df,
  control_samples      = samples_risWt,
  experiment_samples   = samples_risDel,
  metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
  statistic_metrics    = c("mean"),
  mean_threshold       = 0,
  num_cores            = 10
)  -> risWtRisDel_summary_statistics

save(risWtRisDel_summary_statistics, file = "results/risperidone-3q29/risWtRisDel_summary_statistics.RData")
rm(risWtRisDel_summary_statistics)

summarize_and_test_v3(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df,
  control_samples      = samples_wt,
  experiment_samples   = samples_del,
  metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
  statistic_metrics    = c("mean"),
  mean_threshold       = 0,
  num_cores            = 10
)  -> wtDel_summary_statistics

save(wtDel_summary_statistics, file = "results/risperidone-3q29/wtDel_summary_statistics.RData")
rm(wtDel_summary_statistics)

summarize_and_test_v3(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df,
  control_samples      = samples_salWt,
  experiment_samples   = samples_risWt,
  metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
  statistic_metrics    = c("mean"),
  mean_threshold       = 0,
  num_cores            = 10
) -> risWtSalWt_summary_statistics

save(risWtSalWt_summary_statistics, file = "results/risperidone-3q29/risWtSalWt_summary_statistics.RData")


summarize_and_test_v3(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df,
  control_samples      = samples_salDel,
  experiment_samples   = samples_risDel,
  metrics      = c("mean", "median", "skewness", "kurtosis", "sum"),
  statistic_metrics    = c("mean"),
  mean_threshold       = 0,
  num_cores            = 10
) -> risDelSalDel_summary_statistics

save(risDelSalDel_summary_statistics, file = "results/risperidone-3q29/risDelSalDel_summary_statistics.RData")


# risWtSalWt_summary_statistics$raw_data$resolution_0.05$cluster_0$experiment$sum

# perform_statistical_tests

save(samples_salWt,
     samples_risWt,
     ris3q29_integrate,
     ris3q29_st_data,
     risWtSalWt_summary_statistics, 
     risDelSalDel_summary_statistics,
     file = "results/risperidone-3q29/risWtSalWt-risDelsalDel.RData")

