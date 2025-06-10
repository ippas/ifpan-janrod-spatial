# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")

# executes seurat analysis for risperidone
path_to_data <- "data/risperidone-3q29/spaceranger-corrected/"
nfeatures <- 2000
dims <- 1:30

#setwd("/home/mateusz/projects/ifpan-janrod-spatial/")
#getwd()

# read metadata for risperidone
metadata_ris3q29 <- read_metadata(file_path = "data/risperidone-3q29/metadata-3q29-ris.tsv", 
                                           treatments = c("saline", "risperidone")) 

# create sample_names vector contain id samples to analysis
sample_names <- metadata_ris3q29$sample_ID

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

# read peaks for risperidone
info_peaks_ris3q29 <- read_gene_annotation(file_path = "data/risperidone-3q29/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for risperidone
ris3q29_integrate <- integrate_data_seurat(path_to_data, samples = metadata_ris3q29$sample_ID, nfeatures, dims)

# read images to spatial transcriptoms for risperidone
images_ris3q29 <- create_images_tibble(path_to_data, sample_names)


# read barcode data for risperidone
barcode_ris3q29 <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_ris3q29,
  spaceranger_version = "3.1.3"
)

################################################################################
ris3q29_st_data <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_ris3q29,
                                                barcode_info = barcode_ris3q29,
                                                images_info = images_ris3q29,
                                                integrated_data = ris3q29_integrate,
                                                peaks_info = info_peaks_ris3q29)

ris3q29_st_data <- add_seurat_data(spatial_data = ris3q29_st_data,
                                            integrated_data = ris3q29_integrate)

ris3q29_st_data <- add_filtered_data(spatial_data = ris3q29_st_data,
                                              mean_expression_threshold = 0.5)

ris3q29_st_data <- add_colfilt_data(spatial_data = ris3q29_st_data,
                                             min_spot_threshold = 0,
                                             expression_threshold = 2)

ris3q29_st_data <- add_range_normalize_data(spatial_data = ris3q29_st_data,
                                                     range = 1500,
                                                     flatten = 1,
                                                     threshold = 500)

ris3q29_st_data <- add_clusters_data(spatial_data = ris3q29_st_data,
                                              integrated_data = ris3q29_integrate,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

ris3q29_st_data <- evaluate_clustering_stability(spatial_data = ris3q29_st_data,
                                                          seurat_object = ris3q29_integrate,
                                                          resolution_start = 0.05,
                                                          resolution_end = 2,
                                                          resolution_step = 0.05)

################################################################################
# POPRAWIĆ NAZWY

for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)){
  ris3q29_st_data <-
    add_quantile_norm_data(
      spatial_data = ris3q29_st_data,
      resolution = {{resolution}},
      num_cores = 16,
      data_type = "raw_data"
    )
}


ris3q29_st_data$stability_results %>%
  ggplot(aes(x = resolution, y = silhouette_score)) +
  geom_line()



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



summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05, 
                   num_cores = 16,
                   data_params_df = data_params_df,
                   control_samples = samples_salWt,
                   experiment_samples = samples_salDel,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> salWtSalDel_summary_statistics

save(salWtSalDel_summary_statistics, file = "results/risperidone-3q29/salWtSalDel_summary_statistics.RData")
rm(salWtSalDel_summary_statistics)


summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05, 
                   num_cores = 16,
                   data_params_df = data_params_df,
                   control_samples = samples_risWt,
                   experiment_samples = samples_risDel,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> risWtRisDel_summary_statistics

save(risWtRisDel_summary_statistics, file = "results/risperidone-3q29/risWtRisDel_summary_statistics.RData")
rm(risWtRisDel_summary_statistics)


summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05,
                   num_cores = 16,
                   data_params_df = data_params_df,
                   control_samples = samples_wt,
                   experiment_samples = samples_del,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> wtDel_summary_statistics

save(wtDel_summary_statistics, file = "results/risperidone-3q29/wtDel_summary_statistics.RData")
rm(wtDel_summary_statistics)

summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05, 
                   num_cores = 16,
                   data_params_df = data_params_df,
                   control_samples = samples_salWt,
                   experiment_samples = samples_risWt,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> risWtSalWt_summary_statistics

save(risWtSalWt_summary_statistics, file = "results/risperidone-3q29/risWtSalWt_summary_statistics.RData")

summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05, 
                   num_cores = 16,
                   data_params_df = data_params_df,
                   control_samples = samples_salDel,
                   experiment_samples = samples_risDel,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("mean", "median", "skewness", "kurtosis", "sum")) -> risDelSalDel_summary_statistics

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

