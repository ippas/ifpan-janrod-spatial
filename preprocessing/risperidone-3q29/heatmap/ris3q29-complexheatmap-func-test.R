spatial_generate_cluster_heatmap(
  data            = wtDel_summary_statistics,
  genes           = genes_del3q29,
  data_type       = "raw_data",
  resolution      = "resolution_0.4",
  summary_metric  = "mean",
  scale_range     = c(-4, 4),
  cluster_name    = "cluster_5",
  metadata        = sample_info,
  group1          = "mouse_genotype",
  group2          = "treatment",
  colors_group1   = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2   = c("saline" = "orange", "risperidone" = "green"),
  title           = "quantile normalize; res0.4, mean, cluster 1", 
  remove_empty_rows = FALSE,
  colors_on_heatmap = c("#6a89b1ff", "white", "#832524ff"),
  z_scale = TRUE,
  show_values = FALSE,
  verbose = T
)


summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0,
                   num_cores = 16,
                   data_params_df = data_params_df[5,],
                   control_samples = samples_wt,
                   experiment_samples = samples_del,
                   mean_threshold = 0,
                   statistic_metrics = c("mean"),
                   metrics = c("mean")) -> tmp

ris3q29_st_data$clusters %>% 
  filter(cluster_resolution_0.4 == 4) %>%
  select(sample, barcode, experiment, mouse_genotype, treatment, cluster_resolution_0.4) %>% dim


ris3q29_summary_global_expression$quantile_normalize_resolution_0.4 <- compute_global_expression_summary(
  spatial_data       = ris3q29_st_data,
  data_type          = "quantile_normalize_resolution_0.4",
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  trim = 0,
  metric             = c("mean", "sum"),
  verbose            = TRUE
)



ris3q29_summary_global_expression <- compute_global_expression_multi(
  spatial_data       = ris3q29_st_data,
  data_type          = c("raw_data",  "quantile_normalize_resolution_0.1",  
                         "quantile_normalize_resolution_0.2", "quantile_normalize_resolution_0.4"),
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  trim = 0,
  metric             = c("mean", "sum"),
  verbose            = TRUE
)



spatial_generate_global_heatmap(
  data               = ris2q29_summary_global_expression,
  genes              = genes_del3q29,
  data_type          = "quantile_normalize_resolution_0.4",
  summary_metric     = "mean",
  metadata           = sample_info,
  group1             = "mouse_genotype",
  group2             = "treatment",
  colors_group1      = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2      = c("saline" = "orange", "risperidone" = "green"),
  title              = "quantile normalize; mean",
  scale_range        = c(-2, 2),
  remove_empty_rows  = F,
  colors_on_heatmap  = c("#6a89b1ff", "white", "#832524ff"),
  z_scale            = TRUE,
  show_values        = FALSE,
  verbose            = TRUE
)

