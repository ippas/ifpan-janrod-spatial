filter_cluster_summary_statistics <- function(data,
                                    gene_vector,
                                    data_type = "quantile_normalize",
                                    resolution = "resolution_0.4",
                                    cluster = "cluster_0",
                                    summary_metric = "sum") {
  # Extract control and experiment data
  control_df <- data[[data_type]][[resolution]][[cluster]]$control[[summary_metric]]
  experiment_df <- data[[data_type]][[resolution]][[cluster]]$experiment[[summary_metric]]
  
  # Combine, convert to df, and extract gene names
  combined_df <- cbind(control_df, experiment_df) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "peak") %>%
    mutate(gene = stringr::word(peak, 3, sep = "-")) %>%
    select(peak, gene, everything()) %>%
    filter(gene %in% gene_vector)
  
  return(combined_df)
}

################################################################################
# summarize_and_test
summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05,
                   num_cores = 10,
                   data_params_df = data_params_df[3,],
                   control_samples = samples_wt,
                   experiment_samples = samples_del,
                   mean_threshold = 0,
                   statistic_metrics = c("mean", "median", "sum"),
                   metrics = c("expression_spot", "mean", "median", "sum"))   -> tmp


tmp$raw_data$resolution_0.4 <- result_all

spatial_generate_cluster_heatmap(
  data            = tmp,
  genes           = genes_del3q29,
  data_type       = "raw_data",
  resolution      = "resolution_0.4",
  summary_metric  = "mean",
  scale_range     = c(-2, 2),
  cluster_name    = "cluster_0",
  metadata        = sample_info,
  group1          = "mouse_genotype",
  group2          = "treatment",
  colors_group1   = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2   = c("saline" = "orange", "risperidone" = "green"),
  title           = "quantile normalize; res0.4, mean, cluster 1", 
  remove_empty_rows = F,
  colors_on_heatmap = c("#6a89b1ff", "white", "#832524ff"),
  z_scale = T,
  show_values = F,
  verbose = T
)

filter_cluster_summary_statistics(
  data = tmp,
  data_type = "raw_data", 
  resolution = "resolution_0.4",
  summary_metric = "sum",
  gene_vector = genes_del3q29
)%>%  head(1)


tmp$raw_data$resolution_0.4$cluster_0$control$expression_spot$S13839Nr1["peak-160966-Bdh1",] -> vector1

################################################################################

summarize_and_test_v2(
  spatial_data         = ris3q29_st_data,
  trim                 = 0.05,
  data_params_df       = data_params_df[3,],
  control_samples      = samples_wt,
  experiment_samples   = samples_del,
  summary_metrics      = c("expression_spot", "mean", "median", "sum"),
  statistic_metrics    = c("mean", "median", "sum"),
  mean_threshold       = 0,
  num_cores            = 16,
  verbose              = TRUE
) -> tmp_v2

tmp_v2$raw_data$resolution_0.4$cluster_0$control$expression_spot$S13839Nr1["peak-160966-Bdh1",] -> vector2
filter_cluster_summary_statistics(
  data = tmp_v2,
  data_type = "raw_data", 
  resolution = "resolution_0.4",
  summary_metric = "sum",
  gene_vector = genes_del3q29
) %>% head(1)


spatial_generate_cluster_heatmap(
  data            = tmp_v2,
  genes           = genes_del3q29,
  data_type       = "raw_data",
  resolution      = "resolution_0.4",
  summary_metric  = "mean",
  scale_range     = c(-2, 2),
  cluster_name    = "cluster_0",
  metadata        = sample_info,
  group1          = "mouse_genotype",
  group2          = "treatment",
  colors_group1   = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2   = c("saline" = "orange", "risperidone" = "green"),
  title           = "quantile normalize; res0.4, mean, cluster 1", 
  remove_empty_rows = F,
  colors_on_heatmap = c("#6a89b1ff", "white", "#832524ff"),
  z_scale = T,
  show_values = F,
  verbose = T
)

tmp$raw_data$resolution_0.4$cluster_2

df1 <- data.frame(cell = names(vector1), val1 = vector1)
df2 <- data.frame(cell = names(vector2), val2 = vector2)

df1 <- data.frame(cell = names(v1), val1 = v1)
df2 <- data.frame(cell = names(v2), val2 = v2)

merged <- merge(df1, df2, by = "cell", all = TRUE)merged <- merge(df1, df2, by = "cell", all = TRUE)


#######################################3
tmp_v3 <- compute_data_summary_v3(
  spatial_data = ris3q29_st_data,
  resolution = data_params_df[3, "resolution"],
  trim = 0.05,
  num_cores = 16,
  control_samples = samples_wt,
  experiment_samples = samples_del,
  data_type = data_params_df[3, "data_type"],
  min_number_spots = 20,  # Można dopasować, jeśli potrzebujesz innego progu
  metrics = c("expression_spot", "mean", "median", "sum")
)

