ris3q29_summary_global_expression <- compute_global_expression_multi(
  spatial_data       = ris3q29_st_data,
  data_type          = c("raw_data", "quantile_normalize_resolution_0.4"),
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  trim = 0,
  metric             = c("mean", "sum"),
  verbose            = TRUE
)


# ustaw ścieżkę i nazwę pliku
svg("results/risperidone-3q29/figures/heatmap-deletion-genes/heatmapGlobal_quantileNormalize04_meanExpression.svg", width = 8, height = 6)

# wywołanie funkcji generującej wykres
spatial_generate_global_heatmap(
  data               = ris3q29_summary_global_expression,
  genes              = genes_del3q29,
  summary_metric     = "mean",
  data_type          = "quantile_normalize_resolution_0.4",
  metadata           = sample_info,
  group1             = "mouse_genotype",
  group2             = "treatment",
  colors_group1      = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2      = c("saline" = "orange", "risperidone" = "green"),
  title              = "Bulk expression; deletion genes; quantile normalize; mean",
  scale_range        = c(-2, 2),
  remove_empty_rows  = FALSE,
  colors_on_heatmap  = c("#6a89b1ff", "white", "#832524ff"),
  z_scale            = TRUE,
  show_values        = FALSE,
  threshold_expr     = 0.01, 
  select_max_peak    = FALSE,
  verbose            = TRUE
)

# zakończ zapis
dev.off()


################################################################################
# heatmap per clusters
summarize_and_test(spatial_data = ris3q29_st_data,
                   trim = 0.05,
                   num_cores = 16,
                   data_params_df = data_params_df[9,],
                   control_samples = samples_wt,
                   experiment_samples = samples_del,
                   mean_threshold = 0,
                   statistic_metrics = c("mean"),
                   metrics = c("mean", "median", "sum")) -> tmp

sample_info <- ris3q29_st_data$sample_information
# result_summary_statistics_v3$quantile_normalize$resolution_0.4$cluster_20 <- NULL


"results/risperidone-3q29/figures/heatmap-deletion-genes/perCluster"
spatial_generate_cluster_heatmap(
  data            = tmp,
  # data = tmp,
  genes           = genes_del3q29,
  data_type       = "quantile_normalize",
  resolution      = "resolution_0.4",
  summary_metric  = "mean",
  scale_range     = c(-2, 2),
  cluster_name    = "cluster_0",
  metadata        = sample_info,
  group1          = "mouse_genotype",
  group2          = "treatment",
  colors_group1   = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2   = c("saline" = "orange", "risperidone" = "green"),
  title           = "quantile normalize; res0.4, mean, cluster 0", 
  remove_empty_rows = T,
  colors_on_heatmap = c("#6a89b1ff", "white", "#832524ff"),
  z_scale = TRUE,
  show_values = FALSE,
  verbose = T
)

output_dir <- "results/risperidone-3q29/figures/heatmap-deletion-genes/perCluster"
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 0:19) {
  cluster_name <- paste0("cluster_", i)
  plot_title   <- paste0("quantile normalize; res0.4, mean, ", cluster_name)
  file_name    <- paste0("quantileNormalize_res0.4_mean_cluster", i, ".png")
  file_path    <- file.path(output_dir, file_name)
  
  message(">>> Generating heatmap for: ", cluster_name, " -> ", file_name)
  
  png(filename = file_path, width = 2400, height = 2000, res = 300)  # 300 dpi
  
  spatial_generate_cluster_heatmap(
    data               = tmp,
    genes              = genes_del3q29,
    data_type          = "quantile_normalize",
    resolution         = "resolution_0.4",
    summary_metric     = "mean",
    scale_range        = c(-2, 2),
    cluster_name       = cluster_name,
    metadata           = sample_info,
    group1             = "mouse_genotype",
    group2             = "treatment",
    colors_group1      = c("wtwt" = "gray", "wtdel" = "blue"),
    colors_group2      = c("saline" = "orange", "risperidone" = "green"),
    title              = plot_title,
    remove_empty_rows  = TRUE,
    colors_on_heatmap  = c("#6a89b1ff", "white", "#832524ff"),
    z_scale            = TRUE,
    show_values        = FALSE,
    verbose            = FALSE
  )
  
  dev.off()
}

