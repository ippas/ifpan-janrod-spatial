# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")

# executes seurat analysis for ldopa
path_to_data <- "data/ldopa/spaceranger-corrected/"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
sample_names <- list.files(path = path_to_data)

# read metadata for ldopa
metadata_ldopa <- read_metadata(file_path = "data/samples-spatial-metadata.tsv", 
                                      treatments = c("saline", "ldopa"))

# visualization test
samples_saline <- metadata_ldopa %>% 
  filter(treatment == "saline" & experiment == "tif-ldopa") %>%
  .[, 1]

samples_ldopa <- metadata_ldopa %>% 
  filter(treatment == "ldopa") %>%
  # filter(sample_ID != "S5295Nr2") %>% 
  .[, 1]

sample_names <- c(samples_saline, samples_ldopa)

# read peaks for ldopa
info_peaks_ldopa <- read_gene_annotation(file_path = "data/ldopa/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for ldopa
ldopa_integrate <- integrate_data_seurat(path_to_data, nfeatures, dims)

# read images to spatial transcriptoms for ldopa
images_ldopa <- create_images_tibble(path_to_data, sample_names)

# read barcode data for ldopa
barcode_ldopa <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_ldopa
)

# poprawić do 
ldopa_st_data <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_ldopa,
                                                barcode_info = barcode_ldopa,
                                                images_info = images_ldopa,
                                                integrated_data = ldopa_integrate,
                                                peaks_info = info_peaks_ldopa)

ldopa_st_data <- add_seurat_data(spatial_data = ldopa_st_data,
                                            integrated_data = ldopa_integrate)

ldopa_st_data <- add_filtered_data(spatial_data = ldopa_st_data,
                                              mean_expression_threshold = 0.5)

ldopa_st_data <- add_colfilt_data(spatial_data = ldopa_st_data,
                                             min_spot_threshold = 0,
                                             expression_threshold = 2)

ldopa_st_data <- add_range_normalize_data(spatial_data = ldopa_st_data,
                                                     range = 1500,
                                                     flatten = 1,
                                                     threshold = 500)

ldopa_st_data <- add_clusters_data(spatial_data = ldopa_st_data,
                                              integrated_data = ldopa_integrate,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

ldopa_st_data <- evaluate_clustering_stability(spatial_data = ldopa_st_data,
                                                          seurat_object = ldopa_integrate,
                                                          resolution_start = 0.05,
                                                          resolution_end = 2,
                                                          resolution_step = 0.05)




for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8, 1)){
  ldopa_st_data <-
    add_quantile_norm_data(
      spatial_data = ldopa_st_data,
      resolution = {{resolution}},
      num_cores = 24,
      data_type = "raw_data"
    )
}


ldopa_st_data$stability_results %>%
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
    "quantile_normalize_resolution_0.1",
    "quantile_normalize_resolution_0.15",
    "quantile_normalize_resolution_0.2",
    "quantile_normalize_resolution_0.4",
    "quantile_normalize_resolution_0.8",
    "quantile_normalize_resolution_1"
  )
),
resolution = rep(c(0.1, 0.15, 0.2, 0.4, 0.8, 1), 5),
data_type_name = c(
  rep("raw_data", 6),
  rep("quantile_metric", 6),
  rep("range_normalize", 6),
  rep("seurat", 6),
  rep("quantile_normalize", 6)
)) %>%
  mutate(quantile_normalization = ifelse(data_type_name == "quantile_metric", TRUE, FALSE)) %>% 
  filter(data_type_name != "quantile_metric") -> data_params_df


# prepare sample names to statistics analysis
samples_saline <- metadata_ldopa %>% 
  filter(treatment == "saline" & mouse_genotype == "tif-mutant") %>%
  .[, 1]

samples_ldopa <- metadata_ldopa %>% 
  filter(treatment == "ldopa") %>%
  filter(sample_ID != "S5295Nr2") %>% 
  .[, 1]



summarize_and_test(spatial_data = ldopa_st_data,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_ldopa,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> ldopa_summary_statistics

save(samples_saline,
     samples_ldopa,
     ldopa_integrate,
     ldopa_st_data,
     ldopa_summary_statistics, 
     file = "results/ldopa/ldopa.RData")



DimPlot(ldopa_integrate, reduction = "umap",
        split.by = "sample", ncol = 4)

DimPlot(ldopa_integrate, reduction = "umap")

# DimPlot(object = ldopa_integrate, reduction = "umap", group.by = )

# 
# visualize clusters
spatial_cluster(spatial_data = ldopa_st_data,
                resolution = 0.4,
                samples = c(samples_saline, samples_ldopa),
                palette = palette_allen,
                size= 1.0,
                ncol = 4)
# 
# spatial_interest_cluster(cluster = 3,
#                          # seurat_object = integrated_analysis,
#                          spatial_data = ldopa_st_data,
#                          resolution = 0.8,
#                          samples = c(samples_saline, samples_ldopa),
#                          size= 1,
#                          ncol = 4)

# 
spatial_gene_plot(spatial_data = ldopa_st_data,
                  type_data = "quantile_normalize",
                  gene = "Itpk1",
                  samples =  c(samples_saline, samples_ldopa),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)
# 
spatial_gene_plot(spatial_data = ldopa_st_data,
                  type_data = "raw_data",
                  gene = "Sgk1",
                  samples =  c(samples_saline, samples_ldopa),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 1.4,
                  ncol = 4,
                  tif_image = F,
                  normalization = T)

spatial_feature_plot_cluster

# Custom palette function: gray to red
my_custom_palette <- function(n) {
  colors <- colorRampPalette(c("#EEEEEE", "#ff0000", "#bf0000", "#800000", "#400000"))(n)
  return(colors)
}


################################################################################
# reduce data for ldopa
ldopa_st_data[c("samples", "sample_information", "bcs_inforamtion", "images_information",  "clusters", "quantile_normalize_resolution_0.8", "quantile_normalize_1")] -> ldopa_st_data_reduced
ldopa_summary_statistics$quantile_normalize[c("resolution_0.8", "resolution_1")] -> ldopa_summary_statistics_reduced


save(samples_saline,
     samples_ldopa,
     ldopa_integrate,
     ldopa_st_data_reduced,
     ldopa_summary_statistics_reduced,
     file = "results/ldopa/ldopa-reduced.RData")


ldopa_st_data %>% names
  
ldopa_st_data[c("samples", "sample_information", "bcs_inforamtion", "i2mages_information", "raw_data", "clusters", "quantile_normalize_resolution_0.8", "quantile_normalize_1")] %>% object.size()

ldopa_summary_statistics$quantile_normalize[c("resolution_0.8", "resolution_1")] %>% object.size()
