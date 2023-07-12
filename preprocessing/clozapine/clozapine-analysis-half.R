# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")


# executes seurat analysis for clozapine
path_to_data <- "data/clozapine/spaceranger-corrected-half//"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
# sample_names <- list.files(path = path_to_data)

# read metadata for clozapine
metadata_clozapine <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                 treatments = c("saline", "clozapine"))

# read peaks for clozapine
info_peaks_clozapine <- read_gene_annotation(file_path = "data/clozapine/gene-annotation/peaks-annotate-reduction.tsv")


# visualization test
samples_saline <- metadata_clozapine %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_clozapine <- metadata_clozapine %>% 
  filter(treatment == "clozapine" & mouse_genotype == "wt") %>%
  .[, 1]

sample_names <- c(samples_saline, samples_clozapine)


# executes seurat analysis for clozapine
clozapine_integrate_half <- integrate_data_seurat(path_to_data, nfeatures, dims)

# read images to spatial transcriptoms for clozapine
images_clozapine_half <- create_images_tibble(path_to_data, sample_names)

# read barcode data for clozapine
barcode_clozapine_half <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_clozapine_half
)

# poprawić do 
clozapine_st_data_half <- create_spatial_data(sample_names = sample_names,
                                           metadata = metadata_clozapine,
                                           barcode_info = barcode_clozapine_half,
                                           images_info = images_clozapine_half,
                                           integrated_data = clozapine_integrate_half,
                                           peaks_info = info_peaks_clozapine)

clozapine_st_data_half <- add_seurat_data(spatial_data = clozapine_st_data_half,
                                       integrated_data = clozapine_integrate_half)

clozapine_st_data_half <- add_filtered_data(spatial_data = clozapine_st_data_half,
                                         mean_expression_threshold = 0.5)

clozapine_st_data_half <- add_colfilt_data(spatial_data = clozapine_st_data_half,
                                        min_spot_threshold = 0,
                                        expression_threshold = 2)

clozapine_st_data_half <- add_range_normalize_data(spatial_data = clozapine_st_data_half,
                                                range = 1500,
                                                flatten = 1,
                                                threshold = 500)

clozapine_st_data_half <- add_clusters_data(spatial_data = clozapine_st_data_half,
                                         integrated_data = clozapine_integrate_half,
                                         resolution_start = 0.05,
                                         resolution_end = 2,
                                         resolution_step = 0.05)

clozapine_st_data_half <- evaluate_clustering_stability(spatial_data = clozapine_st_data_half,
                                                     seurat_object = clozapine_integrate_half,
                                                     resolution_start = 0.05,
                                                     resolution_end = 2,
                                                     resolution_step = 0.05)

for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)){
  clozapine_st_data_half <-
    add_quantile_norm_data(
      spatial_data = clozapine_st_data_half,
      resolution = {{resolution}},
      num_cores = 24,
      data_type = "raw_data"
    )
}


clozapine_st_data_half$stability_results %>%
  ggplot(aes(x = resolution, y = silhouette_score)) +
  geom_line()


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


summarize_and_test(spatial_data = clozapine_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_clozapine,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> clozapine_summary_statistics_half


save(samples_saline,
     samples_clozapine,
     clozapine_integrate_half,
     clozapine_st_data_half,
     clozapine_summary_statistics_half, 
     file = "results/risperidone/clozapine-half.RData")
