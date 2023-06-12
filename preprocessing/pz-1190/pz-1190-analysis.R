# read function
source("preprocessing/functions/functions-spatial-data.R")


# executes seurat analysis for pz1190
path_to_data <- "data/pz-1190/spaceranger-corrected/"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
sample_names <- list.files(path = path_to_data)

# read metadata for pz1190
metadata_pz1190 <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "pz_1190"))

# read peaks for pz1190
info_peaks_pz1190 <- read_gene_annotation(file_path = "data/pz-1190/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for pz1190
pz1190_integrate_half <- integrate_data_seurat(path_to_data, nfeatures, dims)

# read images to spatial transcriptoms for pz1190
images_pz1190_half <- create_images_tibble(path_to_data, sample_names)

# read barcode data for pz1190
barcode_pz1190_half <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_pz1190_half
)

# poprawić do 
pz1190_st_data_half <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_pz1190,
                                                barcode_info = barcode_pz1190_half,
                                                images_info = images_pz1190_half,
                                                integrated_data = pz1190_integrate_half,
                                                peaks_info = info_peaks_pz1190)

pz1190_st_data_half <- add_seurat_data(spatial_data = pz1190_st_data_half,
                                            integrated_data = pz1190_integrate_half)

pz1190_st_data_half <- add_filtered_data(spatial_data = pz1190_st_data_half,
                                              mean_expression_threshold = 0.5)

pz1190_st_data_half <- add_colfilt_data(spatial_data = pz1190_st_data_half,
                                             min_spot_threshold = 0,
                                             expression_threshold = 2)

pz1190_st_data_half <- add_range_normalize_data(spatial_data = pz1190_st_data_half,
                                                     range = 1500,
                                                     flatten = 1,
                                                     threshold = 500)

pz1190_st_data_half <- add_clusters_data(spatial_data = pz1190_st_data_half,
                                              integrated_data = pz1190_integrate_half,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

pz1190_st_data_half <- evaluate_clustering_stability(spatial_data = pz1190_st_data_half,
                                                          seurat_object = pz1190_integrate_half,
                                                          resolution_start = 0.05,
                                                          resolution_end = 2,
                                                          resolution_step = 0.05)





for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)){
  pz1190_st_data_half <-
    add_quantile_norm_data(
      spatial_data = pz1190_st_data_half,
      resolution = {{resolution}},
      num_cores = 24,
      data_type = "raw_data"
    )
}


pz1190_st_data_half$stability_results %>%
  ggplot(aes(x = resolution, y = silhouette_score)) +
  geom_line()

# visualization test
samples_saline <- metadata_pz1190 %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_pz1190 <- metadata_pz1190 %>% 
  filter(treatment == "pz_1190" & mouse_genotype == "wt") %>%
  .[, 1]



spatial_gene_plot(spatial_data = pz1190_st_data_half,
                  type_data = "raw_data",
                  gene = "Egr4",
                  samples =  c(samples_saline, samples_pz1190),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)

# statistics pz1190
summarize_and_test(spatial_data = pz1190_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_pz1190,
                   mean_threshold = 0,
                   metrics = c("mean", "median")) -> pz1190_summary_statistics



summarize_and_test(spatial_data = pz1190_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df[8,],
                   control_samples = samples_saline,
                   experiment_samples = samples_pz1190,
                   mean_threshold = 0,
                   metrics = c("mean", "median")) -> tmp




