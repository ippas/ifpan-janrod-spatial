# read function
source("preprocessing/function-spatial-data.R")



# create sample_names vector contain id samples to analysis
sample_names <- list.files(path = path_to_data)

# executes seurat analysis for risperidone
path_to_data <- "data/risperidone/spaceranger-corrected//"
nfeatures <- 2000
dims <- 1:30

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))

# read peaks for risperidone
info_peaks_risperidone <- read_gene_annotation(file_path = "data/risperidone/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for risperidone
risperidone_integrate <- integrate_data_seurat(path_to_data, nfeatures, dims)

# read images to spatial transcriptoms for risperidone
images_risperidone <- create_images_tibble(path_to_data, sample_names)

# read barcode data for risperidone
barcode_risperidone <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_risperidone
)


risperidone_st_data <- create_spatial_data(sample_names = sample_names,
                                           metadata = metadata_risperidone,
                                           barcode_info = barcode_risperidone,
                                           images_info = images_risperidone,
                                           integrated_data = risperidone_integrate,
                                           peaks_info = info_peaks_risperidone)

risperidone_st_data <- add_seurat_data(spatial_data = risperidone_st_data,
                                       integrated_data = risperidone_integrate)

risperidone_st_data <- add_filtered_data(spatial_data = risperidone_st_data,
                                         mean_expression_threshold = 0.5)

risperidone_st_data <- add_colfilt_data(spatial_data = risperidone_st_data,
                                        min_spot_threshold = 0,
                                        expression_threshold = 2)

risperidone_st_data <- add_range_normalize_data(spatial_data = risperidone_st_data,
                                                range = 1500,
                                                flatten = 1,
                                                threshold = 500)

risperidone_st_data <- add_clusters_data(spatial_data = risperidone_st_data,
                                         integrated_data = risperidone_integrate,
                                         resolution_start = 0.1,
                                         resolution_end = 2,
                                         resolution_step = 0.1)

risperidone_st_data <- evaluate_clustering_stability(spatial_data = risperidone_st_data,
                                                     seurat_object = risperidone_integrate,
                                                     resolution_start = 0.1,
                                                     resolution_end = 2,
                                                     resolution_step = 0.1)

risperidone_st_data$stability_results %>%
  ggplot(aes(x = resolution, y = silhouette_score)) +
  geom_line()


# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]

resolution = 0.8

DimPlot(risperidone_integrate_half, reduction = "umap", 
        split.by = "sample", ncol = 4)

DimPlot(risperidone_integrate_half, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(spatial_data = risperidone_st_data,
                resolution = 0.8,
                samples = c(samples_saline, samples_risperidone),
                palette = palette_allen, 
                size= 1.0, 
                ncol = 4)

# check number of barcode per cluster
calculate_spots_per_cluster_matrix(seurat_object = risperidone_integrate,
                                   spatial_data = risperidone_st_data,
                                   resolution = 0.5) -> spots_per_cluster

{spots_per_cluster %>% select(-c(sample, treatment)) < 30} %>% 
  apply(., 2, sum)


# visualize interest cluster
spatial_interest_cluster(cluster = 5,
                         # seurat_object = integrated_analysis,
                         spatial_data = risperidone_st_data,
                         resolution = 0.5,
                         samples = c(samples_saline, samples_risperidone),
                         size= 1,
                         ncol = 4)


# visualize gene 
spatial_gene_plot(spatial_data = risperidone_st_data,
                  type_data = "raw_data",
                  gene = "Apod",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)

### developed
# statistics calculation
calculate_gene_expression_stats(spatial_data = risperidone_st_data,
                                stat_test = "t.test",
                                data_type = "raw_data",
                                resolution = 0.5,
                                expression_unit = "sample",
                                control_samples = samples_saline,
                                experiment_samples = samples_risperidone,
                                save_results = F
                                ) -> tmp



tmp %>% 
  as.data.frame %>% 
  apply(., 2, as.numeric) %>%
  as.data.frame() %>%
  filter(as.numeric(cluster_4) < 0.01)


