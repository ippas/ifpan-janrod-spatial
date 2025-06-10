# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")

# executes seurat analysis for risperidone
path_to_data <- "data/risperidone/spaceranger-peakWithGene//"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
sample_names <- list.files(path = path_to_data)

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))


# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]

# read peaks for risperidone
info_peaks_risperidone <- read_gene_annotation(file_path = "data/risperidone/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for risperidone
risperidone_integrate_half <- integrate_data_seurat(path_to_data, sample_names, nfeatures, dims)

# read images to spatial transcriptoms for risperidone
images_risperidone_half <- create_images_tibble(path_to_data, sample_names)

# read barcode data for risperidone
barcode_risperidone_half <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_risperidone_half,
  spaceranger_version = "3.1.3"
)

# poprawić do 
risperidone_st_data_half <- create_spatial_data(sample_names = sample_names,
                                                metadata = metadata_risperidone,
                                                barcode_info = barcode_risperidone_half,
                                                images_info = images_risperidone_half,
                                                integrated_data = risperidone_integrate_half,
                                                peaks_info = info_peaks_risperidone)

risperidone_st_data_half <- add_seurat_data(spatial_data = risperidone_st_data_half,
                                            integrated_data = risperidone_integrate_half)

risperidone_st_data_half <- add_filtered_data(spatial_data = risperidone_st_data_half,
                                              mean_expression_threshold = 0.5)

risperidone_st_data_half <- add_colfilt_data(spatial_data = risperidone_st_data_half,
                                             min_spot_threshold = 0,
                                             expression_threshold = 2)

risperidone_st_data_half <- add_range_normalize_data(spatial_data = risperidone_st_data_half,
                                                     range = 1500,
                                                     flatten = 1,
                                                     threshold = 500)

risperidone_st_data_half <- add_clusters_data(spatial_data = risperidone_st_data_half,
                                              integrated_data = risperidone_integrate_half,
                                              resolution_start = 0.05,
                                              resolution_end = 2,
                                              resolution_step = 0.05)

risperidone_st_data_half <- evaluate_clustering_stability(spatial_data = risperidone_st_data_half,
                                                          seurat_object = risperidone_integrate_half,
                                                          resolution_start = 0.05,
                                                          resolution_end = 2,
                                                          resolution_step = 0.05)



# save to tsv expression and coordinates

# save metadata
risperidone_st_data_half$sample_information %>% 
  write.table("data/risperidone/ifpas-jackcraig-spatialdata/risSal-metadata.tsv", 
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE
             )
# save coordinates of barcodes
# baza danych z przypisaniem do klastrów
clusters_df <- risperidone_st_data_half$clusters %>%
  select(sample, barcode, cluster = cluster_resolution_0.4) %>%
  mutate(cluster = paste0("cluster_", cluster))

# ścieżka bazowa
# baza danych z przypisaniem do klastrów
clusters_df <- risperidone_st_data_half$clusters %>%
  select(sample, barcode, cluster = cluster_resolution_0.4) %>%
  mutate(cluster = paste0("cluster_", cluster))

# ścieżka bazowa
clusters_df <- risperidone_st_data_half$clusters %>%
  select(sample, barcode, cluster = cluster_resolution_0.4) %>%
  mutate(cluster = paste0("cluster_", cluster))

# ścieżka bazowa
base_dir <- "data/risperidone/ifpas-jackcraig-spatialdata"

for (sample_id in unique(barcode_risperidone_half$sample)) {
  barcode_risperidone_half %>%
    filter(sample == sample_id) %>%
    left_join(clusters_df, by = c("sample", "barcode")) %>%
    filter(!is.na(cluster)) %>%  # <-- usuwa barcody bez przypisanego klastra
    rename(
      sample = "sample_id",
      row = "row_index",
      col = "col_index",
      imagerow = "imagerow_pixels",
      imagecol = "imagecol_pixels"
    ) %>%
    dplyr::select(sample_id, barcode, tissue, row_index, col_index, imagerow_pixels, imagecol_pixels, cluster) %>%
    write.table(
      file = file.path(base_dir, sample_id, paste0(sample_id, "_barcode_coordinates.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
}

rm(sample_id, base_dir, clusters_df)
# save expressiont data
risperidone_st_data_half$raw_data$data %>% .[1:10, 1:10] %>% 
  as.data.frame()



# library(Matrix)

# Extract data matrix
expr_mat <- risperidone_st_data_half$raw_data$data

# Ścieżka bazowa
base_dir <- "data/risperidone/ifpas-jackcraig-spatialdata"

# Wyciągnij sample_id z kolumn (np. "S6230Nr3_AAACAAG..." → "S6230Nr3")
column_names <- colnames(expr_mat)
sample_ids <- unique(sub("_.*", "", column_names))

# Pętla po unikalnych próbkach
for (sample_id in sample_ids) {
  # Wybierz kolumny dla danego sample_id
  cols <- grep(paste0("^", sample_id, "_"), column_names, value = TRUE)
  
  # Podmacierz z tymi kolumnami
  sample_expr <- expr_mat[, cols, drop = FALSE]
  
  # Konwersja na data.frame z peak_id jako kolumna
  sample_expr_df <- as.data.frame(as.matrix(sample_expr))
  sample_expr_df <- tibble::rownames_to_column(sample_expr_df, "peak_id")
  
  # Usuń prefiks sample_id z nazw kolumn (czyli z kolumn typu "S6230Nr3_AAACAAG..." → "AAACAAG...")
  colnames(sample_expr_df)[-1] <- sub(paste0("^", sample_id, "_"), "", colnames(sample_expr_df)[-1])
  
  # Ścieżka zapisu
  output_file <- file.path(base_dir, sample_id, paste0(sample_id, "_expression_matrix.tsv"))
  
  # Zapisz do pliku TSV
  write.table(
    sample_expr_df,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

rm(expr_mat, column_names, sample_ids, sample_id, cols, sample_expr, sample_expr_df, base_dir)

#
risperidone_st_data_half$clusters %>% head %>% 
  select(c(sample, barcode, cluster_resolution_0.4))

  