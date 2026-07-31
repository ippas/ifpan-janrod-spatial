risperidone_pca_embeddings_1_30 <- Embeddings(
  risperidone_integrate_half,
  "pca"
)[, 1:30]

cluster_df <- risperidone_st_data_half$clusters %>%
  dplyr::mutate(
    spot_id = paste(sample, barcode, sep = "_")
  )

cluster_matrix <- cluster_df %>%
  dplyr::select(
    spot_id,
    dplyr::starts_with("cluster_resolution_")
  ) %>%
  tibble::column_to_rownames("spot_id")


# =========================
# 2. Dopasowanie PCA i klastrów
# =========================

common_ids <- intersect(
  rownames(risperidone_pca_embeddings_1_30),
  rownames(cluster_matrix)
)

feature_matrix <- risperidone_pca_embeddings_1_30[common_ids, , drop = FALSE]
cluster_matrix_matched <- cluster_matrix[common_ids, , drop = FALSE]

message("Dopasowano obserwacje: ", length(common_ids))

output_dir <- "data/risperidone/silhouette_analysis/"

# PCA embeddings
write.table(
  feature_matrix,
  file = file.path(output_dir, "risperidone_half_pca_embeddings_1_30.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# klastry
write.table(
  cluster_matrix_matched,
  file = file.path(output_dir, "risperidone_half_clusters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

risperidone_st_data_half$bcs_information %>% 
  dplyr::select(-c(height, width)) %>%  
  dplyr::mutate(barcode = paste0(sample, "_", barcode)) %>% 
  dplyr::select(-sample) %>% 
  dplyr::filter(barcode %in% common_ids) %>% 
  dplyr::select(barcode, imagerow, imagecol) %>% 
  write.table(
    file = file.path(output_dir, "risperidone_half_spatial_coords.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )  



