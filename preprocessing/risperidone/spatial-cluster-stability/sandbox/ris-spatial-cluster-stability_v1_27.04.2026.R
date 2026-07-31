load("results/risperidone/risperidone-half.RData")


cluster_df <- risperidone_st_data_half$clusters %>%
  dplyr::mutate(spot_id = paste(sample, barcode, sep = "_"))

cluster_matrix <- cluster_df %>%
  dplyr::select(
    spot_id,
    dplyr::starts_with("cluster_resolution_")
  ) %>%
  tibble::column_to_rownames("spot_id")

risperidone_pca_embeddings_1_30 <- Embeddings(
  risperidone_integrate_half,
  "pca"
)[, 1:30]

common_spots <- intersect(
  rownames(risperidone_pca_embeddings_1_30),
  rownames(cluster_matrix)
)

feature_matrix <- risperidone_pca_embeddings_1_30[common_spots, ]
cluster_matrix <- cluster_matrix[common_spots, ]


evaluate_cluster_metrics_single <- function(
  feature_matrix,
  cluster_matrix,
  cluster_col = "cluster_resolution_0.4",
  dims = NULL
) {
  
  # 1) Walidacja
  if (!cluster_col %in% colnames(cluster_matrix)) {
    stop("Kolumna nie istnieje w cluster_matrix: ", cluster_col)
  }
  
  # opcjonalne przycięcie dims (jeśli chcesz np. tylko 1:30)
  if (!is.null(dims)) {
    feature_matrix <- feature_matrix[, dims, drop = FALSE]
  }
  
  # 2) Dopasowanie wierszy
  common_ids <- intersect(rownames(feature_matrix), rownames(cluster_matrix))
  if (length(common_ids) == 0) {
    stop("Brak wspólnych ID między feature_matrix i cluster_matrix")
  }
  feature_matrix <- feature_matrix[common_ids, , drop = FALSE]
  clusters <- cluster_matrix[common_ids, cluster_col]
  
  message("[1/4] Dane dopasowane: ", length(common_ids), " obserwacji")
  
  # 3) Liczenie macierzy odległości
  message("[2/4] Liczę macierz odległości (to może chwilę potrwać)...")
  dist_mat <- dist(feature_matrix)
  
  # 4) Silhouette
  message("[3/4] Liczę silhouette dla: ", cluster_col)
  sil <- cluster::silhouette(
    as.numeric(as.factor(clusters)),
    dist_mat
  )
  
  # 5) Podsumowanie
  res <- data.frame(
    cluster_col = cluster_col,
    resolution = as.numeric(gsub("cluster_resolution_", "", cluster_col)),
    n_clusters = length(unique(clusters)),
    silhouette_score = mean(sil[, 3], na.rm = TRUE)
  )
  
  message("[4/4] Gotowe. Silhouette = ", round(res$silhouette_score, 4),
          " | n_clusters = ", res$n_clusters)
  
  return(res)
}

cluster_cols <- c(
  "cluster_resolution_0.1",
  "cluster_resolution_0.2",
  "cluster_resolution_0.4",
  "cluster_resolution_0.5",
  "cluster_resolution_0.8"
)






feature_matrix <- risperidone_pca_embeddings_1_30[common_ids, , drop = FALSE]
cluster_matrix_matched <- cluster_matrix[common_ids, , drop = FALSE]

message("Dane dopasowane: ", length(common_ids), " obserwacji")
message("Liczę macierz odległości...")
dist_mat <- dist(feature_matrix)

cluster_metrics_list <- list()

# for (cluster_col in cluster_cols) {
#   
#   message("Liczę silhouette dla: ", cluster_col)
#   
#   clusters <- cluster_matrix_matched[[cluster_col]]
#   
#   sil <- cluster::silhouette(
#     as.numeric(as.factor(clusters)),
#     dist_mat
#   )
#   
#   cluster_metrics_list[[cluster_col]] <- data.frame(
#     cluster_col = cluster_col,
#     resolution = as.numeric(gsub("cluster_resolution_", "", cluster_col)),
#     n_clusters = length(unique(clusters)),
#     silhouette_score = mean(sil[, 3], na.rm = TRUE)
#   )
# }

risperidone_cluster_metrics <- dplyr::bind_rows(cluster_metrics_list)
