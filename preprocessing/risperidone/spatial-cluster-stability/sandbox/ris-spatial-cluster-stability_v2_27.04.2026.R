cluster_cols <- colnames(cluster_matrix_matched)


cluster_metrics_list <- vector("list", length(cluster_cols))
names(cluster_metrics_list) <- cluster_cols

for (i in seq_along(cluster_cols)) {
  
  cluster_col <- cluster_cols[i]
  
  message(
    "[", i, "/", length(cluster_cols), "] Liczę silhouette dla: ",
    cluster_col
  )
  
  clusters <- cluster_matrix_matched[[cluster_col]]
  
  sil <- cluster::silhouette(
    as.numeric(as.factor(clusters)),
    dist_mat
  )
  
  silhouette_score <- mean(sil[, 3], na.rm = TRUE)
  
  cluster_metrics_list[[i]] <- data.frame(
    cluster_col = cluster_col,
    resolution = as.numeric(gsub("cluster_resolution_", "", cluster_col)),
    n_clusters = length(unique(clusters)),
    silhouette_score = silhouette_score
  )
  
  rm(sil, clusters, silhouette_score)
  gc(verbose = FALSE)
}

risperidone_cluster_metrics <- dplyr::bind_rows(cluster_metrics_list)

risperidone_cluster_metrics
