plot_gene_heatmap_group_deltas_from_bundle <- function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("mouse_genotype","treatment"),
  contrasts = list(
    delta_risperidone = c("wtwt_risperidone","wtdel_risperidone"),
    delta_saline      = c("wtwt_saline","wtdel_saline")
  ),
  cluster_order = NULL,
  zscore_colors = c("blue","white","red"),
  scale_by_row = TRUE,
  scale_limits = c(-2, 2),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_simple = TRUE,
  draw_border = TRUE,
  border_col  = "black",
  border_lwd  = 1,
  row_gap = grid::unit(3, "mm"),
  column_gap = grid::unit(4, "mm")
){
  stopifnot(is.list(bundle), "metadata" %in% names(bundle))
  data_type <- match.arg(data_type)
  
  message("📦 Using matrix type: ", data_type)
  mat_list_name <- switch(data_type,
                          raw = "raw_counts_by_cluster",
                          qn  = "qn_counts_by_cluster")
  
  # ---- Extract & merge ----
  df_list <- lapply(seq_len(nrow(peaks_df)), function(i){
    clust <- peaks_df[[cluster_col]][i]
    pid_orig <- peaks_df[[peak_col]][i]
    
    if (!clust %in% names(bundle[[mat_list_name]])) {
      warning("⚠️ Cluster ", clust, " not found in bundle. Skipping...")
      return(NULL)
    }
    
    mat <- bundle[[mat_list_name]][[clust]]
    pid <- pid_orig
    
    # dopasowanie po końcówce, jeśli nie znaleziono dokładnego
    if (!pid %in% rownames(mat)) {
      matches <- grep(paste0("-", pid, "$"), rownames(mat), value = TRUE)
      if (length(matches) == 1) {
        pid <- matches
      } else if (length(matches) > 1) {
        warning(sprintf(
          "⚠️ Peak '%s' matched multiple rows in cluster '%s'. Using first match: %s",
          pid_orig, clust, matches[1]
        ))
        pid <- matches[1]
      } else {
        if (i == 1) {
          warning(sprintf(
            "⚠️ Peak '%s' not found in cluster '%s'. Example available rownames: %s",
            pid_orig, clust, paste(head(rownames(mat), 5), collapse=", ")
          ))
        }
        return(NULL)
      }
    }
    
    vals  <- as.data.frame(mat[pid, , drop = FALSE])
    vals$peak_id <- pid_orig
    vals$cluster <- clust
    vals
  })
  
  df_list <- Filter(Negate(is.null), df_list)
  if (length(df_list) == 0) stop("❌ None of the peaks were found in any cluster!")
  
  df_long <- dplyr::bind_rows(df_list) %>%
    tidyr::pivot_longer(
      cols = -c(peak_id, cluster),
      names_to = "sample_ID",
      values_to = "value"
    ) %>%
    dplyr::left_join(bundle$metadata, by = "sample_ID")
  
  # ---- log2 transform ----
  if (isTRUE(log2_transform)) {
    df_long <- dplyr::mutate(df_long, value = log2(value + pseudocount))
  }
  
  # ---- group label ----
  df_long <- dplyr::mutate(
    df_long,
    group = interaction(dplyr::across(dplyr::all_of(group_cols)),
                        sep = "_", drop = TRUE)
  )
  
  # ---- aggregate by group ----
  df_agg <- df_long %>%
    dplyr::group_by(cluster, peak_id, group) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(peak_label = paste(cluster, peak_id, sep = "|"))
  
  # ---- build matrix ----
  heat_mat <- df_agg %>%
    dplyr::select(peak_label, group, value) %>%
    tidyr::pivot_wider(names_from = group, values_from = value) %>%
    tibble::column_to_rownames("peak_label") %>%
    as.matrix()
  
  available_groups <- colnames(heat_mat)
  
  # ---- compute contrasts ----
  delta_mat <- sapply(names(contrasts), function(delta_name){
    contrast <- contrasts[[delta_name]]
    if (!all(contrast %in% available_groups)) {
      stop("❌ Missing columns for contrast ", delta_name,
           ". Available: ", paste(available_groups, collapse=", "))
    }
    heat_mat[, contrast[1]] - heat_mat[, contrast[2]]
  })
  
  if (is.null(dim(delta_mat))) {
    delta_mat <- matrix(delta_mat, ncol = 1,
                        dimnames = list(rownames(heat_mat), names(contrasts)))
  }
  
  # ---- row labels ----
  row_labels <- if (row_names_simple) {
    vapply(
      strsplit(gsub(".*\\|", "", rownames(delta_mat)), "-", fixed = TRUE),
      function(parts) parts[length(parts)],
      character(1)
    )
  } else rownames(delta_mat)
  
  # ---- z-score or raw ----
  if (scale_by_row) {
    delta_mat <- t(scale(t(delta_mat)))
    delta_mat[delta_mat < scale_limits[1]] <- scale_limits[1]
    delta_mat[delta_mat > scale_limits[2]] <- scale_limits[2]
  } else {
    message("ℹ️ Skipping z-score and keeping raw delta values (no scale limits applied).")
  }
  
  # ---- enforce cluster order ----
  row_split <- gsub("\\|.*", "", rownames(delta_mat))
  if (is.null(cluster_order)) {
    cluster_levels <- paste0("cluster_", sort(unique(as.integer(gsub("cluster_", "", row_split)))))
  } else {
    cluster_levels <- cluster_order
  }
  row_split <- factor(row_split, levels = cluster_levels)
  
  # ---- color palette ----
  library(ComplexHeatmap)
  library(circlize)
  
  if (scale_by_row) {
    breaks <- c(scale_limits[1], 0, scale_limits[2])
  } else {
    val_min <- min(delta_mat, na.rm = TRUE)
    val_max <- max(delta_mat, na.rm = TRUE)
    breaks <- c(val_min, 0, val_max)
  }
  
  col_fun <- circlize::colorRamp2(
    breaks = breaks,
    colors = zscore_colors
  )
  
  message("🔥 Drawing delta heatmap with fixed cluster order...")
  ht <- ComplexHeatmap::Heatmap(
    delta_mat,
    name = if (scale_by_row) "delta (z-score)" else "delta (log2)",
    row_labels = row_labels,
    col = col_fun,
    border = draw_border,
    border_gp = grid::gpar(col = border_col, lwd = border_lwd),
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,           # ❗ wyłącza sortowanie bloków
    cluster_columns = cluster_columns,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    row_split = row_split,
    row_order = NULL,                     # używa kolejności factor levels
    row_gap = row_gap,
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontface = "bold"),
    show_column_names = TRUE,
    column_title = sprintf("Group-delta Heatmap (%s counts)", data_type),
    row_names_gp = grid::gpar(fontsize = 8)
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    merge_legends = TRUE
  )
  
  invisible(ht)
}
