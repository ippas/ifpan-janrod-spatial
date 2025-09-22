plot_gene_heatmap_group_means_from_bundle <-function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("mouse_genotype","treatment"),
  group_order = NULL,  # 👈 NOWY argument
  group_palette = NULL,
  palette = NULL,
  zscore_colors = c("blue","white","red"),
  cluster_order = NULL,
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
    pid   <- peaks_df[[peak_col]][i]
    mat   <- bundle[[mat_list_name]][[clust]]
    vals  <- as.data.frame(mat[pid, , drop = FALSE])
    vals$peak_id <- pid
    vals$cluster <- clust
    vals
  })
  
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
  
  # ---- detect available group levels ----
  group_levels <- df_long %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(group = interaction(dplyr::across(dplyr::all_of(group_cols)),
                                      sep = "_", drop = TRUE)) %>%
    dplyr::pull(group)
  
  if (is.null(group_order)) {
    group_order <- group_levels
  } else {
    missing_groups <- setdiff(group_order, group_levels)
    if (length(missing_groups) > 0) {
      warning("⚠️ Grupy z group_order nie znalezione w danych: ",
              paste(missing_groups, collapse = ", "))
    }
  }
  
  # ---- aggregate ----
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
  
  # ---- force column order ----
  available_cols <- colnames(heat_mat)
  common_cols <- intersect(group_order, available_cols)
  if (length(common_cols) == 0) {
    stop("❌ Żadna z kolumn z group_order nie pasuje do danych!")
  }
  heat_mat <- heat_mat[, common_cols, drop = FALSE]
  
  # ---- drop empty rows ----
  heat_mat <- heat_mat[rowSums(!is.na(heat_mat)) > 0, , drop = FALSE]
  mode(heat_mat) <- "numeric"
  
  # ---- row labels ----
  row_labels <- if (row_names_simple) {
    vapply(
      strsplit(gsub(".*\\|", "", rownames(heat_mat)), "-", fixed = TRUE),
      function(parts) parts[length(parts)],
      character(1)
    )
  } else rownames(heat_mat)
  
  # ---- row-wise zscore ----
  if (scale_by_row) {
    heat_mat <- t(scale(t(heat_mat)))
    heat_mat[heat_mat < scale_limits[1]] <- scale_limits[1]
    heat_mat[heat_mat > scale_limits[2]] <- scale_limits[2]
  }
  
  # ---- color palettes ----
  pal_final <- if (!is.null(group_palette)) group_palette else if (!is.null(palette)) palette else
    list(
      treatment = c("saline" = "#1f77b4", "risperidone" = "#ff7f0e"),
      mouse_genotype = c("wtwt" = "#2ca02c", "wtdel" = "#d62728")
    )
  
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- circlize::colorRamp2(
    breaks = c(scale_limits[1], 0, scale_limits[2]),
    colors = zscore_colors
  )
  
  # ---- column annotation ----
  annot_df <- data.frame(group = colnames(heat_mat)) %>%
    tidyr::separate(group, into = group_cols, sep = "_", remove = FALSE)
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df  = annot_df[group_cols],
    col = pal_final,
    annotation_height = grid::unit(4, "mm"),
    gp = grid::gpar(col = NA)
  )
  
  # ---- row split by cluster ----
  row_split <- gsub("\\|.*", "", rownames(heat_mat))
  
  message("🔥 Drawing group-mean heatmap...")
  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = "z-score",
    row_labels = row_labels,
    top_annotation = ha,
    col = col_fun,
    border = draw_border,
    border_gp = grid::gpar(col = border_col, lwd = border_lwd),
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,
    cluster_columns = cluster_columns,
    column_gap = column_gap,
    row_split = row_split,
    row_gap = row_gap,
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontface = "bold"),
    show_column_names = TRUE,
    column_title = sprintf("Group-mean Heatmap (%s counts)", data_type),
    row_names_gp = grid::gpar(fontsize = 8)
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE
  )
  
  invisible(ht)
}
