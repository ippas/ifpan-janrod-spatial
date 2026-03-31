plot_gene_heatmap_group_means_refZscore_from_bundle <- function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("mouse_genotype","treatment"),
  group_order = NULL,
  group_palette = NULL,
  palette = NULL,
  zscore_colors = c("blue","white","red"),
  cluster_order = NULL,
  scale_by_row = TRUE,
  scale_limits = c(-2, 2),
  reference_group = "wtwt_saline",
  reference_mode  = c("zscore", "diff", "ratio"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_simple = TRUE,
  draw_border = TRUE,
  border_col  = "black",
  border_lwd  = 1,
  row_gap = grid::unit(3, "mm"),
  column_gap = grid::unit(4, "mm"),
  return_matrix = FALSE
){
  stopifnot(is.list(bundle), "metadata" %in% names(bundle))
  data_type <- match.arg(data_type)
  reference_mode <- match.arg(reference_mode)
  
  message("📦 Using matrix type: ", data_type)
  mat_list_name <- switch(data_type,
                          raw = "raw_counts_by_cluster",
                          qn  = "qn_counts_by_cluster")
  
  # ---- Extract & merge ----
  df_list <- lapply(seq_len(nrow(peaks_df)), function(i){
    clust <- peaks_df[[cluster_col]][i]
    pid   <- peaks_df[[peak_col]][i]
    if (!clust %in% names(bundle[[mat_list_name]])) {
      warning("⚠️ Cluster ", clust, " not found in bundle. Skipping...")
      return(NULL)
    }
    mat   <- bundle[[mat_list_name]][[clust]]
    if (!pid %in% rownames(mat)) {
      warning("⚠️ Peak ", pid, " not found in cluster ", clust, ". Skipping...")
      return(NULL)
    }
    vals  <- as.data.frame(mat[pid, , drop = FALSE])
    vals$peak_id <- pid
    vals$cluster <- clust
    vals
  })
  
  df_list <- Filter(Negate(is.null), df_list)
  if (length(df_list) == 0) stop("❌ Żaden z podanych peaków nie został znaleziony w klastrach!")
  
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
  
  if (return_matrix) {
    message("📤 Returning preprocessed group-mean matrix (no z-score).")
    return(heat_mat)
  }
  
  # ---- row labels ----
  row_labels <- if (row_names_simple) {
    vapply(
      strsplit(gsub(".*\\|", "", rownames(heat_mat)), "-", fixed = TRUE),
      function(parts) parts[length(parts)],
      character(1)
    )
  } else rownames(heat_mat)
  
  # ---- reference-based scaling ----
  if (scale_by_row) {
    if (!reference_group %in% colnames(heat_mat)) {
      stop("❌ Kolumna referencyjna '", reference_group, "' nie istnieje w macierzy!")
    }
    ref_values <- heat_mat[, reference_group]
    
    if (reference_mode == "zscore") {
      ref_mean <- ref_values
      ref_sd   <- apply(heat_mat, 1, sd, na.rm = TRUE)
      ref_sd[ref_sd == 0 | is.na(ref_sd)] <- 1
      heat_mat <- sweep(heat_mat, 1, ref_mean, "-")
      heat_mat <- sweep(heat_mat, 1, ref_sd, "/")
    } else if (reference_mode == "diff") {
      heat_mat <- sweep(heat_mat, 1, ref_values, "-")
    } else if (reference_mode == "ratio") {
      heat_mat <- sweep(heat_mat, 1, ref_values, "/")
    }
    
    heat_mat[heat_mat < scale_limits[1]] <- scale_limits[1]
    heat_mat[heat_mat > scale_limits[2]] <- scale_limits[2]
  }
  
  # ---- row split ----
  row_split <- gsub("\\|.*", "", rownames(heat_mat))
  if (is.null(cluster_order)) {
    row_split_num <- as.integer(gsub("cluster_", "", row_split))
    cluster_levels <- paste0("cluster_", sort(unique(row_split_num)))
  } else {
    cluster_levels <- cluster_order
  }
  row_split <- factor(row_split, levels = cluster_levels, ordered = TRUE)
  
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
  
  message("🔥 Drawing reference-based heatmap...")
  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = sprintf("ref-%s", reference_mode),
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
    column_title = sprintf("Group-mean Heatmap (ref: %s, %s counts)", reference_group, data_type),
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



plot_gene_heatmap_group_means_refZscore_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("mouse_genotype", "treatment"),
  group_order = c("wtwt_saline", "wtwt_risperidone", "wtdel_saline", "wtdel_risperidone"),
  # cluster_order = paste0("cluster_", 0:19),
  palette = list(
    mouse_genotype = c("wtwt" = "#8da0cb", "wtdel" = "#e78ac3"),
    treatment      = c("saline" = "#66c2a5", "risperidone" = "#fc8d62")
  ),
  zscore_colors = c("navy","white","firebrick3"),
  scale_by_row = TRUE,
  scale_limits = c(-4, 4),
  row_names_simple = TRUE,
  reference_group = "wtwt_saline",   # 👈 nowy argument — kolumna odniesienia
  reference_mode  = "zscore"         # 👈 sposób przeliczenia względem referencji
)
