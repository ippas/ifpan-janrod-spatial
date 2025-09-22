plot_gene_heatmap_from_bundle <- function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("treatment","mouse_genotype"),
  group_palette = NULL,
  palette = NULL,
  zscore_colors = c("blue","white","red"),
  cluster_order = NULL,
  scale_by_row = TRUE,
  scale_limits = c(-2, 2),
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  divide_by_cluster = TRUE,
  divide_columns_by_groups = TRUE,
  column_gap = grid::unit(4, "mm"),
  row_gap = grid::unit(3, "mm"),
  row_title_rot = 0,
  row_dend_width = grid::unit(15, "mm"),
  draw_border = TRUE,
  border_col  = "black",
  border_lwd  = 1,
  row_names_simple = TRUE
){
  #' Plot a cluster-split heatmap from a bundle object
  #'
  #' This function extracts peak-level expression values from a `bundle` object (containing
  #' cluster-wise matrices), optionally applies log2 transformation and row-wise z-score
  #' normalization, and visualizes results using `ComplexHeatmap::Heatmap`. Rows can be
  #' split by cluster, columns by metadata-defined groups, and the heatmap can include
  #' colored annotations, custom palettes, and optional borders.
  #'
  #' @param bundle List containing:
  #'   * `metadata`: data.frame with at least columns matching `group_cols`
  #'   * `raw_counts_by_cluster` or `qn_counts_by_cluster`: named list of matrices,
  #'     one per cluster (rownames = peaks, colnames = samples)
  #' @param peaks_df Data.frame containing peaks to be visualized. Must include
  #'   columns specified in `peak_col` and `cluster_col`.
  #' @param peak_col Column name with peak identifiers (default = `"peak_id"`).
  #' @param cluster_col Column name with cluster identifiers (default = `"cluster"`).
  #' @param data_type Either `"raw"` or `"qn"`, selects which matrix list from `bundle` to use.
  #' @param log2_transform Logical, apply log2 transform to values (default TRUE).
  #' @param pseudocount Numeric, value added before log2 transform (default = 1).
  #' @param group_cols Character vector, metadata columns to annotate columns by (default: treatment + genotype).
  #' @param group_palette Named list with colors for `group_cols` values.
  #' @param palette Alias for `group_palette` (kept for backward compatibility).
  #' @param zscore_colors Vector of 3 colors used for z-score scale: low, mid, high.
  #' @param cluster_order Character vector defining row-split order. If NULL, tries to sort numerically.
  #' @param scale_by_row Logical, if TRUE z-score scale each row before plotting.
  #' @param scale_limits Numeric length-2 vector with min and max for clamping z-scores.
  #' @param show_colnames Logical, whether to show column names.
  #' @param cluster_rows Logical, whether to perform hierarchical clustering of rows (within each split).
  #' @param cluster_columns Logical, whether to cluster columns.
  #' @param divide_by_cluster Logical, whether to split rows by `cluster`.
  #' @param divide_columns_by_groups Logical, whether to split columns by `group_cols`.
  #' @param column_gap grid::unit specifying gap between column groups (default 4mm).
  #' @param row_gap grid::unit specifying gap between row-split blocks (default 3mm).
  #' @param row_title_rot Numeric rotation for row-split titles (0 = horizontal).
  #' @param row_dend_width grid::unit, width of the row dendrogram.
  #' @param draw_border Logical, whether to draw borders around heatmap blocks.
  #' @param border_col Character, color of block border (default black).
  #' @param border_lwd Numeric, line width of border (default 1).
  #' @param row_names_simple Logical, if TRUE shows only gene symbols (last part after `"-"`).
  #'   If FALSE, shows full `cluster|peak_id` identifiers.
  #'
  #' @return Invisibly returns a ComplexHeatmap object.
  #' @examples
  #' plot_gene_heatmap_from_bundle(
  #'   bundle = my_bundle,
  #'   peaks_df = selected_peaks,
  #'   peak_col = "name_id",
  #'   cluster_col = "cluster",
  #'   data_type = "qn",
  #'   row_names_simple = TRUE
  #' )
  
  stopifnot(is.list(bundle), "metadata" %in% names(bundle))
  data_type <- match.arg(data_type)
  
  message("📦 Using matrix type: ", data_type)
  mat_list_name <- switch(data_type,
                          raw = "raw_counts_by_cluster",
                          qn  = "qn_counts_by_cluster")
  if (!(mat_list_name %in% names(bundle))) {
    stop(sprintf("Matrix '%s' not found in bundle", mat_list_name))
  }
  if (length(zscore_colors) != 3) {
    stop("`zscore_colors` must have exactly 3 colors: c(low, mid, high).")
  }
  
  # ---- Extract values for selected peaks ----
  message("📥 Collecting values for ", nrow(peaks_df), " peaks from clusters...")
  df_list <- lapply(seq_len(nrow(peaks_df)), function(i){
    clust <- peaks_df[[cluster_col]][i]
    pid   <- peaks_df[[peak_col]][i]
    if (!(clust %in% names(bundle[[mat_list_name]]))) {
      stop(sprintf("Cluster '%s' not found in bundle", clust))
    }
    mat <- bundle[[mat_list_name]][[clust]]
    if (!(pid %in% rownames(mat))) {
      stop(sprintf("Peak_id '%s' not found in cluster '%s'", pid, clust))
    }
    vals <- as.data.frame(mat[pid, , drop = FALSE])
    vals$peak_id <- pid
    vals$cluster <- clust
    vals
  })
  
  # ---- Long-format data ----
  df_long <- dplyr::bind_rows(df_list, .id = "row_id") %>%
    tidyr::pivot_longer(
      cols = -c(row_id, peak_id, cluster),
      names_to = "sample_ID",
      values_to = "value"
    ) %>%
    dplyr::left_join(bundle$metadata, by = "sample_ID")
  
  # ---- Log2 transform ----
  if (isTRUE(log2_transform)) {
    message("🔢 Applying log2 transformation with pseudocount = ", pseudocount)
    df_long <- dplyr::mutate(df_long, value = log2(value + pseudocount))
  }
  
  # ---- Column grouping ----
  if (length(group_cols) > 0) {
    df_long <- dplyr::mutate(
      df_long,
      group = interaction(dplyr::across(dplyr::all_of(group_cols)),
                          sep = "_", drop = TRUE)
    )
  } else {
    df_long$group <- "all"
  }
  
  # ---- Cluster order ----
  if (is.null(cluster_order)) {
    cluster_levels <- unique(peaks_df[[cluster_col]])
    nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", cluster_levels)))
    cluster_order <- if (all(!is.na(nums))) cluster_levels[order(nums)] else cluster_levels
  }
  
  df_long$cluster <- factor(df_long$cluster, levels = cluster_order)
  df_long <- df_long %>%
    dplyr::mutate(peak_label = paste(cluster, peak_id, sep = "|")) %>%
    dplyr::arrange(cluster, peak_id)
  
  # ---- Build heatmap matrix ----
  heat_mat <- df_long %>%
    dplyr::select(peak_label, sample_ID, value) %>%
    dplyr::group_by(peak_label, sample_ID) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_ID, values_from = value) %>%
    tibble::column_to_rownames("peak_label") %>%
    as.matrix()
  
  mode(heat_mat) <- "numeric" # ensure numeric
  
  # ---- Row labels (mapper, if requested) ----
  if (isTRUE(row_names_simple)) {
    row_labels <- vapply(
      strsplit(gsub(".*\\|", "", rownames(heat_mat)), "-", fixed = TRUE),
      function(parts) parts[length(parts)],
      character(1)
    )
  } else {
    row_labels <- rownames(heat_mat)
  }
  
  # ---- Column annotation ----
  annot_col <- bundle$metadata %>%
    dplyr::mutate(
      treatment = factor(treatment, levels = c("saline", "risperidone")),
      mouse_genotype = factor(mouse_genotype, levels = c("wtwt", "wtdel"))
    ) %>%
    tibble::column_to_rownames("sample_ID")
  
  sample_order <- annot_col %>%
    dplyr::arrange(treatment, mouse_genotype) %>%
    rownames()
  
  common_samples <- intersect(sample_order, colnames(heat_mat))
  heat_mat  <- heat_mat[, common_samples, drop = FALSE]
  annot_col <- annot_col[common_samples, , drop = FALSE]
  
  # ---- Row-wise scaling ----
  if (scale_by_row) {
    message("📏 Scaling rows (z-score) and clipping to ", paste(scale_limits, collapse = " ... "))
    heat_mat <- t(scale(t(heat_mat)))
    heat_mat[heat_mat < scale_limits[1]] <- scale_limits[1]
    heat_mat[heat_mat > scale_limits[2]] <- scale_limits[2]
  }
  
  # ---- Color palettes ----
  if (!is.null(group_palette) && !is.null(palette)) {
    warning("Both `group_palette` and `palette` provided. Using `group_palette`.")
  }
  pal_final <- if (!is.null(group_palette)) group_palette else if (!is.null(palette)) palette else
    list(treatment = c("saline" = "#1f77b4", "risperidone" = "#ff7f0e"),
         mouse_genotype = c("wtwt" = "#2ca02c", "wtdel" = "#d62728"))
  
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- circlize::colorRamp2(
    breaks = c(scale_limits[1], 0, scale_limits[2]),
    colors = zscore_colors
  )
  
  ann_legend_param <- lapply(group_cols, function(x) list(nrow = 1))
  names(ann_legend_param) <- group_cols
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df  = annot_col[group_cols],
    col = pal_final,
    annotation_legend_param = ann_legend_param,
    annotation_height = grid::unit(4, "mm"),
    gp = grid::gpar(col = NA)
  )
  
  row_split <- if (isTRUE(divide_by_cluster)) gsub("\\|.*", "", rownames(heat_mat)) else NULL
  column_split <- if (isTRUE(divide_columns_by_groups) && length(group_cols) > 0)
    interaction(annot_col[, group_cols, drop = FALSE], sep = "_", drop = TRUE) else NULL
  
  message("🔥 Drawing heatmap with ComplexHeatmap...")
  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = "z-score",
    row_labels = row_labels,
    top_annotation = ha,
    col = col_fun,
    border = draw_border,
    border_gp = grid::gpar(col = border_col, lwd = border_lwd),
    heatmap_legend_param = list(direction = "horizontal"),
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,
    cluster_columns = cluster_columns,
    column_split = column_split,
    column_gap = column_gap,
    row_split = row_split,
    row_gap = if (isTRUE(divide_by_cluster)) row_gap else grid::unit(0, "mm"),
    row_dend_width = row_dend_width,
    row_title_rot = row_title_rot,
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_title_side = "left",
    show_column_names = show_colnames,
    show_row_names = TRUE,
    column_title = sprintf("Heatmap (%s counts)", data_type),
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
