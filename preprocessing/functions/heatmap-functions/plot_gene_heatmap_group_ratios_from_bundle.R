plot_gene_heatmap_group_ratios_from_bundle <- function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  ratio_eps = 1e-9,                    # używane tylko gdy log2_transform = FALSE
  baseline_treatment = "saline",       # referencja do dzielenia w obrębie genotypu
  group_cols = c("mouse_genotype","treatment"),
  group_order = c("wtwt_saline","wtwt_risperidone",
                  "wtdel_saline","wtdel_risperidone"),
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
  column_gap = grid::unit(4, "mm"),
  return_matrix = FALSE
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
    tidyr::pivot_longer(cols = -c(peak_id, cluster),
                        names_to = "sample_ID", values_to = "value") %>%
    dplyr::left_join(bundle$metadata, by = "sample_ID")
  
  # ---- log2 transform (opcjonalnie przed agregacją) ----
  if (isTRUE(log2_transform)) {
    df_long <- dplyr::mutate(df_long, value = log2(value + pseudocount))
  }
  
  # ---- group label ----
  df_long <- dplyr::mutate(
    df_long,
    group = interaction(dplyr::across(dplyr::all_of(group_cols)),
                        sep = "_", drop = TRUE)
  )
  
  # ---- aggregate (średnia per peak × cluster × grupa) ----
  df_agg <- df_long %>%
    dplyr::group_by(cluster, peak_id, group,
                    dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # ---- compute ratio/log2FC względem baseline_treatment w obrębie genotypu ----
  if (isTRUE(log2_transform)) {
    # mamy już log2(qn + pc) → log2FC = log2(x) - log2(baseline)
    df_ratio <- df_agg %>%
      dplyr::group_by(cluster, peak_id, mouse_genotype) %>%
      dplyr::mutate(
        baseline = value[treatment == baseline_treatment],
        value = value - baseline   # log2(qn/saline)
      ) %>%
      dplyr::ungroup()
    legend_name <- "log2(qn/saline)"
  } else {
    # bez log2: licz prawdziwy stosunek qn/saline (z eps, by uniknąć /0)
    df_ratio <- df_agg %>%
      dplyr::group_by(cluster, peak_id, mouse_genotype) %>%
      dplyr::mutate(
        baseline = value[treatment == baseline_treatment],
        value = (value + ratio_eps) / (baseline + ratio_eps)
      ) %>%
      dplyr::ungroup()
    legend_name <- "qn/saline"
  }
  
  # komunikat o brakach baseline
  if (any(!is.finite(df_ratio$baseline))) {
    warning("⚠️ Brak wartości baseline (", baseline_treatment,
            ") dla części wierszy — te wiersze będą pominięte.")
  }
  
  df_ratio <- df_ratio %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::mutate(peak_label = paste(cluster, peak_id, sep = "|"))
  
  # ---- build matrix ----
  heat_mat <- df_ratio %>%
    dplyr::select(peak_label, group, value) %>%
    tidyr::pivot_wider(names_from = group, values_from = value) %>%
    tibble::column_to_rownames("peak_label") %>%
    as.matrix()
  
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
    message("📤 Returning matrix: ", legend_name)
    return(heat_mat)
  }
  
  # ---- row labels ----
  row_labels <- if (row_names_simple) {
    vapply(strsplit(gsub(".*\\|", "", rownames(heat_mat)), "-", fixed = TRUE),
           function(parts) parts[length(parts)], character(1))
  } else rownames(heat_mat)
  
  # ---- optional row-wise z-score ----
  if (scale_by_row) {
    heat_mat <- t(scale(t(heat_mat)))
    heat_mat[heat_mat < scale_limits[1]] <- scale_limits[1]
    heat_mat[heat_mat > scale_limits[2]] <- scale_limits[2]
  }
  
  # ---- row split (kolejność klastrów) ----
  row_split <- gsub("\\|.*", "", rownames(heat_mat))
  if (is.null(cluster_order)) {
    row_split_num <- suppressWarnings(as.integer(gsub("cluster_", "", row_split)))
    cluster_levels <- paste0("cluster_", sort(unique(row_split_num)))
  } else {
    cluster_levels <- cluster_order
  }
  row_split <- factor(row_split, levels = cluster_levels, ordered = TRUE)
  
  # ---- palety do adnotacji kolumn ----
  pal_final <- if (!is.null(group_palette)) group_palette else if (!is.null(palette)) palette else
    list(
      treatment = c("saline" = "#1f77b4", "risperidone" = "#ff7f0e"),
      mouse_genotype = c("wtwt" = "#2ca02c", "wtdel" = "#d62728")
    )
  
  library(ComplexHeatmap)
  library(circlize)
  
  # ---- color function (zależna od trybu) ----
  if (scale_by_row) {
    # z-score → symetryczne limity wokół 0
    col_fun <- circlize::colorRamp2(
      breaks = c(scale_limits[1], 0, scale_limits[2]),
      colors = zscore_colors
    )
    legend_title <- "z-score"
  } else if (isTRUE(log2_transform)) {
    # log2(qn/saline) → centrum 0, zakres dynamiczny symetryczny
    max_abs <- max(abs(heat_mat), na.rm = TRUE)
    col_fun <- circlize::colorRamp2(
      breaks = c(-max_abs, 0, max_abs),
      colors = zscore_colors
    )
    legend_title <- legend_name  # "log2(qn/saline)"
  } else {
    # qn/saline → centrum 1, zakres dynamiczny (może być niesymetryczny)
    min_v <- min(heat_mat, na.rm = TRUE)
    max_v <- max(heat_mat, na.rm = TRUE)
    # zabezpieczenie, gdy wszystko == 1
    if (!is.finite(min_v) || !is.finite(max_v) || min_v == max_v) {
      min_v <- 0.5; max_v <- 1.5
    }
    col_fun <- circlize::colorRamp2(
      breaks = c(min_v, 1, max_v),
      colors = zscore_colors
    )
    legend_title <- legend_name  # "qn/saline"
  }
  
  # ---- column annotation ----
  annot_df <- data.frame(group = colnames(heat_mat)) %>%
    tidyr::separate(group, into = group_cols, sep = "_", remove = FALSE)
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = annot_df[group_cols],
    col = pal_final,
    annotation_height = grid::unit(4, "mm"),
    gp = grid::gpar(col = NA)
  )
  
  message("🔥 Drawing heatmap (", legend_title, ") ...")
  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = legend_title,
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
    column_title = sprintf("Group-wise %s Heatmap (%s counts)",
                           legend_title, data_type),
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
