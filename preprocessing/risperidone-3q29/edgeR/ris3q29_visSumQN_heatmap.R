# ##############################################################################
# ---- preapare data ----
# ##############################################################################

ris3q29_bundle_data_to_visualization$raw_counts_by_cluster$cluster_0 %>% head

read.delim(
  file = "data/risperidone-3q29/data_from_magda/3q29-ris-spatial-edger-gene-list-pInteraction0.01v2.csv",
  sep = ",",
  header = TRUE,
  row.names = NULL
) %>% 
  select(-X) %>% 
  mutate(name_id = paste0(peak_id, '-', gene_symbol)) %>% 
  select(c(name_id, everything())) -> df

df <- read_excel("data/risperidone-3q29/tmp/3q29-ris-spatial-edger-gene-list-pInteraction0.01-all.xlsx")


# ##############################################################################
# ---- function ----
# ##############################################################################
plot_gene_heatmap_from_bundle <- function(
  bundle,
  peaks_df,
  peak_col    = "peak_id",
  cluster_col = "cluster",
  data_type   = c("raw","qn"),
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("treatment","mouse_genotype"),
  palette = NULL,
  cluster_order = NULL,
  scale_by_row = TRUE,
  scale_limits = c(-2, 2),
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  divide_by_cluster = TRUE,
  row_gap = grid::unit(3, "mm"),
  row_title_rot = 0,
  row_dend_width = grid::unit(15, "mm")
){
  stopifnot(is.list(bundle), "metadata" %in% names(bundle))
  data_type <- match.arg(data_type)
  
  mat_list_name <- switch(
    data_type,
    raw = "raw_counts_by_cluster",
    qn  = "qn_counts_by_cluster"
  )
  if (!(mat_list_name %in% names(bundle))) {
    stop(sprintf("Brak macierzy '%s' w bundle", mat_list_name))
  }
  
  message("📥 Przygotowuję dane dla ", nrow(peaks_df), " peaków...")
  df_list <- lapply(seq_len(nrow(peaks_df)), function(i){
    clust <- peaks_df[[cluster_col]][i]
    pid   <- peaks_df[[peak_col]][i]
    
    if (!(clust %in% names(bundle[[mat_list_name]]))) {
      stop(sprintf("Cluster '%s' nie istnieje w bundle", clust))
    }
    mat <- bundle[[mat_list_name]][[clust]]
    if (!(pid %in% rownames(mat))) {
      stop(sprintf("Peak_id '%s' nie znaleziony w %s", pid, clust))
    }
    
    vals <- mat[pid, , drop = FALSE] %>% as.data.frame()
    vals$peak_id <- pid
    vals$cluster <- clust
    vals
  })
  
  df_long <- dplyr::bind_rows(df_list, .id = "row_id") %>%
    tidyr::pivot_longer(
      cols = -c(row_id, peak_id, cluster),
      names_to = "sample_ID",
      values_to = "value"
    ) %>%
    dplyr::left_join(bundle$metadata, by = "sample_ID")
  
  if (isTRUE(log2_transform)) {
    df_long <- dplyr::mutate(df_long, value = log2(value + pseudocount))
  }
  
  if (length(group_cols) > 0) {
    df_long <- dplyr::mutate(df_long,
                             group = interaction(dplyr::across(dplyr::all_of(group_cols)),
                                                 sep = "_", drop = TRUE))
  } else {
    df_long$group <- "all"
  }
  
  # 🔢 DOMYŚLNE SORTOWANIE KLASRÓW NUMERYCZNIE
  if (is.null(cluster_order)) {
    cluster_levels <- unique(peaks_df[[cluster_col]])
    cluster_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", cluster_levels)))
    if (all(!is.na(cluster_nums))) {
      cluster_order <- cluster_levels[order(cluster_nums)]
    } else {
      cluster_order <- cluster_levels
    }
  }
  
  df_long$cluster <- factor(df_long$cluster, levels = cluster_order)
  
  df_long <- df_long %>%
    dplyr::mutate(peak_label = paste(cluster, peak_id, sep = "|"))
  
  # 📏 Sortujemy wiersze wg cluster_order + peak_id
  df_long <- dplyr::arrange(df_long, cluster, peak_id)
  
  heat_mat <- df_long %>%
    dplyr::select(peak_label, sample_ID, value) %>%
    tidyr::pivot_wider(names_from = sample_ID, values_from = value) %>%
    tibble::column_to_rownames("peak_label") %>%
    as.matrix()
  
  # Upewniamy się, że wiersze są w żądanej kolejności
  rownames_order <- rownames(heat_mat)
  heat_mat <- heat_mat[rownames_order, , drop = FALSE]
  
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
  heat_mat <- heat_mat[, common_samples, drop = FALSE]
  annot_col <- annot_col[common_samples, , drop = FALSE]
  
  # 🔥 Skalowanie ręczne
  if (scale_by_row) {
    heat_mat <- t(scale(t(heat_mat)))  # z-score wierszami
    heat_mat[heat_mat < scale_limits[1]] <- scale_limits[1]
    heat_mat[heat_mat > scale_limits[2]] <- scale_limits[2]
  }
  
  if (is.null(palette)) {
    palette <- list(
      treatment = c("saline" = "#1f77b4", "risperidone" = "#ff7f0e"),
      mouse_genotype = c("wtwt" = "#2ca02c", "wtdel" = "#d62728")
    )
  }
  
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- circlize::colorRamp2(
    breaks = c(scale_limits[1], 0, scale_limits[2]),
    colors = c("blue", "white", "red")
  )
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = annot_col[group_cols],
    col = palette,
    annotation_height = grid::unit(4, "mm"),
    gp = grid::gpar(col = NA)
  )
  
  row_split <- NULL
  if (isTRUE(divide_by_cluster)) {
    row_split <- gsub("\\|.*", "", rownames(heat_mat))
    row_split <- factor(row_split, levels = cluster_order)
  }
  
  # 📌 Wyłączamy klastrowanie wierszy, jeśli cluster_order jest podane jawnie
  cluster_rows_final <- if (!is.null(cluster_order)) FALSE else cluster_rows
  
  message("🔥 Rysuję heatmapę z ComplexHeatmap...")
  ComplexHeatmap::Heatmap(
    heat_mat,
    name = "z-score",
    top_annotation = ha,
    col = col_fun,
    cluster_rows = cluster_rows_final,
    cluster_columns = cluster_columns,
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
}


# ##############################################################################
# ---- test function ----
# ##############################################################################
# przykładowa tabela z interesującymi peakami
peaks_df <- data.frame(
  peak_id = c("peak-305-Xkr4", "peak-575-Mrpl15", "peak-601-Lypla1", "peak-622-Tcea1"),
  cluster = c("cluster_0",     "cluster_0",      "cluster_3",       "cluster_5")
)

big_tbl %>% 
  filter(FDR_allSalAllRis < 0.01 | FDR_allWtAllDel < 0.01) -> df


plot_gene_heatmap_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = TRUE,
  pseudocount = 1,
  cluster_order = paste0("cluster_", 0:19),
  palette = list(                               # alias — OK
    treatment = c("saline" = "#66c2a5", "risperidone" = "#fc8d62"),
    mouse_genotype = c("wtwt" = "#8da0cb", "wtdel" = "#e78ac3")
  ),
  zscore_colors = c("navy","white","firebrick3"),
  scale_by_row = TRUE,
  scale_limits = c(-2, 2)
)



