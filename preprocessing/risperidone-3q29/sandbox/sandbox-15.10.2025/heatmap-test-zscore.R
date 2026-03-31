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
  scale_limits = c(-3, 3),
  row_names_simple = TRUE,
  reference_group = "wtwt_saline",   # 👈 nowy argument — kolumna odniesienia
  reference_mode  = "zscore"         # 👈 sposób przeliczenia względem referencji
)


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
  scale_limits = c(-3, 3),
  row_names_simple = TRUE,
  reference_group = "wtwt_saline",   # 👈 nowy argument — kolumna odniesienia
  reference_mode  = "zscore"         # 👈 sposób przeliczenia względem referencji
)
