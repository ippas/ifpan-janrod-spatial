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


ris3q29_bundle_data_to_visualization$raw_counts_by_cluster$cluster_0 %>% head

plot_gene_heatmap_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  return_matrix = T,
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
  scale_limits = c(-2, 2),
  row_names_simple = T
) -> tmp


plot_gene_heatmap_group_means_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = TRUE,
  pseudocount = 1,
  group_cols = c("treatment","mouse_genotype"),
  # group_order = c("wtwt_saline", "wtwt_risperidone", "wtdel_saline", "wtdel_risperidone"),
  # cluster_order = paste0("cluster_", 0:19),
  palette = list(
    mouse_genotype = c("wtwt" = "#8da0cb", "wtdel" = "#e78ac3"),
    treatment      = c("saline" = "#66c2a5", "risperidone" = "#fc8d62")
  ),
  zscore_colors = c("navy","white","firebrick3"),
  scale_by_row = TRUE,
  scale_limits = c(-2, 2),
  row_names_simple = TRUE
)

plot_gene_heatmap_group_deltas_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = F,
  pseudocount = 1,
  group_cols = c("mouse_genotype","treatment"),
  contrasts = list(
    delta_risperidone = c("wtwt_risperidone","wtdel_risperidone"),
    delta_saline      = c("wtwt_saline","wtdel_saline")
  ),
  # cluster_order = paste0("cluster_", 0:19),  # wymuszona kolejność
  scale_by_row = F,
  zscore_colors = c("navy","white","firebrick3"),
  scale_limits = c(-2, 2)
)

plot_gene_heatmap_group_deltas_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = F,
  pseudocount = 1,
  group_cols = c("mouse_genotype","treatment"),
  contrasts = list(
    delta_wt = c("wtwt_saline","wtwt_risperidone"),
    delta_del      = c("wtdel_saline","wtdel_risperidone")
  ),
  # cluster_order = paste0("cluster_", 0:19),  # wymuszona kolejność
  scale_by_row = F,
  zscore_colors = c("navy","white","firebrick3"),
  scale_limits = c(-2, 2)
)


plot_gene_heatmap_group_ratios_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  log2_transform = F,          # → legenda: "log2(qn/saline)"
  group_cols = c("mouse_genotype","treatment"),
  group_order = c("wtwt_saline","wtwt_risperidone","wtdel_saline","wtdel_risperidone"),
  scale_by_row = FALSE,           # bez z-score
  zscore_colors = c("navy","white","firebrick3"),
  palette = list(
    mouse_genotype = c("wtwt" = "#8da0cb", "wtdel" = "#e78ac3"),
    treatment      = c("saline" = "#66c2a5", "risperidone" = "#fc8d62")
  ),
  row_names_simple = TRUE
)
