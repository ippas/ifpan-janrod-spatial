
data_to_heatmap_test <- plot_gene_heatmap_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  peaks_df    = df,
  peak_col    = "name_id",
  cluster_col = "cluster",
  data_type   = "qn",
  group_cols  = c("mouse_genotype","treatment"),
  group_palette = list(
    treatment      = c("saline"="#66c2a5","risperidone"="#fc8d62"),
    mouse_genotype = c("wtwt"="#8da0cb","wtdel"="#e78ac3")
  ),
  cluster_order = paste0("cluster_", 0:19),
  zscore_colors = c("navy", "white", "firebrick3"),
  impute_by_group = T,
  add_marker_for_imputed_data = TRUE,
  log2_transform = F,
  
  return_matrix = T,
  scale_by_row = F
)


data_to_heatmap_test %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(rowname, into = c("cluster", "name_id"), sep = "\\|") %>% 
  separate(name_id, into = c("peak_id", "gene_symbol"), sep = "-(?=[^-]+$)", remove = FALSE) %>% 
  gather(key = "sample_ID", value = "value", -cluster, -name_id, -peak_id, -gene_symbol) %>% 
  left_join(., metadata_ris3q29[,c("sample_ID", "mouse_genotype", "treatment")]) %>% 
  mutate(group = paste0(mouse_genotype, "_", treatment)) %>% 
  group_by(cluster, name_id, peak_id, gene_symbol, group) %>%
  nest() %>% 
  mutate(mean_nameID_cluster = map(data, ~ .x %>% .$value %>% mean)) %>%
  mutate(sd_nameID_cluster = map(data, ~ .x %>% .$value %>% sd)) %>%
  unnest() %>% 
  ungroup() %>% 
  group_by(cluster, name_id, peak_id, gene_symbol) %>%
  nest() %>% 
  mutate(mean_nameID_cluster_wtwtSaline = map(data, ~ .x %>% filter(group == "wtwt_saline") %>% .$value %>% mean)) %>% 
  mutate(sd_nameID_cluster_wtwtSaline = map(data, ~ .x %>% filter(group == "wtwt_saline") %>% .$value %>% sd)) %>% 
  unnest() %>% 
  group_by(cluster, name_id, group) %>% 
  mutate(zscore = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(zscore_vs_saline = (value - mean_nameID_cluster_wtwtSaline) / sd_nameID_cluster_wtwtSaline) %>% 
  mutate(cluster_nameID = paste0(cluster, "_", name_id)) %>% 
  mutate(value_by_wtwtSaline = value/mean_nameID_cluster_wtwtSaline) %>% 
  group_by(cluster, name_id, group) %>% 
  mutate(zscore_normWtwtSaline = (value_by_wtwtSaline - mean(value_by_wtwtSaline, na.rm = TRUE)) / sd(value_by_wtwtSaline, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(column_ID = paste0(sample_ID, "_", group))



  # 


