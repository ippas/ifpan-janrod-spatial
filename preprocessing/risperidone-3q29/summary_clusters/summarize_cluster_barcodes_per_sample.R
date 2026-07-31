ris3q29_st_data$clusters %>%
  select(sample, barcode, mouse_genotype, treatment, cluster_resolution_0.4) %>%
  mutate(
    group_order = case_when(
      treatment == "saline" & mouse_genotype == "wtwt" ~ 1,
      treatment == "saline" & mouse_genotype != "wtwt" ~ 2,
      treatment != "saline" & mouse_genotype == "wtwt" ~ 3,
      treatment != "saline" & mouse_genotype != "wtwt" ~ 4
    ),
    sample_id = paste0(
      ifelse(treatment == "saline", "sal", "ris"),
      "_",
      ifelse(mouse_genotype == "wtwt", "wt", "del"),
      "_",
      sample
    ),
    cluster_id = paste0("cluster_", cluster_resolution_0.4),
    cluster_name = ifelse(
      cluster_resolution_0.4 == 0,
      "cluster_22",
      paste0("cluster_", cluster_resolution_0.4)
    )
  ) %>%
  dplyr::count(group_order, sample_id, cluster_id, cluster_name, name = "n_barcodes") %>%
  arrange(
    as.numeric(gsub("cluster_", "", cluster_id)),
    group_order,
    sample_id
  ) %>%
  select(-group_order) %>%
  tidyr::pivot_wider(
    id_cols = c(cluster_id, cluster_name),
    names_from = sample_id,
    values_from = n_barcodes,
    values_fill = 0
  ) %>%
  arrange(as.numeric(gsub("cluster_", "", cluster_id))) %>%
  writexl::write_xlsx(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/summary_clusters/ris3q29_cluster0.4_barcode_counts_per_sample_17.06.2026.xlsx"
  )
