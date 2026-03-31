risperidone_integrate_half <- FindClusters(risperidone_integrate_half, resolution = 0.4)

markers_res0.4 <- FindAllMarkers(
  # object = clozapine_integrate_half,
  object = risperidone_integrate_half,
  assay  = "RNA",      # albo "RNA"
  slot   = "data",
  test.use = "wilcox",
  logfc.threshold = 0.5,
  min.pct = 0.25,
  only.pos = TRUE
)



## ===============================
## 1) Prepare data
## ===============================
markers_tbl <- markers_res0.4 %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  rownames_to_column(var = "rownames") %>% 
  select(-rownames) %>%
  rename(gene = "peak_id") %>% 
  select(peak_id, everything()) %>% 
  left_join(., risperidone_st_data_half$raw_data$annotate %>% 
              select(peak_id, gene_name),
            by = "peak_id")

## ===============================
## 2) Split into ordered list
## ===============================
cluster_levels <- markers_tbl %>%
  distinct(cluster) %>%
  mutate(cluster_id = as.integer(str_remove(cluster, "cluster_"))) %>%
  arrange(cluster_id) %>%
  pull(cluster)

markers_list <- markers_tbl %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  arrange(cluster) %>%
  group_split(cluster, .keep = TRUE)

names(markers_list) <- cluster_levels

## ===============================
## 3) Write to XLSX (1 sheet per cluster)
## ===============================
out_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone/markersByClustersRisperidone_res0.4_12.02.2026.xlsx"  # <- podmień na swoją ścieżkę jeśli chcesz

wb <- createWorkbook()

walk(names(markers_list), function(cl) {
  df <- markers_list[[cl]]
  
  addWorksheet(wb, cl)             # sheet name = "cluster_0", "cluster_1", ...
  writeData(wb, cl, df)
  
  ## freeze first row
  freezePane(wb, cl, firstRow = TRUE)
  
  ## (opcjonalnie) dopasuj szerokości kolumn
  # setColWidths(wb, cl, cols = 1:ncol(df), widths = "auto")
})

saveWorkbook(wb, out_file, overwrite = TRUE)



# ##############################################################################
# ---- clusters log2FC > 0.3 & min.pct > 30% ----
# ##############################################################################

markers_res0.4 <- FindAllMarkers(
  object = risperidone_integrate_half,
  assay  = "RNA",
  slot   = "data",
  test.use = "wilcox",
  logfc.threshold = 0.3,  # <- zmiana
  min.pct = 0.3,          # <- zmiana
  only.pos = TRUE
)

## ===============================
## 1) Prepare data
## ===============================
markers_tbl <- markers_res0.4 %>%
  mutate(cluster = paste0("cluster_", cluster)) %>%
  rownames_to_column(var = "rownames") %>%
  select(-rownames) %>%
  rename(gene = "peak_id") %>%
  select(peak_id, everything()) %>%
  left_join(
    risperidone_st_data_half$raw_data$annotate %>%
      select(peak_id, gene_name),
    by = "peak_id"
  )

## ===============================
## 2) Split into ordered list
## ===============================
cluster_levels <- markers_tbl %>%
  distinct(cluster) %>%
  mutate(cluster_id = as.integer(str_remove(cluster, "cluster_"))) %>%
  arrange(cluster_id) %>%
  pull(cluster)

markers_list <- markers_tbl %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  arrange(cluster) %>%
  group_split(cluster, .keep = TRUE)

names(markers_list) <- cluster_levels

## ===============================
## 3) Write to XLSX (1 sheet per cluster)
## ===============================
out_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone/markersByClustersRisperidone_res0.4_log2FC0.3minPct0.3_16.02.2026.xlsx"

wb <- createWorkbook()

walk(names(markers_list), function(cl) {
  df <- markers_list[[cl]]
  
  addWorksheet(wb, cl)
  writeData(wb, cl, df)
  
  ## freeze first row
  freezePane(wb, cl, firstRow = TRUE)
  
  ## (optional) auto-width
  # setColWidths(wb, cl, cols = 1:ncol(df), widths = "auto")
})

saveWorkbook(wb, out_file, overwrite = TRUE)

