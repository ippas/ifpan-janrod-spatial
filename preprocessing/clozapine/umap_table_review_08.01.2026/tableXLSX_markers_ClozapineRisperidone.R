load("results/clozapine/clozapine-half.RData")


clozapine_integrate_half <- FindClusters(clozapine_integrate_half, resolution = 0.4)

markers_res0.4 <- FindAllMarkers(
  object = clozapine_integrate_half,
  assay  = "RNA",      # albo "RNA"
  slot   = "data",
  test.use = "wilcox",
  logfc.threshold = 0.5,
  min.pct = 0.25,
  only.pos = TRUE
)

clozapine_st_data_half$raw_data$annotate %>% 
  select(peak_id, gene_name)



## ===============================
## 1) Prepare data
## ===============================
markers_tbl <- markers_res0.4 %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  rownames_to_column(var = "rownames") %>% 
  select(-rownames) %>%
  rename(gene = "peak_id") %>% 
  select(peak_id, everything()) %>% 
  left_join(., clozapine_st_data_half$raw_data$annotate %>% 
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
out_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/markersByClustersClozapine_res0.4_09.01.2025.xlsx"  # <- podmień na swoją ścieżkę jeśli chcesz

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
# ---- res = 0.45 ----
# ##############################################################################
clozapine_integrate_half <- FindClusters(clozapine_integrate_half, resolution = 0.45)

markers_res0.45 <- FindAllMarkers(
  object = clozapine_integrate_half,
  # assay  = "SCT",      # albo "RNA"
  slot   = "data",
  test.use = "wilcox",
  logfc.threshold = 0.5,
  min.pct = 0.25,
  only.pos = TRUE
)



## ===============================
## 1) Prepare data
## ===============================
markers_tbl <- markers_res0.45 %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  rownames_to_column(var = "rownames") %>% 
  select(-rownames) %>%
  rename(gene = "peak_id") %>% 
  select(peak_id, everything()) %>% 
  left_join(., clozapine_st_data_half$raw_data$annotate %>% 
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
out_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/markersByClustersClozapine_res0.45_09.01.2025.xlsx"  # <- podmień na swoją ścieżkę jeśli chcesz

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



tmp_env$risperidone_integrate_half


# ##############################################################################
# ---- risperidone ----
# ##############################################################################

## ============================================================
## 1) clusters (res 0.4) – DIRECTLY
## ============================================================
tmp_env$risperidone_integrate_half <- FindClusters(
  tmp_env$risperidone_integrate_half,
  resolution = 0.4
)

## (opcjonalnie, jeśli chcesz mieć pewność, że neighbors istnieją / są aktualne)
# tmp_env$risperidone_integrate_half <- FindNeighbors(tmp_env$risperidone_integrate_half, dims = 1:30)
# tmp_env$risperidone_integrate_half <- FindClusters(tmp_env$risperidone_integrate_half, resolution = 0.4)

## ============================================================
## 2) markers
## ============================================================
markers_res0.4_ris <- FindAllMarkers(
  object = tmp_env$risperidone_integrate_half,
  slot   = "data",
  test.use = "wilcox",
  logfc.threshold = 0.5,
  min.pct = 0.25,
  only.pos = TRUE
)

## ============================================================
## 3) annotations (peak_id -> gene_name)
## ============================================================
annot_ris <- tmp_env$risperidone_st_data_half$raw_data$annotate %>%
  select(peak_id, gene_name) %>%
  distinct(peak_id, .keep_all = TRUE)

## ============================================================
## 4) Prepare table
## ============================================================
markers_tbl_ris <- markers_res0.4_ris %>%
  mutate(cluster = paste0("cluster_", cluster)) %>%
  rownames_to_column(var = "rownames") %>% 
  select(-rownames) %>%
  rename(gene = "peak_id") %>% 
  select(peak_id, everything()) %>% 
  left_join(annot_ris, by = "peak_id")

## ============================================================
## 5) Split into ordered list
## ============================================================
cluster_levels_ris <- markers_tbl_ris %>%
  distinct(cluster) %>%
  mutate(cluster_id = as.integer(str_remove(cluster, "cluster_"))) %>%
  arrange(cluster_id) %>%
  pull(cluster)

markers_list_ris <- markers_tbl_ris %>%
  mutate(cluster = factor(cluster, levels = cluster_levels_ris)) %>%
  arrange(cluster) %>%
  group_split(cluster, .keep = TRUE)

names(markers_list_ris) <- cluster_levels_ris

## ============================================================
## 6) Write to XLSX
## ============================================================
out_file_ris <- "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/markersByClustersRisperidone_res0.4_09.01.2025.xlsx"

wb <- createWorkbook()

walk(names(markers_list_ris), function(cl) {
  df <- markers_list_ris[[cl]]
  addWorksheet(wb, cl)
  writeData(wb, cl, df)
  freezePane(wb, cl, firstRow = TRUE)   # freeze header row
})

saveWorkbook(wb, out_file_ris, overwrite = TRUE)

