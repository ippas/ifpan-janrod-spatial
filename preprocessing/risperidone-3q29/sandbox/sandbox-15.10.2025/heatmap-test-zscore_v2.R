# --- 1️⃣ Przygotowanie danych i obliczenie klasycznego z-score
df <- data_to_heatmap_test %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(rowname, into = c("cluster", "name_id"), sep = "\\|") %>% 
  separate(name_id, into = c("peak_id", "gene_symbol"),
           sep = "-(?=[^-]+$)", remove = FALSE) %>% 
  gather(key = "sample_ID", value = "value",
         -cluster, -name_id, -peak_id, -gene_symbol) %>% 
  left_join(metadata_ris3q29[, c("sample_ID", "mouse_genotype", "treatment")],
            by = "sample_ID") %>% 
  mutate(group = paste0(mouse_genotype, "_", treatment)) %>%
  
  # --- klasyczny z-score per (cluster + name_id)
  group_by(cluster, name_id, peak_id, gene_symbol) %>%
  mutate(
    mean_nameID_cluster = mean(value, na.rm = TRUE),
    sd_nameID_cluster   = sd(value, na.rm = TRUE),
    zscore = (value - mean_nameID_cluster) / sd_nameID_cluster
  ) %>%
  ungroup() %>%
  
  # --- dodatkowe kolumny identyfikacyjne
  mutate(
    cluster_nameID = paste0(cluster, "_", name_id),
    column_ID = paste0(sample_ID, "_", group)
  )

# --- 2️⃣ Przygotowanie macierzy do heatmapy
mat <- df %>%
  select(cluster_nameID, column_ID, zscore) %>%
  pivot_wider(names_from = column_ID, values_from = zscore) %>%
  column_to_rownames("cluster_nameID") %>%
  as.matrix()

mat[is.na(mat)] <- 0  # uzupełnij brakujące wartości

# --- 3️⃣ Ustalona kolejność grup (wtwt_saline po lewej)
group_order <- c("wtwt_saline", "wtwt_risperidone", "wtdel_saline", "wtdel_risperidone")

groups <- gsub(".*_(wtwt_saline|wtwt_risperidone|wtdel_saline|wtdel_risperidone)$",
               "\\1", colnames(mat))
group_colors <- c(
  "wtwt_saline" = "#66c2a5",
  "wtwt_risperidone" = "#fc8d62",
  "wtdel_saline" = "#8da0cb",
  "wtdel_risperidone" = "#e78ac3"
)

# --- posortuj kolumny wg grup (kontrola po lewej)
col_order <- order(factor(groups, levels = group_order))
mat <- mat[, col_order]
groups <- groups[col_order]

# --- adnotacje nad kolumnami
ha_col <- HeatmapAnnotation(
  Group = groups,
  col = list(Group = group_colors),
  annotation_height = unit(4, "mm")
)

# --- 4️⃣ Kolory i klastrowanie kolumn
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# --- wykonaj klastrowanie kolumn
hc_cols <- hclust(dist(t(mat)), method = "ward.D2")

# --- zamień dwie główne gałęzie stronami (dla efektu wizualnego)
hc_cols$merge[1, ] <- hc_cols$merge[1, c(2, 1)]

# --- 5️⃣ Rysowanie heatmapy
ht <- Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = as.dendrogram(hc_cols),  # 👈 używamy odwróconego drzewa
  show_row_names = TRUE,           
  show_column_names = TRUE,        
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12, rot = 90),
  top_annotation = ha_col,
  column_dend_height = unit(30, "mm"),
  row_dend_width = unit(30, "mm"),
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  column_title = "Sample ID and Group",
  row_title = "Cluster | Name ID",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# --- 6️⃣ Rysuj z legendą na dole
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
