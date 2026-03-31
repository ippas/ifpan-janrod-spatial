# ============================================================
# 🔥 Heatmapa z Z-score względem grupy kontrolnej wtwt_saline
# ============================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# --- 1️⃣ Przygotowanie danych i obliczenie z-score względem wtwt_saline
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
  
  # --- 🔹 Z-score względem średniej i SD z grupy kontrolnej (wtwt_saline)
  group_by(cluster, name_id, peak_id, gene_symbol) %>%
  mutate(
    mean_ctrl = mean(value[group == "wtwt_saline"], na.rm = TRUE),
    sd_ctrl   = sd(value[group == "wtwt_saline"], na.rm = TRUE),
    # zabezpieczenie przed dzieleniem przez 0
    sd_ctrl   = ifelse(is.na(sd_ctrl) | sd_ctrl == 0, 1, sd_ctrl),
    zscore    = (value - mean_ctrl) / sd_ctrl
  ) %>%
  ungroup() %>%
  
  # --- Dodatkowe kolumny identyfikacyjne
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

mat[is.na(mat)] <- 0  # uzupełnij brakujące wartości zerami

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
col_fun <- colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick3"))

# --- wykonaj klastrowanie kolumn
hc_cols <- hclust(dist(t(mat)), method = "ward.D2")

# --- zamień dwie główne gałęzie stronami (dla efektu wizualnego)
hc_cols$merge[1, ] <- hc_cols$merge[1, c(2, 1)]

# --- 5️⃣ Rysowanie heatmapy
ht <- Heatmap(
  mat,
  name = "Z-score (vs wtwt_saline)",
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
    title = "Z-score (vs control)",
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
