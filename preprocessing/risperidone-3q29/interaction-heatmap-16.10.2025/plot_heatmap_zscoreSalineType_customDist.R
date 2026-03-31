plot_heatmap_zscore_salinetype <- function(
  data_to_heatmap_test,
  metadata_ris3q29,
  base_distance = "2 - r",        # "2 - r" lub "1 - r"
  r_scale = 1,                    # siła wpływu r
  same_cluster_bonus = 0.2,       # 0..1 (procent); 0.2 = ±20%
  hclust_method = "complete",
  scale_limits = c(-3, 3),
  clip_outliers = TRUE,
  verbose_limit = 6,              # ile pierwszych wierszy/kolumn logować (diagnostycznie)
  row_dend_width = 25,            # 📏 szerokość dendrogramu w mm
  rev_heatmap = FALSE             # 🔁 odwrócenie kolejności na heatmapie
) {
  # ===== Pakiety
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(tibble)
  
  message("📘 [1] Przygotowanie danych i obliczenie Z-score względem saline...")
  
  # ===== 1. Dane + z-score per genotype vs saline
  df <- data_to_heatmap_test %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    tidyr::separate(rowname, into = c("cluster", "name_id"), sep = "\\|") %>%
    tidyr::separate(name_id, into = c("peak_id", "gene_symbol"),
                    sep = "-(?=[^-]+$)", remove = FALSE) %>%
    tidyr::gather(key = "sample_ID", value = "value",
                  -cluster, -name_id, -peak_id, -gene_symbol) %>%
    left_join(metadata_ris3q29[, c("sample_ID", "mouse_genotype", "treatment")],
              by = "sample_ID") %>%
    mutate(group = paste0(mouse_genotype, "_", treatment)) %>%
    group_by(cluster, name_id, peak_id, gene_symbol, mouse_genotype) %>%
    mutate(
      mean_ctrl = mean(value[treatment == "saline"], na.rm = TRUE),
      sd_ctrl   = sd(value[treatment == "saline"], na.rm = TRUE),
      sd_ctrl   = ifelse(is.na(sd_ctrl) | sd_ctrl == 0, 1, sd_ctrl),
      zscore    = (value - mean_ctrl) / sd_ctrl
    ) %>%
    ungroup() %>%
    mutate(
      cluster_nameID = paste0(cluster, "_", name_id),
      column_ID = paste0(sample_ID, "_", group)
    )
  
  # ===== 2. Macierz do heatmapy
  message("📘 [2] Tworzenie macierzy Z-score do heatmapy...")
  mat <- df %>%
    select(cluster_nameID, column_ID, zscore) %>%
    tidyr::pivot_wider(names_from = column_ID, values_from = zscore) %>%
    column_to_rownames("cluster_nameID") %>%
    as.matrix()
  mat[is.na(mat)] <- 0
  
  # ===== 3. Przycinanie wartości
  if (clip_outliers) {
    message("📘 [3] Przycinanie wartości do skali...")
    mat[mat > scale_limits[2]] <- scale_limits[2]
    mat[mat < scale_limits[1]] <- scale_limits[1]
  }
  
  # ===== 4. Kolejność kolumn (grup)
  message("📘 [4] Ustalanie kolejności grup próbek...")
  group_order <- c("wtwt_saline", "wtwt_risperidone", "wtdel_saline", "wtdel_risperidone")
  groups <- gsub(".*_(wtwt_saline|wtwt_risperidone|wtdel_saline|wtdel_risperidone)$",
                 "\\1", colnames(mat))
  group_colors <- c(
    "wtwt_saline"      = "#66c2a5",
    "wtwt_risperidone" = "#fc8d62",
    "wtdel_saline"     = "#8da0cb",
    "wtdel_risperidone"= "#e78ac3"
  )
  col_order <- order(factor(groups, levels = group_order))
  mat <- mat[, col_order, drop = FALSE]
  groups <- groups[col_order]
  ha_col <- HeatmapAnnotation(
    Group = factor(groups, levels = group_order),
    col = list(Group = group_colors),
    annotation_height = unit(4, "mm")
  )
  
  # ===== 5. Dystanse
  message("📘 [5] Obliczanie niestandardowego dystansu...")
  cor_mat <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
  rn <- rownames(mat)
  
  # Etykiety klastrów
  m <- regexpr("^cluster_[0-9]+", rn)
  cluster_labels <- rep(NA_character_, length(rn))
  cluster_labels[m > 0] <- regmatches(rn, m)
  if (any(is.na(cluster_labels))) {
    cluster_labels[is.na(cluster_labels)] <- sub("^(([^_]+_[0-9]+)).*$", "\\1", rn[is.na(cluster_labels)])
  }
  
  # Bazowy dystans
  base_dist <- if (base_distance == "2 - r") 2 - (cor_mat * r_scale)
  else if (base_distance == "1 - r") 1 - (cor_mat * r_scale)
  else stop("base_distance musi być '2 - r' lub '1 - r'")
  
  same_mat <- outer(cluster_labels, cluster_labels, FUN = "==")
  dist_mat <- base_dist
  dist_mat[same_mat]  <- base_dist[same_mat]  * (1 - same_cluster_bonus)
  dist_mat[!same_mat] <- base_dist[!same_mat] * (1 + same_cluster_bonus)
  diag(dist_mat) <- 0
  dist_mat[dist_mat < 0] <- 0
  
  # ===== 5j. Log diagnostyczny – podgląd przykładowych dystansów
  n_show <- min(verbose_limit, nrow(dist_mat))
  message("🔎 Przykładowe dystanse (", n_show, "×", n_show, ") z logiem zmian:")
  for (i in seq_len(n_show)) {
    for (j in seq_len(n_show)) {
      gi <- rownames(dist_mat)[i]; gj <- colnames(dist_mat)[j]
      same <- same_mat[i, j]
      if (same) {
        message(sprintf("🟢 Ten sam klaster (%s): %s vs %s | %.3f → %.3f",
                        cluster_labels[i], gi, gj,
                        base_dist[i, j], dist_mat[i, j]))
      } else {
        message(sprintf("🔴 Różne klastry (%s vs %s): %s vs %s | %.3f → %.3f",
                        cluster_labels[i], cluster_labels[j],
                        gi, gj, base_dist[i, j], dist_mat[i, j]))
      }
    }
  }
  message(sprintf("📊 Zakres dystansów: min=%.3f, max=%.3f",
                  suppressWarnings(min(dist_mat[is.finite(dist_mat)])),
                  suppressWarnings(max(dist_mat[is.finite(dist_mat)]))))
  
  # ===== 6. Klasteryzacja
  message("📘 [6] Klasteryzacja metodą ", hclust_method, "...")
  hc_rows <- hclust(as.dist(dist_mat), method = hclust_method)
  hc_rows <- as.dendrogram(hc_rows)
  
  # 🔁 Odwrócenie kolejności (jeśli TRUE)
  if (rev_heatmap) {
    message("🔄 [Opcja] Odwracanie kolejności wierszy heatmapy (rev_heatmap = TRUE)...")
    hc_rows <- rev(hc_rows)
  }
  
  # ===== 7. Heatmapa
  message("📘 [7] Rysowanie heatmapy...")
  col_fun <- colorRamp2(
    c(scale_limits[1], 0, scale_limits[2]),
    c("navy", "white", "firebrick3")
  )
  legend_breaks <- c(scale_limits[1], 0, scale_limits[2])
  
  ht <- Heatmap(
    mat,
    name = "Z-score (vs saline per genotype)",
    col = col_fun,
    cluster_rows = hc_rows,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10, rot = 90),
    top_annotation = ha_col,
    border = "black",
    column_split = factor(groups, levels = group_order),
    column_gap = unit(0, "mm"),
    row_gap = unit(0, "mm"),
    rect_gp = gpar(col = NA),
    row_dend_width = unit(row_dend_width, "mm"),   # 📏 szerokość dendrogramu
    heatmap_legend_param = list(
      title = "Z-score (vs saline)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      at = legend_breaks,
      labels = as.character(legend_breaks)
    ),
    column_title = "Sample ID and Group (Normalized to saline per genotype)",
    row_title = "Cluster | Name ID",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
  
  message("✅ Zakończono generowanie heatmapy.")
}


plot_heatmap_zscore_salinetype <- function(
  data_to_heatmap_test,
  metadata_ris3q29,
  base_distance = "2 - r",        # "2 - r" lub "1 - r"
  r_scale = 1,                    # siła wpływu r
  same_cluster_bonus = 0.2,       # 0..1 (procent); 0.2 = ±20%
  hclust_method = "complete",
  scale_limits = c(-3, 3),
  clip_outliers = TRUE,
  verbose_limit = 6,              # ile pierwszych wierszy/kolumn logować (diagnostycznie)
  row_dend_width = 25,            # szerokość dendrogramu w mm
  rev_heatmap = FALSE,            # odwrócenie kolejności na heatmapie
  heatmap_colors = c("navy", "white", "firebrick3")  # ✅ NOWY ARGUMENT: low, mid, high
) {
  # ===== Pakiety
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(tibble)
  library(grid)
  
  # ===== Walidacja kolorów
  if (!is.character(heatmap_colors) || length(heatmap_colors) != 3) {
    stop("heatmap_colors musi być wektorem 3 kolorów: c(low, mid, high).")
  }
  
  message("📘 [1] Przygotowanie danych i obliczenie Z-score względem saline...")
  
  # ===== 1. Dane + z-score per genotype vs saline
  df <- data_to_heatmap_test %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    tidyr::separate(rowname, into = c("cluster", "name_id"), sep = "\\|") %>%
    tidyr::separate(name_id, into = c("peak_id", "gene_symbol"),
                    sep = "-(?=[^-]+$)", remove = FALSE) %>%
    tidyr::gather(key = "sample_ID", value = "value",
                  -cluster, -name_id, -peak_id, -gene_symbol) %>%
    left_join(metadata_ris3q29[, c("sample_ID", "mouse_genotype", "treatment")],
              by = "sample_ID") %>%
    mutate(group = paste0(mouse_genotype, "_", treatment)) %>%
    group_by(cluster, name_id, peak_id, gene_symbol, mouse_genotype) %>%
    mutate(
      mean_ctrl = mean(value[treatment == "saline"], na.rm = TRUE),
      sd_ctrl   = sd(value[treatment == "saline"], na.rm = TRUE),
      sd_ctrl   = ifelse(is.na(sd_ctrl) | sd_ctrl == 0, 1, sd_ctrl),
      zscore    = (value - mean_ctrl) / sd_ctrl
    ) %>%
    ungroup() %>%
    mutate(
      cluster_nameID = paste0(cluster, "_", name_id),
      column_ID = paste0(sample_ID, "_", group)
    )
  
  # ===== 2. Macierz do heatmapy
  message("📘 [2] Tworzenie macierzy Z-score do heatmapy...")
  mat <- df %>%
    select(cluster_nameID, column_ID, zscore) %>%
    tidyr::pivot_wider(names_from = column_ID, values_from = zscore) %>%
    column_to_rownames("cluster_nameID") %>%
    as.matrix()
  mat[is.na(mat)] <- 0
  
  # ===== 3. Przycinanie wartości
  if (clip_outliers) {
    message("📘 [3] Przycinanie wartości do skali...")
    mat[mat > scale_limits[2]] <- scale_limits[2]
    mat[mat < scale_limits[1]] <- scale_limits[1]
  }
  
  # ===== 4. Kolejność kolumn (grup)
  message("📘 [4] Ustalanie kolejności grup próbek...")
  group_order <- c("wtwt_saline", "wtwt_risperidone", "wtdel_saline", "wtdel_risperidone")
  groups <- gsub(".*_(wtwt_saline|wtwt_risperidone|wtdel_saline|wtdel_risperidone)$",
                 "\\1", colnames(mat))
  group_colors <- c(
    "wtwt_saline"       = "#66c2a5",
    "wtwt_risperidone"  = "#fc8d62",
    "wtdel_saline"      = "#8da0cb",
    "wtdel_risperidone" = "#e78ac3"
  )
  col_order <- order(factor(groups, levels = group_order))
  mat <- mat[, col_order, drop = FALSE]
  groups <- groups[col_order]
  
  ha_col <- HeatmapAnnotation(
    Group = factor(groups, levels = group_order),
    col = list(Group = group_colors),
    annotation_height = unit(4, "mm")
  )
  
  # ===== 5. Dystanse
  message("📘 [5] Obliczanie niestandardowego dystansu...")
  cor_mat <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
  rn <- rownames(mat)
  
  # Etykiety klastrów
  m <- regexpr("^cluster_[0-9]+", rn)
  cluster_labels <- rep(NA_character_, length(rn))
  cluster_labels[m > 0] <- regmatches(rn, m)
  if (any(is.na(cluster_labels))) {
    cluster_labels[is.na(cluster_labels)] <- sub("^(([^_]+_[0-9]+)).*$", "\\1", rn[is.na(cluster_labels)])
  }
  
  # Bazowy dystans
  base_dist <- if (base_distance == "2 - r") {
    2 - (cor_mat * r_scale)
  } else if (base_distance == "1 - r") {
    1 - (cor_mat * r_scale)
  } else {
    stop("base_distance musi być '2 - r' lub '1 - r'")
  }
  
  same_mat <- outer(cluster_labels, cluster_labels, FUN = "==")
  dist_mat <- base_dist
  dist_mat[same_mat]  <- base_dist[same_mat]  * (1 - same_cluster_bonus)
  dist_mat[!same_mat] <- base_dist[!same_mat] * (1 + same_cluster_bonus)
  diag(dist_mat) <- 0
  dist_mat[dist_mat < 0] <- 0
  
  # ===== 5j. Log diagnostyczny – podgląd przykładowych dystansów
  n_show <- min(verbose_limit, nrow(dist_mat))
  message("🔎 Przykładowe dystanse (", n_show, "×", n_show, ") z logiem zmian:")
  for (i in seq_len(n_show)) {
    for (j in seq_len(n_show)) {
      gi <- rownames(dist_mat)[i]; gj <- colnames(dist_mat)[j]
      same <- same_mat[i, j]
      if (same) {
        message(sprintf("🟢 Ten sam klaster (%s): %s vs %s | %.3f → %.3f",
                        cluster_labels[i], gi, gj,
                        base_dist[i, j], dist_mat[i, j]))
      } else {
        message(sprintf("🔴 Różne klastry (%s vs %s): %s vs %s | %.3f → %.3f",
                        cluster_labels[i], cluster_labels[j],
                        gi, gj, base_dist[i, j], dist_mat[i, j]))
      }
    }
  }
  message(sprintf("📊 Zakres dystansów: min=%.3f, max=%.3f",
                  suppressWarnings(min(dist_mat[is.finite(dist_mat)])),
                  suppressWarnings(max(dist_mat[is.finite(dist_mat)]))))
  
  # ===== 6. Klasteryzacja
  message("📘 [6] Klasteryzacja metodą ", hclust_method, "...")
  hc_rows <- hclust(as.dist(dist_mat), method = hclust_method)
  hc_rows <- as.dendrogram(hc_rows)
  
  # Odwrócenie kolejności (jeśli TRUE)
  if (rev_heatmap) {
    message("🔄 [Opcja] Odwracanie kolejności wierszy heatmapy (rev_heatmap = TRUE)...")
    hc_rows <- rev(hc_rows)
  }
  
  # ===== 7. Heatmapa
  message("📘 [7] Rysowanie heatmapy...")
  col_fun <- colorRamp2(
    c(scale_limits[1], 0, scale_limits[2]),
    heatmap_colors
  )
  legend_breaks <- c(scale_limits[1], 0, scale_limits[2])
  
  ht <- Heatmap(
    mat,
    name = "Z-score (vs saline per genotype)",
    col = col_fun,
    cluster_rows = hc_rows,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10, rot = 90),
    top_annotation = ha_col,
    border = "black",
    column_split = factor(groups, levels = group_order),
    column_gap = unit(0, "mm"),
    row_gap = unit(0, "mm"),
    rect_gp = gpar(col = NA),
    row_dend_width = unit(row_dend_width, "mm"),
    heatmap_legend_param = list(
      title = "Z-score (vs saline)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      at = legend_breaks,
      labels = as.character(legend_breaks)
    ),
    column_title = "Sample ID and Group (Normalized to saline per genotype)",
    row_title = "Cluster | Name ID",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
  
  message("✅ Zakończono generowanie heatmapy.")
}

# ##############################################################################

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
  impute_by_group = TRUE,
  add_marker_for_imputed_data = TRUE,
  log2_transform = FALSE,
  return_matrix = TRUE,
  scale_by_row = FALSE
)

data_to_heatmap_test %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  separate(id, into = c("cluster", "rest"), sep = "\\|") %>% 
  separate(rest, into = c("peak", "gene_symbol"), sep = "-(?=[^-]+$)") %>%  filter(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE))
  .$gene_symbol %>% duplicated()




plot_heatmap_zscore_salinetype(
  data_to_heatmap_test,
  metadata_ris3q29,
  same_cluster_bonus = 0.25,
  r_scale = 1,
  row_dend_width = 30,
  scale_limits = c(-6, 6),
  rev_heatmap = T
)

svg("results/risperidone-3q29/figures/heatmap-interactionEdgeR-16.10.2025/heatmap_interaction_zscoreRefSaline_dist1rSameCluster0.25.svg", width = 10, height = 16)
plot_heatmap_zscore_salinetype(
  data_to_heatmap_test,
  metadata_ris3q29,
  same_cluster_bonus = 0.25,
  r_scale = 1,
  row_dend_width = 30,
  scale_limits = c(-3, 3),
  rev_heatmap = TRUE
)
dev.off()
