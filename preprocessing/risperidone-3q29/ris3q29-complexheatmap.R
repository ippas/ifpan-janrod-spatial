# 0) Przygotowanie surowych danych
tmp <- wtDel_summary_statistics$quantile_normalize$resolution_0.4

genes_del3q29 <- c(
  "Bdh1", "Dlg1", "Mfi2", "Pigz", "Ncbp2", "Senp5", "Pak2", "Pigx",
  "Cep19", "Nrros", "Bex6",  "Fbxo45", "Wdr53", "Smco1", "Rnf168",
  "Ubxn7", "Tm4sf19", "Tctex1d2", "Pcyt1a", "Slc51a", "Zdhhc19", "Tfrc"
)

cluster_names <- names(tmp) %>%
  tibble::tibble(name = .) %>%
  dplyr::mutate(index = as.integer(stringr::str_extract(name, "\\d+$"))) %>%
  dplyr::arrange(index) %>%
  dplyr::pull(name)

heatmap_df <- purrr::map_dfr(cluster_names, function(cl) {
  cbind(
    tmp[[cl]]$control$mean,
    tmp[[cl]]$experiment$mean
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(
      gene    = stringr::str_extract(peak, "[^-]+$"),
      cluster = cl
    ) %>%
    dplyr::filter(gene %in% genes_del3q29)
}) %>%
  tidyr::pivot_longer(
    cols      = -c(cluster, peak, gene),
    names_to  = "sample_ID",
    values_to = "expression"
  ) %>%
  dplyr::left_join(sample_info, by = "sample_ID")


df_cluster0 <- heatmap_df %>%
  dplyr::filter(cluster == "cluster_0") %>%
  dplyr::group_by(cluster, peak) %>%
  dplyr::mutate(
    scaled_expression = stats::scale(expression)[, 1]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sample_ID      = forcats::fct_inorder(sample_ID),
    peak           = forcats::fct_rev(peak),
    mouse_genotype = base::factor(mouse_genotype, levels = c("wtwt", "wtdel"))
  ) %>%
  dplyr::group_by(peak) %>%
  dplyr::filter(
    !base::all(base::is.na(scaled_expression))
  ) %>%
  dplyr::ungroup()



#################################################################################
# 1) meta2 – już masz uporządkowane wg genotypu i treatment:
meta2 <- df_cluster0 %>%
  distinct(sample_ID, mouse_genotype, treatment) %>%
  mutate(
    mouse_genotype = factor(mouse_genotype, levels = c("wtwt","wtdel")),
    treatment      = factor(treatment,      levels = c("saline","risperidone"))
  ) %>%
  arrange(mouse_genotype, treatment)

# 2) Przeporzadkuj macierz ekspresji:
mat2 <- mat[, meta2$sample_ID, drop = FALSE]

# 3) Zbuduj top_annotation w tej samej kolejności:
col_ann2 <- ComplexHeatmap::HeatmapAnnotation(
  df = data.frame(
    Genotype  = meta2$mouse_genotype,
    Treatment = meta2$treatment
  ),
  which                  = "column",
  col = list(
    Genotype  = c(wtwt = "#1f77b4", wtdel = "#ff7f0e"),
    Treatment = c(saline = "#6baed6", risperidone = "#fd8d3c")
  ),
  annotation_name_side   = "left",
  annotation_legend_param= list(direction = "horizontal", nrow = 1)
)

# 4) Narysuj heatmapę bez żadnego column_order:
ht <- ComplexHeatmap::Heatmap(
  mat2,
  name                 = "Z-score",
  col                  = circlize::colorRamp2(c(-4, 0, 4),
                                              c("blue","white","red")),
  na_col               = "grey80",
  top_annotation       = col_ann2,
  show_row_names       = TRUE,   row_names_side    = "left",
  show_column_names    = TRUE,   column_names_side = "top",
  cluster_rows         = FALSE,  cluster_columns   = FALSE,
  heatmap_legend_param = list(direction = "horizontal", nrow = 1)
)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  merge_legend           = TRUE
)


