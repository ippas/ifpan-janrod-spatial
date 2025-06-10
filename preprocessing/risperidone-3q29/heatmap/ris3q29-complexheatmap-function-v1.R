generate_cluster_heatmap <- function(data, genes, cluster_name, metadata) {
  # 1) Przygotowanie surowych danych
  tmp <- data$quantile_normalize$resolution_0.4
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
      dplyr::filter(gene %in% genes)
  }) %>%
    tidyr::pivot_longer(
      cols      = -c(cluster, peak, gene),
      names_to  = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")

  # 2) Filtracja clusteru i skalowanie ekspresji
  df_cluster <- heatmap_df %>%
    dplyr::filter(cluster == cluster_name) %>%
    dplyr::group_by(cluster, peak) %>%
    dplyr::mutate(
      scaled_expression = base::scale(expression)[, 1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID      = forcats::fct_inorder(sample_ID),
      peak           = forcats::fct_rev(peak),
      mouse_genotype = base::factor(mouse_genotype, levels = c("wtwt", "wtdel"))
    ) %>%
    dplyr::group_by(peak) %>%
    dplyr::filter(!base::all(base::is.na(scaled_expression))) %>%
    dplyr::ungroup()

  # 3) Przygotowanie meta2
  meta2 <- df_cluster %>%
    dplyr::distinct(sample_ID, mouse_genotype, treatment) %>%
    dplyr::mutate(
      mouse_genotype = base::factor(mouse_genotype, levels = c("wtwt", "wtdel")),
      treatment      = base::factor(treatment,      levels = c("saline", "risperidone"))
    ) %>%
    dplyr::arrange(mouse_genotype, treatment)

  # 4) Budowa macierzy ekspresji
  mat <- df_cluster %>%
    dplyr::select(peak, sample_ID, scaled_expression) %>%
    tidyr::pivot_wider(
      names_from   = sample_ID,
      values_from  = scaled_expression
    ) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]

  # 5) Tworzenie adnotacji kolumn
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

  # 6) Rysowanie heatmapy
  ht <- ComplexHeatmap::Heatmap(
    mat2,
    name                 = "Z-score",
    col                  = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
    na_col               = "grey80",
    top_annotation       = col_ann2,
    show_row_names       = TRUE,  row_names_side    = "left",
    show_column_names    = TRUE,  column_names_side = "top",
    cluster_rows         = FALSE, cluster_columns   = FALSE,
    heatmap_legend_param = list(direction = "horizontal", nrow = 1)
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )

  invisible(ht)
}


# Wywołanie funkcji dla klastra 0
generate_cluster_heatmap(
  data          = wtDel_summary_statistics,
  genes         = genes_del3q29,
  cluster_name  = "cluster_7",
  metadata      = sample_info
)

