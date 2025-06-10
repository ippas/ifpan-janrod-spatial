generate_cluster_heatmap <- function(data, genes, cluster_name, metadata,
                                     group1 = NULL, group2 = NULL,
                                     colors_group1 = NULL, colors_group2 = NULL,
                                     data_type = "raw_data", resolution = "resolution_0.4",
                                     summary_metric = "mean", title = NULL) {
  # 1) Przygotowanie surowych danych
  tmp <- data[[data_type]][[resolution]]
  cluster_names <- names(tmp) %>%
    tibble::tibble(name = .) %>%
    dplyr::mutate(index = as.integer(stringr::str_extract(name, "\\d+$"))) %>%
    dplyr::arrange(index) %>%
    dplyr::pull(name)
  
  heatmap_df <- purrr::map_dfr(cluster_names, function(cl) {
    cbind(
      tmp[[cl]]$control[[summary_metric]],
      tmp[[cl]]$experiment[[summary_metric]]
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
    dplyr::mutate(scaled_expression = base::scale(expression)[,1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID = forcats::fct_inorder(sample_ID),
      peak      = forcats::fct_rev(peak)
    ) %>%
    dplyr::group_by(peak) %>%
    dplyr::filter(!base::all(base::is.na(scaled_expression))) %>%
    dplyr::ungroup()
  
  # 3) Przygotowanie meta2 z dynamicznym grupowaniem
  group_vars <- na.omit(c(group1, group2))
  meta2 <- df_cluster %>%
    dplyr::distinct(sample_ID, dplyr::across(dplyr::all_of(group_vars))) %>%
    {
      df <- .
      if (!is.null(group1)) df[[group1]] <- forcats::fct_inorder(df[[group1]])
      if (!is.null(group2)) df[[group2]] <- forcats::fct_inorder(df[[group2]])
      df
    } %>%
    {
      if (length(group_vars) == 2) {
        dplyr::arrange(., .data[[group1]], .data[[group2]])
      } else if (length(group_vars) == 1) {
        dplyr::arrange(., .data[[group_vars]])
      } else {
        .
      }
    }
  
  # 4) Budowa macierzy ekspresji w porządku meta2
  mat <- df_cluster %>%
    dplyr::select(peak, sample_ID, scaled_expression) %>%
    tidyr::pivot_wider(
      names_from  = sample_ID,
      values_from = scaled_expression
    ) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]
  
  # 5) Przygotowanie Heatmap arguments z możliwością tytułu
  hm_args <- list(
    mat2,
    name                 = paste0(toupper(substring(summary_metric,1,1)), substring(summary_metric,2), " Z-score"),
    col                  = circlize::colorRamp2(c(-4,0,4), c("blue","white","red")),
    na_col               = "grey80",
    show_row_names       = TRUE,
    row_names_side       = "left",
    show_column_names    = TRUE,
    column_names_side    = "top",
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    column_title         = title,
    heatmap_legend_param = list(direction = "horizontal", nrow = 1)
  )
  # 6) Dodanie adnotacji z kolorami grup
  if (length(group_vars) > 0) {
    ann_df <- as.data.frame(meta2[dplyr::all_of(group_vars)])
    col_list <- list()
    if (!is.null(group1) && !is.null(colors_group1)) col_list[[group1]] <- colors_group1
    if (!is.null(group2) && !is.null(colors_group2)) col_list[[group2]] <- colors_group2
    hm_args$top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df                   = ann_df,
      which                = "column",
      col                  = col_list,
      annotation_name_side = "left"
    )
  }
  
  # 7) Rysowanie heatmapy
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )
  
  invisible(ht)
}


generate_cluster_heatmap(
  data         = wtDel_summary_statistics,
  genes        = genes_del3q29,
  data_type = "quantile_normalize",
  resolution = "resolution_0.4",
  summary_metric = "sum",
  cluster_name = "cluster_0",
  metadata     = sample_info,
  group1       = "mouse_genotype",
  group2       = "treatment",
  colors_group1 = c("wtwt" = "red", "wtdel" = "blue"),
  colors_group2 = c("saline" = "orange", "risperidone" =  "green"),
  title = "Cluster 0; mean"
)
