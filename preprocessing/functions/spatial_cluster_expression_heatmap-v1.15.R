spatial_generate_cluster_heatmap <- function(data, genes, cluster_name, metadata,
                                     group1 = NULL, group2 = NULL,
                                     colors_group1 = NULL, colors_group2 = NULL,
                                     data_type = "raw_data", resolution = "resolution_0.4",
                                     summary_metric = "mean", title = NULL,
                                     scale_range = c(-4, 4),
                                     remove_empty_rows = TRUE,
                                     colors_on_heatmap = c("blue", "white", "red"),
                                     z_scale = TRUE,
                                     show_values = FALSE,
                                     verbose = FALSE) {
  if (verbose) message("[1] Loading raw data...")
  tmp <- data[[data_type]][[resolution]]
  cluster_names <- names(tmp) %>%
    tibble::tibble(name = .) %>%
    dplyr::mutate(index = as.integer(stringr::str_extract(name, "\\d+$"))) %>%
    dplyr::arrange(index) %>%
    dplyr::pull(name)
  
  if (verbose) message("[2] Building heatmap dataframe...")
  heatmap_df <- purrr::map_dfr(cluster_names, function(cl) {
    df <- cbind(
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
    if (verbose) message(sprintf("  - Processed %s: %d rows", cl, nrow(df)))
    df
  }) %>%
    tidyr::pivot_longer(
      cols      = -c(cluster, peak, gene),
      names_to  = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")
  if (verbose) print(utils::head(heatmap_df))
  
  if (verbose) message("[3] Scaling values and filtering for cluster...")
  df_cluster <- heatmap_df %>%
    dplyr::filter(cluster == cluster_name) %>%
    dplyr::group_by(cluster, peak) %>%
    dplyr::mutate(
      value = if (z_scale) base::scale(expression)[,1] else expression
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID = forcats::fct_inorder(sample_ID),
      peak      = forcats::fct_rev(peak)
    )
  if (verbose) {
    # Summary of NA counts after scaling
    peak_counts <- df_cluster %>%
      dplyr::group_by(peak) %>%
      dplyr::summarise(n_non_na = sum(!is.na(value)),
                       total = dplyr::n())
    message("[3a] Value counts per peak after scaling:")
    print(utils::head(peak_counts))
    zero_peaks <- peak_counts$peak[peak_counts$n_non_na == 0]
    if (length(zero_peaks) > 0) {
      message("Peaks with no non-NA values after scaling: ", paste(zero_peaks, collapse = ", "))
    }
  }
  
  if (verbose) message("[4] Removing empty rows...")
  df_cluster <- df_cluster %>%
    dplyr::group_by(peak) %>%
    dplyr::filter(
      !(remove_empty_rows && base::all(base::is.na(value)))
    ) %>%
    dplyr::ungroup()
  
  if (verbose) message("[5] Preparing metadata annotation...")
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
      } else {.}
    }
  if (verbose) print(utils::head(meta2))
  
  if (verbose) message("[6] Constructing matrix...")
  mat <- df_cluster %>%
    dplyr::select(peak, sample_ID, value) %>%
    tidyr::pivot_wider(
      names_from  = sample_ID,
      values_from = value
    ) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]
  if (verbose) message(sprintf("  Matrix dimensions: %d x %d", nrow(mat2), ncol(mat2)))
  if (nrow(mat2) == 0 || ncol(mat2) == 0) {
    if (verbose) message(sprintf("No data available for cluster '%s'. Exiting.", cluster_name))
    return(invisible(NULL))
  }
  
  if (verbose) message("[7] Setting up heatmap parameters...")
  if (z_scale) {
    breaks <- seq(scale_range[1], scale_range[2], length.out = length(colors_on_heatmap))
    legend_name <- paste0(toupper(substring(summary_metric, 1, 1)), substring(summary_metric, 2), " Z-score")
  } else {
    data_range <- range(mat2, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = length(colors_on_heatmap))
    legend_name <- summary_metric
  }
  hm_args <- list(
    mat2,
    name                 = legend_name,
    col                  = circlize::colorRamp2(breaks, colors_on_heatmap),
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
  
  if (show_values) {
    hm_args$cell_fun <- function(j, i, x, y, width, height, fill) {
      val <- mat2[i, j]
      txt <- if (is.na(val)) "NA" else sprintf("%.2f", val)
      grid::grid.text(txt, x, y, gp = grid::gpar(fontsize = 6))
    }
  }
  
  if (length(group_vars) > 0) {
    if (verbose) message("[8] Adding group annotations...")
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
  
  if (verbose) message("[9] Drawing heatmap...")
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )
  
  invisible(ht)
}
