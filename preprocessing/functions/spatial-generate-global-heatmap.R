
spatial_generate_global_heatmap <- function(data,
                                    genes,
                                    metadata,
                                    summary_metric = "mean",
                                    data_type = "raw_data",
                                    group1 = NULL,
                                    group2 = NULL,
                                    colors_group1 = NULL,
                                    colors_group2 = NULL,
                                    title = NULL,
                                    scale_range = c(-4, 4),
                                    remove_empty_rows = TRUE,
                                    colors_on_heatmap = c("blue", "white", "red"),
                                    z_scale = TRUE,
                                    show_values = FALSE,
                                    verbose = FALSE) {
  
  #' Generate a heatmap of gene expression across samples
  #'
  #' This function visualizes global summary expression data (e.g., mean expression)
  #' for a given set of genes across samples using a heatmap. It supports optional z-score scaling,
  #' annotations for sample groups, custom color schemes, and clustered or fixed row/column order.
  #'
  #' @param data A nested list, e.g. from `compute_global_expression_multi()`, where
  #'   `data[[data_type]][[summary_metric]]` is a data.frame with expression summaries.
  #' @param genes Character vector of gene names to display on the heatmap.
  #' @param metadata A data.frame with at least `sample_ID` and optional grouping columns.
  #' @param summary_metric Character. The name of the summary metric to visualize (e.g., `"mean"`).
  #' @param data_type Character. Which `data_type` to use from the `data` list. Defaults to `"raw_data"`.
  #' @param group1 Character. First grouping variable (e.g., "genotype").
  #' @param group2 Character. Second grouping variable (e.g., "treatment").
  #' @param colors_group1 Named vector of colors for `group1` levels.
  #' @param colors_group2 Named vector of colors for `group2` levels.
  #' @param title Character. Heatmap title.
  #' @param scale_range Numeric vector of length 2. Min and max z-score range.
  #' @param remove_empty_rows Logical. Whether to remove genes with all NA values.
  #' @param colors_on_heatmap Color vector of length 3. Passed to colorRamp2.
  #' @param z_scale Logical. Whether to apply z-score scaling per gene (row).
  #' @param show_values Logical. If TRUE, displays numeric expression values in cells.
  #' @param verbose Logical. If TRUE, prints progress messages.
  #'
  #' @return A ComplexHeatmap object (invisible).
  #' @export
  
  if (verbose) message("[1] Filtering by genes...")
  
  summary_df <- data[[data_type]][[summary_metric]]
  
  df_filtered <- summary_df %>%
    dplyr::filter(gene %in% genes)
  
  if (verbose) message("[2] Pivoting to long format...")
  df_long <- df_filtered %>%
    tidyr::pivot_longer(
      cols = -c(peak, gene),
      names_to = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")
  
  if (verbose) message("[3] Scaling and reshaping...")
  df_scaled <- df_long %>%
    dplyr::group_by(peak) %>%
    dplyr::mutate(value = if (z_scale) scale(expression)[,1] else expression) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID = forcats::fct_inorder(sample_ID),
      peak = forcats::fct_rev(peak)
    )
  
  if (remove_empty_rows) {
    df_scaled <- df_scaled %>%
      dplyr::group_by(peak) %>%
      dplyr::filter(!all(is.na(value))) %>%
      dplyr::ungroup()
  }
  
  if (verbose) message("[4] Preparing annotations...")
  group_vars <- na.omit(c(group1, group2))
  meta2 <- df_scaled %>%
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
  
  if (verbose) message("[5] Building matrix...")
  mat <- df_scaled %>%
    dplyr::select(peak, sample_ID, value) %>%
    tidyr::pivot_wider(names_from = sample_ID, values_from = value) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]
  
  if (nrow(mat2) == 0 || ncol(mat2) == 0) {
    message(">>> No data to plot. Exiting.")
    return(invisible(NULL))
  }
  
  if (verbose) message("[6] Drawing heatmap...")
  if (z_scale) {
    breaks <- seq(scale_range[1], scale_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Z-score"
  } else {
    data_range <- range(mat2, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Expression"
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
  
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )
  
  invisible(ht)
}

spatial_generate_global_heatmap <- function(data,
                                            genes,
                                            metadata,
                                            summary_metric = "mean",
                                            data_type = "raw_data",
                                            group1 = NULL,
                                            group2 = NULL,
                                            colors_group1 = NULL,
                                            colors_group2 = NULL,
                                            title = NULL,
                                            scale_range = c(-4, 4),
                                            remove_empty_rows = TRUE,
                                            colors_on_heatmap = c("blue", "white", "red"),
                                            z_scale = TRUE,
                                            show_values = FALSE,
                                            select_max_peak = FALSE,
                                            aggregation_function = "mean",
                                            verbose = FALSE) {
  
  if (verbose) message("[1] Filtering by genes...")
  
  summary_df <- data[[data_type]][[summary_metric]]
  df_filtered <- summary_df %>%
    dplyr::filter(gene %in% genes)
  
  # [1.1] Select max peak per gene (simple vector version)
  if (select_max_peak) {
    if (verbose) message("[1.1] Selecting one peak per gene...")
    
    peak_selection <- df_filtered %>%
      tidyr::pivot_longer(cols = -c(peak, gene), names_to = "sample_ID", values_to = "expression") %>%
      dplyr::group_by(gene, peak) %>%
      dplyr::summarise(
        agg_expr = dplyr::case_when(
          aggregation_function == "mean"   ~ mean(expression, na.rm = TRUE),
          aggregation_function == "sum"    ~ sum(expression, na.rm = TRUE),
          aggregation_function == "median" ~ median(expression, na.rm = TRUE),
          TRUE ~ mean(expression, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_max(order_by = agg_expr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::pull(peak)
    
    df_filtered <- df_filtered %>%
      dplyr::filter(peak %in% peak_selection)
  }
  
  if (verbose) message("[2] Pivoting to long format...")
  df_long <- df_filtered %>%
    tidyr::pivot_longer(
      cols = -c(peak, gene),
      names_to = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")
  
  if (verbose) message("[3] Scaling and reshaping...")
  
  safe_scale <- function(x) {
    if (sum(!is.na(x)) <= 1) return(rep(NA_real_, length(x)))
    scale(x)[, 1]
  }
  
  df_scaled <- df_long %>%
    dplyr::group_by(peak) %>%
    dplyr::mutate(value = if (z_scale) safe_scale(expression) else expression) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID = forcats::fct_inorder(sample_ID),
      peak = forcats::fct_rev(peak)
    )
  
  if (remove_empty_rows) {
    df_scaled <- df_scaled %>%
      dplyr::group_by(peak) %>%
      dplyr::filter(!all(is.na(value))) %>%
      dplyr::ungroup()
  }
  
  if (verbose) message("[4] Preparing annotations...")
  group_vars <- na.omit(c(group1, group2))
  meta2 <- df_scaled %>%
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
  
  if (verbose) message("[5] Building matrix...")
  mat <- df_scaled %>%
    dplyr::select(peak, sample_ID, value) %>%
    tidyr::pivot_wider(names_from = sample_ID, values_from = value) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]
  
  if (nrow(mat2) == 0 || ncol(mat2) == 0) {
    message(">>> No data to plot. Exiting.")
    return(invisible(NULL))
  }
  
  if (verbose) message("[6] Drawing heatmap...")
  if (z_scale) {
    breaks <- seq(scale_range[1], scale_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Z-score"
  } else {
    data_range <- range(mat2, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Expression"
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
  
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )
  
  invisible(ht)
}


spatial_generate_global_heatmap <- function(data,
                                            genes,
                                            metadata,
                                            summary_metric = "mean",
                                            data_type = "raw_data",
                                            group1 = NULL,
                                            group2 = NULL,
                                            colors_group1 = NULL,
                                            colors_group2 = NULL,
                                            title = NULL,
                                            scale_range = c(-4, 4),
                                            remove_empty_rows = TRUE,
                                            colors_on_heatmap = c("blue", "white", "red"),
                                            z_scale = TRUE,
                                            show_values = FALSE,
                                            select_max_peak = FALSE,
                                            threshold_expr = NULL,
                                            aggregation_function = "mean",
                                            verbose = FALSE) {
  
  if (verbose) message("[1] Filtering by genes...")
  
  summary_df <- data[[data_type]][[summary_metric]]
  df_filtered <- summary_df %>%
    dplyr::filter(gene %in% genes)
  
  # [1.1] Select max peak per gene (exclusive mode)
  if (select_max_peak) {
    if (verbose) message("[1.1] Selecting one peak per gene...")
    peak_selection <- df_filtered %>%
      tidyr::pivot_longer(cols = -c(peak, gene), names_to = "sample_ID", values_to = "expression") %>%
      dplyr::group_by(gene, peak) %>%
      dplyr::summarise(
        agg_expr = dplyr::case_when(
          aggregation_function == "mean"   ~ mean(expression, na.rm = TRUE),
          aggregation_function == "sum"    ~ sum(expression, na.rm = TRUE),
          aggregation_function == "median" ~ median(expression, na.rm = TRUE),
          TRUE ~ mean(expression, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_max(order_by = agg_expr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::pull(peak)
    
    df_filtered <- df_filtered %>%
      dplyr::filter(peak %in% peak_selection)
  }
  
  # [1.2] Select all peaks passing threshold (only if not using max-peak)
  if (!select_max_peak && !is.null(threshold_expr)) {
    if (verbose) message("[1.2] Selecting peaks with aggregated expression above threshold...")
    
    peak_selection <- df_filtered %>%
      tidyr::pivot_longer(cols = -c(peak, gene), names_to = "sample_ID", values_to = "expression") %>%
      dplyr::group_by(gene, peak) %>%
      dplyr::summarise(
        agg_expr = dplyr::case_when(
          aggregation_function == "mean"   ~ mean(expression, na.rm = TRUE),
          aggregation_function == "sum"    ~ sum(expression, na.rm = TRUE),
          aggregation_function == "median" ~ median(expression, na.rm = TRUE),
          TRUE ~ mean(expression, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      dplyr::filter(agg_expr >= threshold_expr) %>%
      dplyr::pull(peak)
    
    df_filtered <- df_filtered %>%
      dplyr::filter(peak %in% peak_selection)
  }
  
  if (verbose) message("[2] Pivoting to long format...")
  df_long <- df_filtered %>%
    tidyr::pivot_longer(
      cols = -c(peak, gene),
      names_to = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")
  
  if (verbose) message("[3] Scaling and reshaping...")
  
  safe_scale <- function(x) {
    if (sum(!is.na(x)) <= 1) return(rep(NA_real_, length(x)))
    scale(x)[, 1]
  }
  
  df_scaled <- df_long %>%
    dplyr::group_by(peak) %>%
    dplyr::mutate(value = if (z_scale) safe_scale(expression) else expression) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_ID = forcats::fct_inorder(sample_ID),
      peak = forcats::fct_rev(peak)
    )
  
  if (remove_empty_rows) {
    df_scaled <- df_scaled %>%
      dplyr::group_by(peak) %>%
      dplyr::filter(!all(is.na(value))) %>%
      dplyr::ungroup()
  }
  
  if (verbose) message("[4] Preparing annotations...")
  group_vars <- na.omit(c(group1, group2))
  meta2 <- df_scaled %>%
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
  
  if (verbose) message("[5] Building matrix...")
  mat <- df_scaled %>%
    dplyr::select(peak, sample_ID, value) %>%
    tidyr::pivot_wider(names_from = sample_ID, values_from = value) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  mat2 <- mat[, meta2$sample_ID, drop = FALSE]
  
  if (nrow(mat2) == 0 || ncol(mat2) == 0) {
    message(">>> No data to plot. Exiting.")
    return(invisible(NULL))
  }
  
  if (verbose) message("[6] Drawing heatmap...")
  if (z_scale) {
    breaks <- seq(scale_range[1], scale_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Z-score"
  } else {
    data_range <- range(mat2, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Expression"
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
  
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    merge_legend           = TRUE
  )
  
  invisible(ht)
}

