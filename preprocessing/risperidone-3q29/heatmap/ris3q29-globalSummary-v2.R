compute_global_expression_summary <- function(spatial_data,
                                              data_type = "raw_data",
                                              control_samples,
                                              experiment_samples,
                                              metric = c("mean", "sum"),
                                              trim = 0,
                                              verbose = FALSE) {
  expr_matrix <- spatial_data[[data_type]]$data
  metadata <- spatial_data[[data_type]]$metadata
  metadata$sample_barcode <- paste(metadata$sample, metadata$barcode, sep = "_")
  
  selected_samples <- c(control_samples, experiment_samples)
  selected_barcodes <- metadata %>%
    filter(sample %in% selected_samples) %>%
    pull(sample_barcode)
  
  expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% selected_barcodes, drop = FALSE]
  
  col_sample_map <- metadata %>%
    filter(sample_barcode %in% colnames(expr_matrix)) %>%
    select(sample_barcode, sample)
  
  sample_levels <- unique(col_sample_map$sample)
  
  if (verbose) {
    cat(">>> Computing metrics:", paste(metric, collapse = ", "), "\n")
  }
  
  results_list <- list()
  
  for (m in metric) {
    if (verbose) {
      cat(">>", m, "\n")
      pb <- txtProgressBar(min = 0, max = length(sample_levels), style = 3)
    }
    
    metric_df <- data.frame(
      peak = spatial_data[[data_type]]$annotate$peak_id,
      gene = spatial_data[[data_type]]$annotate$gene_name,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(sample_levels)) {
      s <- sample_levels[i]
      barcodes <- col_sample_map %>%
        filter(sample == s) %>%
        pull(sample_barcode)
      
      mat <- expr_matrix[, barcodes, drop = FALSE]
      
      values <- switch(m,
                       sum        = rowSums(mat, na.rm = TRUE),
                       mean       = apply(mat, 1, function(x) mean(x, trim = trim, na.rm = TRUE)),
                       median     = apply(mat, 1, function(x) median(x, na.rm = TRUE)),
                       var        = apply(mat, 1, function(x) var(x, na.rm = TRUE)),
                       IQR        = apply(mat, 1, function(x) IQR(x, na.rm = TRUE)),
                       diff_range = apply(mat, 1, function(x) diff(range(x, na.rm = TRUE))),
                       skewness   = apply(mat, 1, function(x) e1071::skewness(x, na.rm = TRUE)),
                       kurtosis   = apply(mat, 1, function(x) e1071::kurtosis(x, na.rm = TRUE)),
                       stop("Unsupported metric: ", m)
      )
      
      metric_df[[s]] <- values
      
      if (verbose) setTxtProgressBar(pb, i)
    }
    
    if (verbose) close(pb)
    
    results_list[[m]] <- metric_df
  }
  
  return(results_list)
}

compute_global_expression_multi <- function(spatial_data,
                                            data_types,
                                            control_samples,
                                            experiment_samples,
                                            metrics = c("mean", "sum"),
                                            trim = 0,
                                            verbose = TRUE) {
  #' Wrapper to compute global expression summaries across multiple data types
  #'
  #' @param spatial_data List with assays in spatial format (must contain named entries matching `data_types`)
  #' @param data_types Character vector with names of data_type slots to iterate over
  #' @param control_samples Character vector of control sample IDs
  #' @param experiment_samples Character vector of experimental sample IDs
  #' @param metrics Character vector of metrics to compute (e.g., "mean", "sum")
  #' @param trim Trim fraction for computing trimmed means (default 0)
  #' @param verbose Print progress
  #'
  #' @return Named list with structure:
  #'   result[[data_type]][[metric]] -> dataframe of peak/gene vs sample
  
  result_list <- list()
  
  for (dt in data_types) {
    if (verbose) {
      cat("\n=====================================\n")
      cat(">>> Processing data_type:", dt, "\n")
      cat("=====================================\n")
    }
    
    summary_res <- compute_global_expression_summary(
      spatial_data       = spatial_data,
      data_type          = dt,
      control_samples    = control_samples,
      experiment_samples = experiment_samples,
      metric             = metrics,
      trim               = trim,
      verbose            = verbose
    )
    
    result_list[[dt]] <- summary_res
  }
  
  return(result_list)
}

generate_global_heatmap <- function(data,
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


ris2q29_summary_global_expression <- compute_global_expression_multi(
  spatial_data       = ris3q29_st_data,
  data_types         = c("raw_data", "quantile_normalize_resolution_0.1", 
                         "quantile_normalize_resolution_0.2",
                         "quantile_normalize_resolution_0.4"),
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  metrics            = c("mean", "sum"),
  trim               = 0
)






generate_global_heatmap(
  data               = ris2q29_summary_global_expression,
  genes              = genes_del3q29,
  summary_metric     = "mean",
  data_type          = "quantile_normalize_resolution_0.4",
  metadata           = sample_info,
  group1             = "mouse_genotype",
  group2             = "treatment",
  colors_group1      = c("wtwt" = "gray", "wtdel" = "blue"),
  colors_group2      = c("saline" = "orange", "risperidone" = "green"),
  title              = "quantile normalize; mean",
  scale_range        = c(-2, 2),
  remove_empty_rows  = TRUE,
  colors_on_heatmap  = c("#6a89b1ff", "white", "#832524ff"),
  z_scale            = TRUE,
  show_values        = FALSE,
  verbose            = TRUE
)
