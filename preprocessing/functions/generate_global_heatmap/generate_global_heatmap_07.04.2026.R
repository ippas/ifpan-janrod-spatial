generate_global_heatmap <- function(
  data,
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
  verbose = FALSE,
  aggregate_by_gene = FALSE,
  aggregate_genes = c("mean", "sum", "max"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_order = NULL,
  sample_order = NULL,
  show_gene_peak_mapping = TRUE,
  save_svg = NULL,
  svg_width = 8,
  svg_height = 6
) {
  
  aggregate_genes <- match.arg(aggregate_genes)
  
  if (verbose) message("[1] Selecting summary data...")
  summary_df <- data[[data_type]][[summary_metric]]
  
  if (verbose) message("[2] Filtering genes...")
  df_filtered <- summary_df %>%
    dplyr::filter(gene %in% genes)
  
  if (verbose) message("[3] Pivoting to long format...")
  df_long <- df_filtered %>%
    tidyr::pivot_longer(
      cols = -c(peak, gene),
      names_to = "sample_ID",
      values_to = "expression"
    ) %>%
    dplyr::left_join(metadata, by = "sample_ID")
  
  if (aggregate_by_gene) {
    
    if (verbose) {
      message(sprintf(
        "[4] Aggregating peaks to gene-level signal (%s)...",
        aggregate_genes
      ))
    }
    
    df_long <- df_long %>%
      dplyr::mutate(gene_symbol = sub("^.*-", "", peak))
    
    if (isTRUE(show_gene_peak_mapping)) {
      gene_peak_mapping_df <- df_long %>%
        dplyr::distinct(gene_symbol, peak_id = peak) %>%
        dplyr::arrange(gene_symbol, peak_id)
      
      message("[gene-peak mapping]")
      message(
        paste(
          capture.output(print(gene_peak_mapping_df, n = Inf)),
          collapse = "\n"
        )
      )
    }
    
    df_long <- df_long %>%
      dplyr::group_by(gene_symbol, sample_ID) %>%
      dplyr::summarise(
        expression = dplyr::case_when(
          aggregate_genes == "mean" ~ mean(expression, na.rm = TRUE),
          aggregate_genes == "sum"  ~ sum(expression, na.rm = TRUE),
          aggregate_genes == "max"  ~ max(expression, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      dplyr::rename(peak = gene_symbol) %>%
      dplyr::left_join(metadata, by = "sample_ID")
  }
  
  if (verbose) message("[5] Scaling values...")
  df_scaled <- df_long %>%
    dplyr::group_by(peak) %>%
    dplyr::mutate(
      value = if (isTRUE(z_scale)) as.numeric(scale(expression)[, 1]) else expression
    ) %>%
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
  
  if (verbose) message("[6] Preparing column annotations...")
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
      } else {
        .
      }
    }
  
  if (!is.null(sample_order)) {
    if (is.numeric(sample_order)) {
      if (any(sample_order < 1 | sample_order > nrow(meta2))) {
        stop("sample_order contains indices outside 1..n_samples (nrow(meta2)).")
      }
      meta2 <- meta2[sample_order, , drop = FALSE]
    } else {
      sample_order <- as.character(sample_order)
      
      missing_samples <- setdiff(sample_order, meta2$sample_ID)
      if (length(missing_samples) > 0) {
        stop(
          "sample_order contains sample_IDs not present in the data/metadata: ",
          paste(missing_samples, collapse = ", ")
        )
      }
      
      remaining <- setdiff(meta2$sample_ID, sample_order)
      final_order <- c(sample_order, remaining)
      meta2 <- meta2[match(final_order, meta2$sample_ID), , drop = FALSE]
    }
  }
  
  if (verbose) message("[7] Building heatmap matrix...")
  mat <- df_scaled %>%
    dplyr::select(peak, sample_ID, value) %>%
    tidyr::pivot_wider(
      names_from = sample_ID,
      values_from = value
    ) %>%
    tibble::column_to_rownames("peak") %>%
    as.matrix()
  
  mat <- mat[, meta2$sample_ID, drop = FALSE]
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message(">>> No data to plot. Exiting.")
    return(invisible(NULL))
  }
  
  if (isTRUE(z_scale)) {
    mat <- pmax(pmin(mat, scale_range[2]), scale_range[1])
  }
  
  row_order_idx <- NULL
  if (!isTRUE(cluster_rows) && !is.null(row_order)) {
    
    if (is.numeric(row_order)) {
      if (any(row_order < 1 | row_order > nrow(mat))) {
        stop("row_order contains indices outside 1..nrow(mat).")
      }
      row_order_idx <- row_order
    } else {
      row_order <- as.character(row_order)
      missing_rows <- setdiff(row_order, rownames(mat))
      if (length(missing_rows) > 0) {
        stop(
          "row_order contains rows not present in the matrix: ",
          paste(missing_rows, collapse = ", ")
        )
      }
      remaining <- setdiff(rownames(mat), row_order)
      final_order <- c(row_order, remaining)
      row_order_idx <- match(final_order, rownames(mat))
    }
  }
  
  if (verbose) message("[8] Preparing color mapping / legend...")
  
  if (isTRUE(z_scale)) {
    breaks <- seq(scale_range[1], scale_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Z-score"
    legend_at <- seq(scale_range[1], scale_range[2], by = 1)
  } else {
    data_range <- range(mat, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = length(colors_on_heatmap))
    legend_name <- "Expression"
    legend_at <- NULL
  }
  
  hm_args <- list(
    mat,
    name = legend_name,
    col = circlize::colorRamp2(breaks, colors_on_heatmap),
    na_col = "grey80",
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    column_names_side = "top",
    cluster_rows = isTRUE(cluster_rows),
    cluster_columns = isTRUE(cluster_columns),
    row_order = row_order_idx,
    column_title = title,
    heatmap_legend_param = list(
      direction = "horizontal",
      nrow = 1,
      at = legend_at
    )
  )
  
  if (isTRUE(show_values)) {
    hm_args$cell_fun <- function(j, i, x, y, width, height, fill) {
      val <- mat[i, j]
      lab <- if (is.na(val)) "NA" else sprintf("%.2f", val)
      grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 6))
    }
  }
  
  if (length(group_vars) > 0) {
    ann_df <- as.data.frame(meta2[dplyr::all_of(group_vars)])
    col_list <- list()
    
    if (!is.null(group1) && !is.null(colors_group1)) {
      col_list[[group1]] <- colors_group1
    }
    
    if (!is.null(group2) && !is.null(colors_group2)) {
      col_list[[group2]] <- colors_group2
    }
    
    hm_args$top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = ann_df,
      col = col_list,
      which = "column",
      annotation_name_side = "left"
    )
  }
  
  ht <- do.call(ComplexHeatmap::Heatmap, hm_args)
  
  draw_heatmap <- function() {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE
    )
  }
  
  if (!is.null(save_svg)) {
    if (verbose) message(sprintf("[9] Saving SVG -> %s", save_svg))
    grDevices::svg(filename = save_svg, width = svg_width, height = svg_height)
    draw_heatmap()
    grDevices::dev.off()
  } else {
    draw_heatmap()
  }
  
  invisible(ht)
}

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
  # title              = "quantile normalize; sum",
  scale_range        = c(-2, 2),
  remove_empty_rows  = TRUE,
  colors_on_heatmap  = c("navy", "white", "firebrick3"),
  z_scale            = TRUE,
  show_values        = FALSE,
  verbose            = TRUE,
  aggregate_by_gene  = TRUE,
  aggregate_genes    = "sum",
  cluster_rows       = FALSE,
  svg_width          = 10,
  svg_height         = 7,
  show_gene_peak_mapping = TRUE,
  row_order = c(
    "Bdh1", "Dlg1", "Pigz", "Senp5", "Pak2", "Pigx", "Cep19",
    "Nrros", "Fbxo45", "Wdr53", "Rnf168", "Tctex1d2", "Pcyt1a", "Tfrc"
  ),
  sample_order = c(
    "S13839Nr1", "S13839Nr12", "S13839Nr20", "S13839Nr22", "S13839Nr5",
    "S13839Nr10", "S13839Nr13", "S13839Nr16", "S13839Nr24", "S13839Nr25",
    "S13839Nr4", "S13839Nr11", "S13839Nr14", "S13839Nr18", "S13839Nr26",
    "S13839Nr6", "S13839Nr15", "S13839Nr17", "S13839Nr19", "S13839Nr2",
    "S13839Nr21", "S13839Nr23", "S13839Nr7", "S13839Nr8", "S13839Nr9"
  )
)


# ris2q29_summary_global_expression <- compute_global_expression_multi(
#   spatial_data       = ris3q29_st_data,
#   data_types         = c("raw_data", "quantile_normalize_resolution_0.1", 
#                          "quantile_normalize_resolution_0.2",
#                          "quantile_normalize_resolution_0.4"),
#   control_samples    = samples_wt,
#   experiment_samples = samples_del,
#   metrics            = c("mean", "sum"),
#   trim               = 0
# )

