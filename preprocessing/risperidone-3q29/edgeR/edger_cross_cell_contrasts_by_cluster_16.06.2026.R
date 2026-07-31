# ##############################################################################
# edgeR cross-cell contrasts by cluster
#
# Date: 2026-06-16
#
# Purpose:
#   Calculate two cross-cell contrasts in 2x2 design:
#
#   1) risDel vs salWt
#      wtdel_risperidone - wtwt_saline
#
#   2) risWt vs salDel
#      wtwt_risperidone - wtdel_saline
#
# Output:
#   One combined data.frame:
#     df_cross_cell_contrasts
#
# No saving to file.
# No final filtering by FDR/PValue/logFC.
# ##############################################################################


# ==============================================================================
# 1. Libraries
# ==============================================================================

library(edgeR)
library(dplyr)
library(tibble)
library(purrr)


# ==============================================================================
# 2. Helper functions
# ==============================================================================

parse_ids <- function(rn_vec) {
  
  split_tokens <- strsplit(rn_vec, "[ -]+", perl = TRUE)
  
  peak_id <- vapply(split_tokens, function(tok) {
    if (length(tok) >= 2) {
      paste(tok[1:2], collapse = "-")
    } else {
      tok[1]
    }
  }, character(1))
  
  gene_symbol <- vapply(split_tokens, function(tok) {
    if (length(tok) >= 1) {
      tok[length(tok)]
    } else {
      NA_character_
    }
  }, character(1))
  
  data.frame(
    name_id     = rn_vec,
    peak_id     = peak_id,
    gene_symbol = gene_symbol,
    stringsAsFactors = FALSE
  )
}


add_ids <- function(df) {
  
  if (is.null(rownames(df))) {
    stop("Input table has no rownames. Cannot add name_id / peak_id / gene_symbol.")
  }
  
  cbind(
    parse_ids(rownames(df)),
    df,
    row.names = NULL
  )
}


rename_contrast_columns <- function(df, contrast_name) {
  
  old_names <- c(
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR"
  )
  
  new_names <- paste0(
    old_names,
    "_",
    contrast_name
  )
  
  names(new_names) <- old_names
  
  colnames(df) <- ifelse(
    colnames(df) %in% old_names,
    new_names[colnames(df)],
    colnames(df)
  )
  
  df
}


make_contrast_vector <- function(design, numerator_cell, denominator_cell) {
  
  if (!numerator_cell %in% colnames(design)) {
    stop("Numerator cell not found in design: ", numerator_cell)
  }
  
  if (!denominator_cell %in% colnames(design)) {
    stop("Denominator cell not found in design: ", denominator_cell)
  }
  
  contrast <- numeric(ncol(design))
  names(contrast) <- colnames(design)
  
  contrast[numerator_cell]   <-  1
  contrast[denominator_cell] <- -1
  
  contrast
}


run_qlf_contrast <- function(fit, contrast, contrast_name) {
  
  qlf <- edgeR::glmQLFTest(
    fit,
    contrast = contrast
  )
  
  tab <- edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "none"
  )$table
  
  tab <- add_ids(tab)
  
  tab <- as_tibble(tab)
  
  tab <- rename_contrast_columns(
    df = tab,
    contrast_name = contrast_name
  )
  
  tab
}


# ==============================================================================
# 3. Single-cluster function
# ==============================================================================

run_edger_two_cross_cell_contrasts <- function(
  counts,
  sample_meta,
  sample_col        = "sample_ID",
  first_factor_col  = "mouse_genotype",
  second_factor_col = "treatment",
  first_levels      = c("wtwt", "wtdel"),
  second_levels     = c("saline", "risperidone"),
  robust            = TRUE,
  verbose           = TRUE
) {
  
  msg <- function(...) {
    if (isTRUE(verbose)) {
      message(sprintf(...))
    }
  }
  
  counts <- as.matrix(counts)
  smeta  <- as.data.frame(sample_meta)
  
  # ---------------------------------------------------------------------------
  # Match metadata to count matrix columns
  # ---------------------------------------------------------------------------
  
  idx <- match(
    colnames(counts),
    smeta[[sample_col]]
  )
  
  if (anyNA(idx)) {
    
    missing_samples <- colnames(counts)[is.na(idx)]
    
    stop(
      "Missing metadata for samples: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  smeta <- smeta[idx, , drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # Drop all-NA samples only
  # ---------------------------------------------------------------------------
  
  na_cols <- which(
    apply(counts, 2, function(x) all(is.na(x)))
  )
  
  if (length(na_cols) > 0) {
    
    msg(
      "[preprocess] Dropping %d all-NA samples: %s",
      length(na_cols),
      paste(colnames(counts)[na_cols], collapse = ", ")
    )
    
    counts <- counts[, -na_cols, drop = FALSE]
    smeta  <- smeta[-na_cols, , drop = FALSE]
  }
  
  if (ncol(counts) == 0) {
    warning("All samples were dropped.")
    return(tibble())
  }
  
  # ---------------------------------------------------------------------------
  # Set factors
  # ---------------------------------------------------------------------------
  
  smeta[[first_factor_col]] <- factor(
    smeta[[first_factor_col]],
    levels = first_levels
  )
  
  smeta[[second_factor_col]] <- factor(
    smeta[[second_factor_col]],
    levels = second_levels
  )
  
  if (anyNA(smeta[[first_factor_col]])) {
    bad <- smeta[[sample_col]][is.na(smeta[[first_factor_col]])]
    
    stop(
      "Some samples have NA after setting first factor levels: ",
      paste(bad, collapse = ", ")
    )
  }
  
  if (anyNA(smeta[[second_factor_col]])) {
    bad <- smeta[[sample_col]][is.na(smeta[[second_factor_col]])]
    
    stop(
      "Some samples have NA after setting second factor levels: ",
      paste(bad, collapse = ", ")
    )
  }
  
  smeta$edgeR_group <- interaction(
    smeta[[first_factor_col]],
    smeta[[second_factor_col]],
    sep  = "_",
    drop = TRUE
  )
  
  msg(
    "[meta] Samples: %d | groups: %s",
    ncol(counts),
    paste(levels(smeta$edgeR_group), collapse = ", ")
  )
  
  # ---------------------------------------------------------------------------
  # edgeR object
  # ---------------------------------------------------------------------------
  
  dge <- edgeR::DGEList(
    counts  = counts,
    samples = smeta,
    group   = smeta$edgeR_group
  )
  
  # Standard edgeR low-expression filtering.
  # This is not final filtering by significance.
  keep <- edgeR::filterByExpr(
    dge,
    group = smeta$edgeR_group
  )
  
  n_before <- nrow(dge)
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  dge <- edgeR::calcNormFactors(dge)
  
  msg(
    "[edgeR] Kept %d/%d features after filterByExpr.",
    nrow(dge),
    n_before
  )
  
  # ---------------------------------------------------------------------------
  # Design: one coefficient per 2x2 cell
  # ---------------------------------------------------------------------------
  
  design <- stats::model.matrix(
    ~ 0 + edgeR_group,
    data = dge$samples
  )
  
  colnames(design) <- gsub(
    "^edgeR_group",
    "",
    colnames(design)
  )
  
  msg(
    "[design] Coefficients: %s",
    paste(colnames(design), collapse = ", ")
  )
  
  required_cells <- c(
    "wtwt_saline",
    "wtdel_saline",
    "wtwt_risperidone",
    "wtdel_risperidone"
  )
  
  missing_cells <- setdiff(
    required_cells,
    colnames(design)
  )
  
  if (length(missing_cells) > 0) {
    stop(
      "Missing required design cells: ",
      paste(missing_cells, collapse = ", ")
    )
  }
  
  # ---------------------------------------------------------------------------
  # Fit model once
  # ---------------------------------------------------------------------------
  
  dge <- edgeR::estimateDisp(
    dge,
    design
  )
  
  fit <- edgeR::glmQLFit(
    dge,
    design,
    robust = robust
  )
  
  # ---------------------------------------------------------------------------
  # Define contrasts
  # ---------------------------------------------------------------------------
  
  contrast_risDel_vs_salWt <- make_contrast_vector(
    design = design,
    numerator_cell = "wtdel_risperidone",
    denominator_cell = "wtwt_saline"
  )
  
  contrast_risWt_vs_salDel <- make_contrast_vector(
    design = design,
    numerator_cell = "wtwt_risperidone",
    denominator_cell = "wtdel_saline"
  )
  
  msg("[contrast] risDel_vs_salWt: wtdel_risperidone - wtwt_saline")
  msg("[contrast] risWt_vs_salDel: wtwt_risperidone - wtdel_saline")
  
  # ---------------------------------------------------------------------------
  # Run tests
  # ---------------------------------------------------------------------------
  
  tab_risDel_vs_salWt <- run_qlf_contrast(
    fit = fit,
    contrast = contrast_risDel_vs_salWt,
    contrast_name = "risDel_vs_salWt"
  )
  
  tab_risWt_vs_salDel <- run_qlf_contrast(
    fit = fit,
    contrast = contrast_risWt_vs_salDel,
    contrast_name = "risWt_vs_salDel"
  )
  
  # ---------------------------------------------------------------------------
  # Join both contrasts
  # ---------------------------------------------------------------------------
  
  result <- tab_risDel_vs_salWt %>%
    full_join(
      tab_risWt_vs_salDel,
      by = c(
        "name_id",
        "peak_id",
        "gene_symbol"
      )
    ) %>%
    mutate(
      contrast_risDel_vs_salWt = "wtdel_risperidone - wtwt_saline",
      contrast_risWt_vs_salDel = "wtwt_risperidone - wtdel_saline"
    ) %>%
    select(
      name_id,
      peak_id,
      gene_symbol,
      
      logFC_risDel_vs_salWt,
      logCPM_risDel_vs_salWt,
      F_risDel_vs_salWt,
      PValue_risDel_vs_salWt,
      FDR_risDel_vs_salWt,
      
      logFC_risWt_vs_salDel,
      logCPM_risWt_vs_salDel,
      F_risWt_vs_salDel,
      PValue_risWt_vs_salDel,
      FDR_risWt_vs_salDel,
      
      contrast_risDel_vs_salWt,
      contrast_risWt_vs_salDel
    )
  
  result
}


# ==============================================================================
# 4. Multi-cluster function
# ==============================================================================

run_edger_two_cross_cell_contrasts_for_clusters <- function(
  raw_data_by_resolution,
  sample_meta,
  sample_col        = "sample_ID",
  first_factor_col  = "mouse_genotype",
  second_factor_col = "treatment",
  first_levels      = c("wtwt", "wtdel"),
  second_levels     = c("saline", "risperidone"),
  groups            = c("control", "experiment"),
  count_slot        = "sum",
  clusters          = NULL,
  robust            = TRUE,
  verbose           = TRUE
) {
  
  if (is.null(clusters)) {
    clusters <- names(raw_data_by_resolution)
  }
  
  result_list <- list()
  
  for (cl in clusters) {
    
    if (isTRUE(verbose)) {
      message("\n============================================================")
      message("[cluster] ", cl)
      message("============================================================")
    }
    
    cluster_obj <- raw_data_by_resolution[[cl]]
    
    if (is.null(cluster_obj)) {
      warning("Skipping missing cluster: ", cl)
      next
    }
    
    missing_groups <- setdiff(
      groups,
      names(cluster_obj)
    )
    
    if (length(missing_groups) > 0) {
      warning(
        "Skipping ",
        cl,
        " because groups are missing: ",
        paste(missing_groups, collapse = ", ")
      )
      
      next
    }
    
    count_list <- lapply(groups, function(g) {
      
      x <- cluster_obj[[g]][[count_slot]]
      
      if (is.null(x)) {
        stop(
          "Missing count slot: ",
          cl,
          " / ",
          g,
          " / ",
          count_slot
        )
      }
      
      as.matrix(x)
    })
    
    counts <- do.call(
      cbind,
      count_list
    )
    
    if (anyDuplicated(colnames(counts))) {
      
      warning(
        "Duplicated sample names in ",
        cl,
        ". Keeping first occurrence."
      )
      
      counts <- counts[, !duplicated(colnames(counts)), drop = FALSE]
    }
    
    cluster_result <- run_edger_two_cross_cell_contrasts(
      counts            = counts,
      sample_meta       = sample_meta,
      sample_col        = sample_col,
      first_factor_col  = first_factor_col,
      second_factor_col = second_factor_col,
      first_levels      = first_levels,
      second_levels     = second_levels,
      robust            = robust,
      verbose           = verbose
    )
    
    cluster_result <- cluster_result %>%
      mutate(
        cluster = cl
      ) %>%
      select(
        cluster,
        everything()
      )
    
    result_list[[cl]] <- cluster_result
  }
  
  result_df <- bind_rows(result_list)
  
  result_df
}


# ==============================================================================
# 5. Run analysis
# ==============================================================================

df_cross_cell_contrasts <- run_edger_two_cross_cell_contrasts_for_clusters(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            = metadata_ris3q29,
  sample_col             = "sample_ID",
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  first_levels           = c("wtwt", "wtdel"),
  second_levels          = c("saline", "risperidone"),
  groups                 = c("control", "experiment"),
  count_slot             = "sum",
  clusters               = NULL,
  robust                 = TRUE,
  verbose                = TRUE
)


# ==============================================================================
# 6. Quick checks
# ==============================================================================

df_cross_cell_contrasts %>% dim
  head() %>% as.data.frame()


df_cross_cell_contrasts %>%
  filter(gene_symbol == "Arc")


df_cross_cell_contrasts %>%
  count(cluster)


df_cross_cell_contrasts %>%
  summarise(
    n_rows     = n(),
    n_clusters = n_distinct(cluster),
    n_genes    = n_distinct(gene_symbol),
    n_features = n_distinct(name_id)
  )


####
df_edger_marginal_proportional <- readxl::read_excel("/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_withMeans_09.06.2026.xlsx", 
                                    sheet = 1
                                  )

df_edger_marginal_proportional %>% dim


# ==============================================================================
# 1. Prepare new cross-cell contrasts for joining
# ==============================================================================

df_cross_cell_contrasts_for_join <- df_cross_cell_contrasts %>%
  select(
    cluster,
    name_id,
    peak_id,
    gene_symbol,
    
    logFC_risDel_vs_salWt,
    logCPM_risDel_vs_salWt,
    F_risDel_vs_salWt,
    PValue_risDel_vs_salWt,
    FDR_risDel_vs_salWt,
    
    logFC_risWt_vs_salDel,
    logCPM_risWt_vs_salDel,
    F_risWt_vs_salDel,
    PValue_risWt_vs_salDel,
    FDR_risWt_vs_salDel
  )


# ==============================================================================
# 2. Join with the previously generated edgeR summary table
# ==============================================================================

df_edger_marginal_proportional_with_cross_contrasts <- df_edger_marginal_proportional %>%
  left_join(
    df_cross_cell_contrasts_for_join,
    by = c(
      "cluster",
      "name_id",
      "peak_id",
      "gene_symbol"
    )
  )


# ==============================================================================
# 3. Move newly added contrast columns before mean_* columns
# ==============================================================================

added_cols <- c(
  "logFC_risDel_vs_salWt",
  "logCPM_risDel_vs_salWt",
  "F_risDel_vs_salWt",
  "PValue_risDel_vs_salWt",
  "FDR_risDel_vs_salWt",
  
  "logFC_risWt_vs_salDel",
  "logCPM_risWt_vs_salDel",
  "F_risWt_vs_salDel",
  "PValue_risWt_vs_salDel",
  "FDR_risWt_vs_salDel"
)

mean_cols <- grep(
  "^mean_",
  colnames(df_edger_marginal_proportional_with_cross_contrasts),
  value = TRUE
)

non_added_non_mean_cols <- setdiff(
  colnames(df_edger_marginal_proportional_with_cross_contrasts),
  c(added_cols, mean_cols)
)

df_edger_marginal_proportional_with_cross_contrasts <- df_edger_marginal_proportional_with_cross_contrasts %>%
  select(
    all_of(non_added_non_mean_cols),
    all_of(added_cols),
    all_of(mean_cols)
  )


# ==============================================================================
# 4. Quick checks
# ==============================================================================

df_edger_marginal_proportional_with_cross_contrasts %>%
  head() %>%
  as.data.frame()

df_edger_marginal_proportional_with_cross_contrasts %>%
  summarise(
    n_rows = n(),
    n_cols = ncol(.),
    n_missing_risDel_vs_salWt = sum(is.na(logFC_risDel_vs_salWt)),
    n_missing_risWt_vs_salDel = sum(is.na(logFC_risWt_vs_salDel))
  )


# ==============================================================================
# 5. Save updated table
# ==============================================================================

output_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_withMeans_16.06.2026.xlsx"

write.xlsx(
  df_edger_marginal_proportional_with_cross_contrasts,
  file = output_file,
  overwrite = TRUE
)


# ==============================================================================
# 6. Clean environment
# ==============================================================================

rm(
  list = intersect(
    c(
      "df_cross_cell_contrasts_for_join",
      "df_edger_marginal_proportional_with_cross_contrasts",
      "added_cols",
      "mean_cols",
      "non_added_non_mean_cols",
      "output_file"
    ),
    ls()
  )
)
