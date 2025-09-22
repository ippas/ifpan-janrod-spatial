run_edger_twofactor_for_clusters_v2 <- function(
  raw_data_by_resolution,
  sample_meta,
  sample_col         = "sample_id",
  # factors as in your base function
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  # where counts are stored inside the structure
  groups             = c("control","experiment"),
  count_slot         = "sum",         # e.g. $control$sum / $experiment$sum
  clusters           = NULL,          # NULL = process all clusters
  robust             = TRUE,
  verbose            = TRUE,
  return_combined    = TRUE,          # add one combined table with 'cluster' column
  skip_failed        = TRUE,          # skip clusters that fail, instead of stopping
  ...                                # <---- NOWE: przekazywane dalej do edger_twofactor_anova_and_posthoc
){
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # Choose clusters
  all_clusters <- names(raw_data_by_resolution)
  if (is.null(all_clusters) || !length(all_clusters))
    stop("`raw_data_by_resolution` does not contain cluster names.")
  if (is.null(clusters)) clusters <- all_clusters
  clusters <- intersect(clusters, all_clusters)
  if (!length(clusters)) stop("None of the requested clusters found in data.")
  
  results <- vector("list", length(clusters))
  names(results) <- clusters
  combined_list <- list()
  
  for (cl in clusters){
    .msg("[cluster %s] Preparing count matrix...", cl)
    
    # 1) Extract count matrices for declared groups
    mats <- list()
    for (g in groups){
      node <- raw_data_by_resolution[[cl]][[g]]
      if (is.null(node))
        stop(sprintf("Cluster '%s' does not contain group '%s'.", cl, g))
      mat <- node[[count_slot]]
      if (is.null(mat) || !is.matrix(mat))
        stop(sprintf("In cluster '%s' group '%s', slot '%s' is not a matrix.", cl, g, count_slot))
      mats[[g]] <- mat
    }
    
    # 2) Align genes by rownames and cbind
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common_genes))
      stop(sprintf("Cluster '%s': no common genes across groups %s.",
                   cl, paste(groups, collapse=", ")))
    mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
    counts <- do.call(cbind, mats)
    
    if (is.null(colnames(counts)))
      stop(sprintf("Cluster '%s': counts matrix has no sample colnames.", cl))
    
    # 3) Align metadata
    smeta <- as.data.frame(sample_meta)
    if (any(duplicated(smeta[[sample_col]]))) {
      .msg("[cluster %s] Metadata: duplicates in '%s' — keeping first occurrence.", cl, sample_col)
      smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
    }
    smeta_idx <- match(colnames(counts), smeta[[sample_col]])
    if (anyNA(smeta_idx)) {
      miss <- colnames(counts)[is.na(smeta_idx)]
      stop(sprintf("[cluster %s] Missing metadata for samples: %s", cl, paste(miss, collapse=", ")))
    }
    smeta <- smeta[smeta_idx, , drop = FALSE]
    
    # 4) Run base ANOVA function (NOWA WERSJA Z KONTRASTAMI)
    .msg("[cluster %s] Running two-factor edgeR ANOVA...", cl)
    run_ok <- TRUE
    res_cl <- try(
      edger_twofactor_anova_and_posthoc_v2(
        counts            = counts,
        sample_meta       = smeta,
        sample_col        = sample_col,
        first_factor_col  = first_factor_col,
        second_factor_col = second_factor_col,
        first_levels      = first_levels,
        second_levels     = second_levels,
        robust            = robust,
        verbose           = verbose,
        ...               # <<< KLUCZOWE: przekazujemy np. main_effects="marginal"
      ),
      silent = TRUE
    )
    
    if (inherits(res_cl, "try-error")) {
      msg <- paste0(res_cl)
      if (skip_failed){
        .msg("[cluster %s] ERROR — skipped (skip_failed=TRUE). Message: %s", cl, msg)
        run_ok <- FALSE
      } else {
        stop(sprintf("[cluster %s] ERROR: %s", cl, msg))
      }
    }
    
    if (run_ok){
      results[[cl]] <- res_cl
      if (isTRUE(return_combined) && !is.null(res_cl$combined)) {
        tmp <- res_cl$combined
        tmp$cluster <- cl
        combined_list[[cl]] <- tmp
      }
      .msg("[cluster %s] Finished.", cl)
    } else {
      results[[cl]] <- NULL
    }
  }
  
  out <- list(
    results = results,
    meta = list(
      clusters_processed = clusters,
      groups             = groups,
      count_slot         = count_slot,
      first_factor_col   = first_factor_col,
      second_factor_col  = second_factor_col,
      robust             = robust
    )
  )
  if (isTRUE(return_combined) && length(combined_list)){
    # Pad missing columns for safe row-binding
    all_cols <- Reduce(union, lapply(combined_list, colnames))
    combined_padded <- lapply(combined_list, function(df){
      miss <- setdiff(all_cols, colnames(df))
      if (length(miss)) df[miss] <- NA
      df[, all_cols, drop = FALSE]
    })
    out$combined <- do.call(rbind, combined_padded)
    rownames(out$combined) <- NULL
  }
  
  out
}

run_edger_twofactor_for_clusters_v2 <- function(
  raw_data_by_resolution,
  sample_meta,
  sample_col         = "sample_id",
  # factors as in your base function
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  # where counts are stored inside the structure
  groups             = c("control","experiment"),
  count_slot         = "sum",         # e.g. $control$sum / $experiment$sum
  clusters           = NULL,          # NULL = process all clusters
  robust             = TRUE,
  verbose            = TRUE,
  return_combined    = TRUE,          # add one combined table with 'cluster' column
  skip_failed        = TRUE,          # skip clusters that fail, instead of stopping
  return_preprocessing_data = FALSE,  # <<< nowa opcja
  drop_cols_all_na        = TRUE,     # <<< nowa opcja
  replace_partial_na_with = 0,        # <<< nowa opcja
  ...
){
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # Choose clusters
  all_clusters <- names(raw_data_by_resolution)
  if (is.null(all_clusters) || !length(all_clusters))
    stop("`raw_data_by_resolution` does not contain cluster names.")
  if (is.null(clusters)) clusters <- all_clusters
  clusters <- intersect(clusters, all_clusters)
  if (!length(clusters)) stop("None of the requested clusters found in data.")
  
  results <- vector("list", length(clusters))
  names(results) <- clusters
  combined_list <- list()
  
  for (cl in clusters){
    .msg("[cluster %s] Preparing count matrix...", cl)
    
    # 1) Extract count matrices for declared groups
    mats <- list()
    for (g in groups){
      node <- raw_data_by_resolution[[cl]][[g]]
      if (is.null(node))
        stop(sprintf("Cluster '%s' does not contain group '%s'.", cl, g))
      mat <- node[[count_slot]]
      if (is.null(mat) || !is.matrix(mat))
        stop(sprintf("In cluster '%s' group '%s', slot '%s' is not a matrix.", cl, g, count_slot))
      mats[[g]] <- mat
    }
    
    # 2) Align genes by rownames and cbind
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common_genes))
      stop(sprintf("Cluster '%s': no common genes across groups %s.",
                   cl, paste(groups, collapse=", ")))
    mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
    counts <- do.call(cbind, mats)
    
    if (is.null(colnames(counts)))
      stop(sprintf("Cluster '%s': counts matrix has no sample colnames.", cl))
    
    # 3) Align metadata
    smeta <- as.data.frame(sample_meta)
    if (any(duplicated(smeta[[sample_col]]))) {
      .msg("[cluster %s] Metadata: duplicates in '%s' — keeping first occurrence.", cl, sample_col)
      smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
    }
    smeta_idx <- match(colnames(counts), smeta[[sample_col]])
    if (anyNA(smeta_idx)) {
      miss <- colnames(counts)[is.na(smeta_idx)]
      stop(sprintf("[cluster %s] Missing metadata for samples: %s", cl, paste(miss, collapse=", ")))
    }
    smeta <- smeta[smeta_idx, , drop = FALSE]
    
    # 3a) Usuwanie kolumn samych NA
    if (isTRUE(drop_cols_all_na)) {
      na_all_cols <- which(apply(counts, 2, function(x) all(is.na(x))))
      if (length(na_all_cols) > 0) {
        .msg("[cluster %s] Dropping %d columns with all NA: %s",
             cl, length(na_all_cols), paste(colnames(counts)[na_all_cols], collapse=", "))
        counts <- counts[, -na_all_cols, drop = FALSE]
        smeta  <- smeta[-na_all_cols, , drop = FALSE]
      }
    }
    
    # 3b) Zamiana pojedynczych NA
    if (anyNA(counts)) {
      if (is.null(replace_partial_na_with)) {
        stop(sprintf("[cluster %s] Counts matrix still contains NA. Set replace_partial_na_with.", cl))
      } else {
        n_na <- sum(is.na(counts))
        .msg("[cluster %s] Replacing %d NA with %s", cl, n_na, as.character(replace_partial_na_with))
        counts[is.na(counts)] <- replace_partial_na_with
      }
    }
    
    # sanity: negative values?
    if (any(counts < 0, na.rm = TRUE)) {
      stop(sprintf("[cluster %s] Negative counts not allowed.", cl))
    }
    
    # 3c) optional logCPM
    logCPM_mat <- NULL
    if (isTRUE(return_preprocessing_data)) {
      y <- edgeR::DGEList(counts = counts, samples = smeta)
      y <- edgeR::calcNormFactors(y)
      logCPM_mat <- edgeR::cpm(y, log = TRUE, prior.count = 2)
    }
    
    # 4) Run base ANOVA function
    .msg("[cluster %s] Running two-factor edgeR ANOVA...", cl)
    run_ok <- TRUE
    res_cl <- try(
      edger_twofactor_anova_and_posthoc_v2(
        counts            = counts,
        sample_meta       = smeta,
        sample_col        = sample_col,
        first_factor_col  = first_factor_col,
        second_factor_col = second_factor_col,
        first_levels      = first_levels,
        second_levels     = second_levels,
        robust            = robust,
        verbose           = verbose,
        ...
      ),
      silent = TRUE
    )
    
    if (inherits(res_cl, "try-error")) {
      msg <- paste0(res_cl)
      if (skip_failed){
        .msg("[cluster %s] ERROR — skipped (skip_failed=TRUE). Message: %s", cl, msg)
        run_ok <- FALSE
      } else {
        stop(sprintf("[cluster %s] ERROR: %s", cl, msg))
      }
    }
    
    if (run_ok){
      if (isTRUE(return_preprocessing_data)) {
        res_cl$counts            <- counts
        res_cl$sample_meta       <- smeta
        res_cl$logCPM_per_sample <- logCPM_mat
      }
      results[[cl]] <- res_cl
      
      if (isTRUE(return_combined) && !is.null(res_cl$combined)) {
        tmp <- res_cl$combined
        tmp$cluster <- cl
        combined_list[[cl]] <- tmp
      }
      .msg("[cluster %s] Finished.", cl)
    } else {
      results[[cl]] <- NULL
    }
  }
  
  out <- list(
    results = results,
    meta = list(
      clusters_processed = clusters,
      groups             = groups,
      count_slot         = count_slot,
      first_factor_col   = first_factor_col,
      second_factor_col  = second_factor_col,
      robust             = robust,
      return_preprocessing_data = return_preprocessing_data
    )
  )
  if (isTRUE(return_combined) && length(combined_list)){
    all_cols <- Reduce(union, lapply(combined_list, colnames))
    combined_padded <- lapply(combined_list, function(df){
      miss <- setdiff(all_cols, colnames(df))
      if (length(miss)) df[miss] <- NA
      df[, all_cols, drop = FALSE]
    })
    out$combined <- do.call(rbind, combined_padded)
    rownames(out$combined) <- NULL
  }
  
  out
}
