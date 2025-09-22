prepare_logCPM_per_sample_for_clusters <- function(
  raw_data_by_resolution,
  sample_meta,
  sample_col              = "sample_id",
  groups                  = c("control","experiment"),
  count_slot              = "sum",
  clusters                = NULL,          # NULL = wszystkie
  sample_keep             = NULL,          # wektor ID próbek do zachowania (opcjonalnie)
  first_factor_col        = "genotype",    # używane tylko przy filterByExpr
  second_factor_col       = "treatment",
  normalized.lib.sizes    = TRUE,
  prior.count             = 2,
  filter_by_expr          = FALSE,         # TRUE = edgeR::filterByExpr
  return_format           = c("list","long"),
  zscore_rows             = FALSE,         # skala Z po wierszach do wizualizacji
  drop_cols_all_na        = TRUE,          # usuń kolumny, w których wszystkie wartości to NA
  replace_partial_na_with = 0,             # czym zastąpić pojedyncze NA (NULL = nie ruszaj)
  verbose                 = TRUE
){
  stopifnot(requireNamespace("edgeR", quietly = TRUE))
  stopifnot(is.list(raw_data_by_resolution))
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  return_format <- match.arg(return_format)
  
  all_clusters <- names(raw_data_by_resolution)
  if (is.null(all_clusters) || !length(all_clusters))
    stop("`raw_data_by_resolution` nie zawiera nazw klastrów.")
  if (is.null(clusters)) clusters <- all_clusters
  clusters <- intersect(clusters, all_clusters)
  if (!length(clusters)) stop("Nie znaleziono żądanych klastrów w 'raw_data_by_resolution'.")
  
  res_list <- vector("list", length(clusters)); names(res_list) <- clusters
  long_list <- list()
  
  # pomocnicza: melt + merge z metadanymi
  .melt_long <- function(mat, cluster, sample_meta, sample_col){
    df_long <- as.data.frame(mat, check.names = FALSE)
    df_long$name_id <- rownames(mat)
    parts <- strsplit(df_long$name_id, "-")
    df_long$peak_id <- vapply(parts, function(x) if (length(x) >= 2) paste(x[1:2], collapse="-") else NA_character_, character(1))
    df_long$gene_symbol <- vapply(parts, function(x) if (length(x)) x[length(x)] else NA_character_, character(1))
    df_long$cluster <- cluster
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("Package 'reshape2' is required for return_format='long'.")
    
    long <- reshape2::melt(
      df_long,
      id.vars = c("cluster","name_id","peak_id","gene_symbol"),
      variable.name = "sample_id",
      value.name = "logCPM",
      stringsAsFactors = FALSE
    )
    
    # dołącz metadane (np. treatment, mouse_genotype)
    long <- merge(long, sample_meta, by.x = "sample_id", by.y = sample_col, all.x = TRUE)
    long
  }
  
  for (cl in clusters){
    .msg("[cluster %s] Buduję macierz zliczeń…", cl)
    
    # 1) zbierz macierze dla wskazanych grup
    mats <- list()
    for (g in groups){
      node <- raw_data_by_resolution[[cl]][[g]]
      if (is.null(node)) stop(sprintf("Brak grupy '%s' w klastrze '%s'.", g, cl))
      mat <- node[[count_slot]]
      if (is.null(mat) || !is.matrix(mat))
        stop(sprintf("Slot '%s' w grupie '%s' (klaster %s) nie jest macierzą.", count_slot, g, cl))
      mats[[g]] <- mat
    }
    
    # 2) wspólne geny + cbind
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common_genes))
      stop(sprintf("Klaster %s: brak wspólnych genów pomiędzy grupami.", cl))
    mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
    counts <- do.call(cbind, mats)
    
    # (opcjonalnie) filtr próbek po ID
    if (!is.null(sample_keep)) {
      miss <- setdiff(sample_keep, colnames(counts))
      if (length(miss)) warning(sprintf("[cluster %s] sample_keep nie znaleziono: %s", cl, paste(miss, collapse=", ")))
      keep_cols <- intersect(sample_keep, colnames(counts))
      counts <- counts[, keep_cols, drop = FALSE]
    }
    
    if (is.null(colnames(counts)))
      stop(sprintf("Klaster %s: macierz zliczeń bez nazw kolumn (próbek).", cl))
    
    # 3) dopasuj metadane
    smeta <- as.data.frame(sample_meta)
    if (!sample_col %in% names(smeta))
      stop(sprintf("Kolumna '%s' nie istnieje w sample_meta.", sample_col))
    if (any(duplicated(smeta[[sample_col]]))) {
      .msg("[cluster %s] Duplikaty w '%s' — zachowuję pierwsze wystąpienie.", cl, sample_col)
      smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
    }
    idx <- match(colnames(counts), smeta[[sample_col]])
    if (anyNA(idx)) {
      miss <- colnames(counts)[is.na(idx)]
      stop(sprintf("[cluster %s] Brak metadanych dla próbek: %s", cl, paste(miss, collapse=", ")))
    }
    smeta <- smeta[idx, , drop = FALSE]
    
    # sanity: ujemne wartości?
    if (any(counts < 0, na.rm = TRUE)) {
      stop(sprintf("[cluster %s] Wykryto ujemne zliczenia — edgeR tego nie akceptuje.", cl))
    }
    
    # usuń kolumny z samymi NA
    if (isTRUE(drop_cols_all_na)) {
      na_all_cols <- which(apply(counts, 2, function(x) all(is.na(x))))
      if (length(na_all_cols) > 0) {
        .msg("[cluster %s] Usuwam %d kolumn z samymi NA: %s",
             cl, length(na_all_cols), paste(colnames(counts)[na_all_cols], collapse = ", "))
        counts <- counts[, -na_all_cols, drop = FALSE]
        smeta  <- smeta[-na_all_cols, , drop = FALSE]
      }
    }
    
    # zamiana pojedynczych NA
    if (anyNA(counts)) {
      if (is.null(replace_partial_na_with)) {
        stop(sprintf("[cluster %s] W macierzy zliczeń są NA. Ustaw replace_partial_na_with.", cl))
      } else {
        n_na <- sum(is.na(counts))
        .msg("[cluster %s] Zastępuję %d pojedynczych NA wartością %s.", cl, n_na, as.character(replace_partial_na_with))
        counts[is.na(counts)] <- replace_partial_na_with
      }
    }
    
    if (ncol(counts) == 0) {
      .msg("[cluster %s] Po czyszczeniu brak próbek. Pomijam ten klaster.", cl)
      res_list[[cl]] <- NULL
      next
    }
    
    # 4) edgeR → DGEList + (opcjonalnie) filterByExpr + TMM + logCPM
    y <- edgeR::DGEList(counts = counts, samples = smeta)
    
    if (isTRUE(filter_by_expr) &&
        all(c(first_factor_col, second_factor_col) %in% colnames(smeta))) {
      grp <- interaction(smeta[[first_factor_col]], smeta[[second_factor_col]], drop = TRUE)
      keep <- edgeR::filterByExpr(y, group = grp)
      kept_n <- sum(keep)
      .msg("[cluster %s] filterByExpr: %d/%d genów zachowanych.", cl, kept_n, nrow(y))
      if (kept_n == 0) {
        .msg("[cluster %s] Po filterByExpr brak genów. Pomijam ten klaster.", cl)
        res_list[[cl]] <- NULL
        next
      }
      y <- y[keep, , keep.lib.sizes = FALSE]
    }
    
    y <- edgeR::calcNormFactors(y)
    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = prior.count, normalized.lib.sizes = normalized.lib.sizes)
    
    # 5) (opcjonalnie) Z-score po wierszach
    if (isTRUE(zscore_rows) && nrow(logcpm) > 1) {
      logcpm <- t(scale(t(logcpm)))
    }
    
    res_list[[cl]] <- logcpm
    
    if (return_format == "long") {
      long_list[[cl]] <- .melt_long(logcpm, cl, smeta, sample_col)
    }
  }
  
  if (return_format == "list") {
    return(list(logCPM_by_cluster = res_list))
  } else {
    combined_long <- do.call(rbind, long_list)
    rownames(combined_long) <- NULL
    return(combined_long)
  }
}


prepare_logCPM_per_sample_for_clusters <- function(
  raw_data_by_resolution,
  sample_meta,
  sample_col              = "sample_id",
  groups                  = c("control","experiment"),
  count_slot              = "sum",
  clusters                = NULL,
  sample_keep             = NULL,
  first_factor_col        = "genotype",
  second_factor_col       = "treatment",
  normalized.lib.sizes    = TRUE,
  prior.count             = 2,
  filter_by_expr          = FALSE,
  return_format           = c("list","long"),
  zscore_rows             = FALSE,
  drop_cols_all_na        = TRUE,
  replace_partial_na_with = 0,
  verbose                 = TRUE
){
  stopifnot(requireNamespace("edgeR", quietly = TRUE))
  stopifnot(is.list(raw_data_by_resolution))
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  return_format <- match.arg(return_format)
  
  all_clusters <- names(raw_data_by_resolution)
  if (is.null(all_clusters) || !length(all_clusters))
    stop("`raw_data_by_resolution` nie zawiera nazw klastrów.")
  if (is.null(clusters)) clusters <- all_clusters
  clusters <- intersect(clusters, all_clusters)
  if (!length(clusters)) stop("Nie znaleziono żądanych klastrów w 'raw_data_by_resolution'.")
  
  res_list <- vector("list", length(clusters)); names(res_list) <- clusters
  long_list <- list()
  
  # melt + merge z metadanymi + dodanie kolumny logCPM_avg
  .melt_long <- function(mat, cluster, sample_meta, sample_col, logCPM_avg_vec){
    df_long <- as.data.frame(mat, check.names = FALSE)
    df_long$name_id <- rownames(mat)
    parts <- strsplit(df_long$name_id, "-")
    df_long$peak_id <- vapply(parts, function(x) if (length(x) >= 2) paste(x[1:2], collapse="-") else NA_character_, character(1))
    df_long$gene_symbol <- vapply(parts, function(x) if (length(x)) x[length(x)] else NA_character_, character(1))
    df_long$cluster <- cluster
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("Package 'reshape2' is required for return_format='long'.")
    
    long <- reshape2::melt(
      df_long,
      id.vars = c("cluster","name_id","peak_id","gene_symbol"),
      variable.name = "sample_id",
      value.name = "logCPM",
      stringsAsFactors = FALSE
    )
    
    # dołącz metadane (np. treatment, genotype)
    long <- merge(long, sample_meta, by.x = "sample_id", by.y = sample_col, all.x = TRUE)
    
    # dołącz kolumnę logCPM_avg (edgeR-style average logCPM)
    long$logCPM_avg <- logCPM_avg_vec[match(long$name_id, names(logCPM_avg_vec))]
    long
  }
  
  for (cl in clusters){
    .msg("[cluster %s] Buduję macierz zliczeń…", cl)
    
    # 1) zbierz macierze dla wskazanych grup
    mats <- list()
    for (g in groups){
      node <- raw_data_by_resolution[[cl]][[g]]
      if (is.null(node)) stop(sprintf("Brak grupy '%s' w klastrze '%s'.", g, cl))
      mat <- node[[count_slot]]
      if (is.null(mat) || !is.matrix(mat))
        stop(sprintf("Slot '%s' w grupie '%s' (klaster %s) nie jest macierzą.", count_slot, g, cl))
      mats[[g]] <- mat
    }
    
    # 2) wspólne geny + cbind
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common_genes))
      stop(sprintf("Klaster %s: brak wspólnych genów pomiędzy grupami.", cl))
    mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
    counts <- do.call(cbind, mats)
    
    # (opcjonalnie) filtr próbek po ID
    if (!is.null(sample_keep)) {
      miss <- setdiff(sample_keep, colnames(counts))
      if (length(miss)) warning(sprintf("[cluster %s] sample_keep nie znaleziono: %s", cl, paste(miss, collapse=", ")))
      keep_cols <- intersect(sample_keep, colnames(counts))
      counts <- counts[, keep_cols, drop = FALSE]
    }
    
    if (is.null(colnames(counts)))
      stop(sprintf("Klaster %s: macierz zliczeń bez nazw kolumn (próbek).", cl))
    
    # 3) dopasuj metadane
    smeta <- as.data.frame(sample_meta)
    if (!sample_col %in% names(smeta))
      stop(sprintf("Kolumna '%s' nie istnieje w sample_meta.", sample_col))
    if (any(duplicated(smeta[[sample_col]]))) {
      .msg("[cluster %s] Duplikaty w '%s' — zachowuję pierwsze wystąpienie.", cl, sample_col)
      smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
    }
    idx <- match(colnames(counts), smeta[[sample_col]])
    if (anyNA(idx)) {
      miss <- colnames(counts)[is.na(idx)]
      stop(sprintf("[cluster %s] Brak metadanych dla próbek: %s", cl, paste(miss, collapse=", ")))
    }
    smeta <- smeta[idx, , drop = FALSE]
    
    # sanity: ujemne wartości?
    if (any(counts < 0, na.rm = TRUE)) {
      stop(sprintf("[cluster %s] Wykryto ujemne zliczenia — edgeR tego nie akceptuje.", cl))
    }
    
    # usuń kolumny z samymi NA
    if (isTRUE(drop_cols_all_na)) {
      na_all_cols <- which(apply(counts, 2, function(x) all(is.na(x))))
      if (length(na_all_cols) > 0) {
        .msg("[cluster %s] Usuwam %d kolumn z samymi NA: %s",
             cl, length(na_all_cols), paste(colnames(counts)[na_all_cols], collapse = ", "))
        counts <- counts[, -na_all_cols, drop = FALSE]
        smeta  <- smeta[-na_all_cols, , drop = FALSE]
      }
    }
    
    # zamiana pojedynczych NA
    if (anyNA(counts)) {
      if (is.null(replace_partial_na_with)) {
        stop(sprintf("[cluster %s] W macierzy zliczeń są NA. Ustaw replace_partial_na_with.", cl))
      } else {
        n_na <- sum(is.na(counts))
        .msg("[cluster %s] Zastępuję %d pojedynczych NA wartością %s.", cl, n_na, as.character(replace_partial_na_with))
        counts[is.na(counts)] <- replace_partial_na_with
      }
    }
    
    if (ncol(counts) == 0) {
      .msg("[cluster %s] Po czyszczeniu brak próbek. Pomijam ten klaster.", cl)
      res_list[[cl]] <- NULL
      next
    }
    
    # 4) edgeR → DGEList + (opcjonalnie) filterByExpr + TMM + logCPM
    y <- edgeR::DGEList(counts = counts, samples = smeta)
    y <- edgeR::calcNormFactors(y)
    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = prior.count, normalized.lib.sizes = normalized.lib.sizes)
    
    # policz "średnie" logCPM (edgeR-style)
    eff_libsizes <- y$samples$lib.size * y$samples$norm.factors
    mean_lib <- mean(eff_libsizes)
    sum_counts <- rowSums(y$counts)
    cpm_all <- (sum_counts / mean_lib) * 1e6
    logCPM_avg_vec <- log2(cpm_all + prior.count)
    
    # 5) (opcjonalnie) Z-score po wierszach
    if (isTRUE(zscore_rows) && nrow(logcpm) > 1) {
      logcpm <- t(scale(t(logcpm)))
    }
    
    res_list[[cl]] <- logcpm
    
    if (return_format == "long") {
      long_list[[cl]] <- .melt_long(logcpm, cl, smeta, sample_col, logCPM_avg_vec)
    }
  }
  
  if (return_format == "list") {
    return(list(logCPM_by_cluster = res_list))
  } else {
    combined_long <- do.call(rbind, long_list)
    rownames(combined_long) <- NULL
    return(combined_long)
  }
}
