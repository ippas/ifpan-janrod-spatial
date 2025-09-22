edger_twofactor_anova_and_posthoc <- function(
  counts,                    # matrix: genes x samples (rownames = genes, colnames = samples)
  sample_meta,               # data.frame: sample_col, first_factor_col, second_factor_col
  sample_col         = "sample_id",
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  robust             = TRUE,
  verbose            = TRUE,
  gatekeep           = FALSE,          # if TRUE → create posthoc_gate with new FDR
  gatekeep_fdr       = 0.05,
  gatekeep_log2FC    = 0,
  p_adjust_method    = "BH"
){
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---- Helpers ---------------------------------------------------------------
  .add_ids <- function(df){
    rn <- rownames(df)
    ensembl_id  <- sub("^([^_-]+)[_-].*$", "\\1", rn)
    gene_symbol <- ifelse(grepl("^[^_-]+[_-]", rn),
                          sub("^[^_-]+[_-](.*)$", "\\1", rn),
                          NA_character_)
    cbind(
      data.frame(
        name_id     = rn,
        ensembl_id  = ensembl_id,
        gene_symbol = gene_symbol,
        stringsAsFactors = FALSE
      ),
      df,
      row.names = NULL
    )
  }
  .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
  .rename_block <- function(df, block, id_cols){
    x <- .pick(df, c(id_cols, "logFC","logCPM","F","PValue","FDR"))
    colnames(x) <- c(id_cols,
                     paste0(block,"_logFC"),
                     paste0(block,"_logCPM"),
                     paste0(block,"_F"),
                     paste0(block,"_PValue"),
                     paste0(block,"_FDR"))
    x
  }
  
  # ---- Metadata --------------------------------------------------------------
  .msg("[meta] Aligning metadata to count matrix...")
  smeta <- as.data.frame(sample_meta)
  idx <- match(colnames(counts), smeta[[sample_col]])
  if (anyNA(idx)) stop("Missing metadata for samples.")
  smeta <- smeta[idx, , drop = FALSE]
  
  if (is.null(first_levels))  first_levels  <- unique(as.character(smeta[[first_factor_col]]))
  if (is.null(second_levels)) second_levels <- unique(as.character(smeta[[second_factor_col]]))
  smeta[[first_factor_col]]  <- factor(smeta[[first_factor_col]],  levels = first_levels)
  smeta[[second_factor_col]] <- factor(smeta[[second_factor_col]], levels = second_levels)
  f1 <- first_factor_col; f2 <- second_factor_col
  .msg("[meta] Samples: %d | First '%s': %s | Second '%s': %s",
       ncol(counts), f1, paste(levels(smeta[[f1]]), collapse = ", "),
       f2, paste(levels(smeta[[f2]]), collapse = ", "))
  
  # ---- edgeR preprocessing ---------------------------------------------------
  .msg("[edgeR] DGEList, filterByExpr, TMM normalization...")
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = interaction(dge$samples[[f1]], dge$samples[[f2]], drop = TRUE))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  .msg("[edgeR] Kept %d/%d genes after filtering.", nrow(dge), nrow(counts))
  
  # ---- A) ANOVA-like model ---------------------------------------------------
  .msg("[model] Fitting ANOVA-like model: ~ %s * %s ...", f1, f2)
  designA <- stats::model.matrix(stats::as.formula(paste0("~ ", f1, " * ", f2)), data = dge$samples)
  dgeA    <- edgeR::estimateDisp(dge, designA)
  fitA    <- edgeR::glmQLFit(dgeA, designA, robust = robust)
  
  .msg("[test] Omnibus QL F-test (all non-intercept coefficients)...")
  q_global <- edgeR::glmQLFTest(fitA, coef = 2:ncol(designA))
  gtab_raw <- edgeR::topTags(q_global, n = Inf)$table
  tab_global <- .add_ids(gtab_raw)
  
  .msg("[test] Main effect '%s'...", f1)
  f1_idx <- grep(paste0("^", f1), colnames(designA))
  q_f1   <- edgeR::glmQLFTest(fitA, coef = f1_idx)
  tab_f1 <- .add_ids(edgeR::topTags(q_f1, n = Inf)$table)
  
  .msg("[test] Main effect '%s'...", f2)
  f2_idx <- grep(paste0("^", f2), colnames(designA))
  q_f2   <- edgeR::glmQLFTest(fitA, coef = f2_idx)
  tab_f2 <- .add_ids(edgeR::topTags(q_f2, n = Inf)$table)
  
  .msg("[test] Interaction '%s:%s'...", f1, f2)
  inter_idx <- grep(":", colnames(designA))
  q_int   <- edgeR::glmQLFTest(fitA, coef = inter_idx)
  tab_int <- .add_ids(edgeR::topTags(q_int, n = Inf)$table)
  
  # ---- B) Cell-means model ---------------------------------------------------
  .msg("[model] Fitting cell-means model for simple effects: ~ 0 + %s:%s ...", f1, f2)
  designM <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", f1, ":", f2)), data = dge$samples)
  dgeM    <- edgeR::estimateDisp(dge, designM)
  fitM    <- edgeR::glmQLFit(dgeM, designM, robust = robust)
  colsM   <- colnames(designM)
  
  cell_name  <- function(A, B) paste0(f1, A, ":", f2, B)
  which_cell <- function(A, B) match(cell_name(A, B), colsM)
  
  build_simple_effects <- function(){
    out <- list()
    levA <- levels(smeta[[f1]]); levB <- levels(smeta[[f2]])
    # SECOND within FIRST
    for (a in levA){
      for (i in seq_len(length(levB)-1)){
        for (j in (i+1):length(levB)){
          b1 <- levB[i]; b2 <- levB[j]
          out[[sprintf("SECOND_%s-%s|FIRST_%s", b2, b1, a)]] <- c(which_cell(a,b2), which_cell(a,b1))
        }
      }
    }
    # FIRST within SECOND
    for (b in levB){
      for (i in seq_len(length(levA)-1)){
        for (j in (i+1):length(levA)){
          a1 <- levA[i]; a2 <- levA[j]
          out[[sprintf("FIRST_%s-%s|SECOND_%s", a2, a1, b)]] <- c(which_cell(a2,b), which_cell(a1,b))
        }
      }
    }
    out
  }
  contrasts_idx <- build_simple_effects()
  
  run_contrast <- function(ix){
    contr <- numeric(ncol(designM))
    contr[ix[1]] <-  1
    contr[ix[2]] <- -1
    edgeR::glmQLFTest(fitM, contrast = contr)
  }
  
  .msg("[posthoc] Computing %d simple-effect contrasts (raw edgeR)...", length(contrasts_idx))
  posthoc_raw <- lapply(contrasts_idx, function(ix) edgeR::topTags(run_contrast(ix), n = Inf)$table)
  names(posthoc_raw) <- names(contrasts_idx)
  
  # ---- C) Gatekeeping --------------------------------------------------------
  posthoc_gate <- NULL
  if (isTRUE(gatekeep)) {
    .msg("[gatekeep] Applying thresholds: FDR <= %.3f ; |log2FC| >= %.2f ...",
         gatekeep_fdr, gatekeep_log2FC)
    
    gkeep <- rownames(gtab_raw)[gtab_raw$FDR <= gatekeep_fdr]
    
    if (gatekeep_log2FC > 0 && length(posthoc_raw)) {
      all_genes <- Reduce(union, lapply(posthoc_raw, rownames))
      max_abs <- setNames(rep(-Inf, length(all_genes)), all_genes)
      for (df in posthoc_raw) {
        if (!"logFC" %in% colnames(df)) next
        common <- intersect(names(max_abs), rownames(df))
        max_abs[common] <- pmax(max_abs[common], abs(df[common,"logFC"]), na.rm = TRUE)
      }
      big_eff <- names(max_abs)[max_abs >= gatekeep_log2FC]
      gkeep <- intersect(gkeep, big_eff)
    }
    gate_mask <- rownames(dge) %in% gkeep
    .msg("[gatekeep] %d gene(s) pass the global filter.", sum(gate_mask))
    
    posthoc_gate <- lapply(posthoc_raw, function(df){
      out <- df
      ok <- (rownames(df) %in% names(gate_mask)[gate_mask]) & !is.na(df$PValue)
      newFDR <- rep(NA_real_, nrow(df))
      if (any(ok)) newFDR[ok] <- p.adjust(df$PValue[ok], method = p_adjust_method)
      out$FDR_gate <- newFDR
      out
    })
    names(posthoc_gate) <- names(posthoc_raw)
  }
  
  # ---- Wrap with IDs ---------------------------------------------------------
  posthoc     <- lapply(posthoc_raw, .add_ids)
  if (!is.null(posthoc_gate)) posthoc_gate <- lapply(posthoc_gate, .add_ids)
  
  # ---- D) Combined table -----------------------------------------------------
  .msg("[combine] Building combined table...")
  id_cols <- c("name_id","ensembl_id","gene_symbol")
  ANOVA_global <- .pick(.add_ids(gtab_raw), c(id_cols,"logCPM","F","PValue","FDR"))
  colnames(ANOVA_global) <- c(id_cols,"ANOVA_global_logCPM","ANOVA_global_F","ANOVA_global_PValue","ANOVA_global_FDR")
  ANOVA_first <- .pick(tab_f1, c(id_cols,"F","PValue","FDR"))
  colnames(ANOVA_first) <- c(id_cols,"ANOVA_first_F","ANOVA_first_PValue","ANOVA_first_FDR")
  ANOVA_second <- .pick(tab_f2, c(id_cols,"F","PValue","FDR"))
  colnames(ANOVA_second) <- c(id_cols,"ANOVA_second_F","ANOVA_second_PValue","ANOVA_second_FDR")
  ANOVA_inter <- .pick(tab_int, c(id_cols,"F","PValue","FDR"))
  colnames(ANOVA_inter) <- c(id_cols,"ANOVA_interaction_F","ANOVA_interaction_PValue","ANOVA_interaction_FDR")
  
  ph_blocks <- lapply(names(posthoc), function(nm) .rename_block(posthoc[[nm]], nm, id_cols))
  combined <- Reduce(function(a,b) merge(a,b, by = id_cols, all = TRUE),
                     c(list(ANOVA_global, ANOVA_first, ANOVA_second, ANOVA_inter), ph_blocks))
  
  .msg("[combine] Combined table: %d genes x %d columns.", nrow(combined), ncol(combined))
  .msg("[done] Two-factor ANOVA + post-hoc finished.")
  
  # ---- Out -------------------------------------------------------------------
  list(
    anova        = list(global=tab_global, main_first=tab_f1, main_second=tab_f2, interaction=tab_int),
    posthoc      = posthoc,       # original edgeR results
    posthoc_gate = posthoc_gate,  # if gatekeep=TRUE: results with new FDR_gate
    combined     = combined,
    meta = list(
      first_factor  = list(name=f1, levels=levels(smeta[[f1]])),
      second_factor = list(name=f2, levels=levels(smeta[[f2]])),
      robust = robust,
      gatekeep = gatekeep,
      gatekeep_fdr = gatekeep_fdr,
      gatekeep_log2FC = gatekeep_log2FC,
      p_adjust_method = p_adjust_method
    )
  )
}

edger_twofactor_anova_and_posthoc <- function(
  counts,                    
  sample_meta,               
  sample_col         = "sample_id",
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  robust             = TRUE,
  verbose            = TRUE,
  gatekeep           = FALSE,          
  gatekeep_fdr       = 0.05,
  gatekeep_log2FC    = 0,
  p_adjust_method    = "BH"
){
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---- Helpers ---------------------------------------------------------------
  .parse_ids <- function(rn_vec){
    # split tokens by space or hyphen
    split_tokens <- strsplit(rn_vec, "[ -]+", perl = TRUE)
    peak_id <- vapply(split_tokens, function(tok){
      if (length(tok) >= 2) paste(tok[1:2], collapse = "-") else tok[1]
    }, character(1))
    gene_symbol <- vapply(split_tokens, function(tok){
      if (length(tok) >= 1) tok[length(tok)] else NA_character_
    }, character(1))
    data.frame(
      name_id     = rn_vec,
      peak_id     = peak_id,
      gene_symbol = gene_symbol,
      stringsAsFactors = FALSE
    )
  }
  .add_ids <- function(df){
    stopifnot(!is.null(rownames(df)))
    cbind(.parse_ids(rownames(df)), df, row.names = NULL)
  }
  .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
  
  .rename_block <- function(df, block, id_cols, drop_logFC = FALSE){
    x <- .pick(df, c(id_cols, "logFC","logCPM","F","PValue","FDR"))
    if (drop_logFC && "logFC" %in% colnames(x)) x$logFC <- NULL
    colnames(x) <- c(
      id_cols,
      if (!drop_logFC) paste0(block,"_logFC"),
      paste0(block,"_logCPM"),
      paste0(block,"_F"),
      paste0(block,"_PValue"),
      paste0(block,"_FDR")
    )
    x
  }
  
  # ---- Metadata --------------------------------------------------------------
  .msg("[meta] Aligning metadata to count matrix...")
  smeta <- as.data.frame(sample_meta)
  idx <- match(colnames(counts), smeta[[sample_col]])
  if (anyNA(idx)) stop("Missing metadata for samples.")
  smeta <- smeta[idx, , drop = FALSE]
  
  if (is.null(first_levels))  first_levels  <- unique(as.character(smeta[[first_factor_col]]))
  if (is.null(second_levels)) second_levels <- unique(as.character(smeta[[second_factor_col]]))
  smeta[[first_factor_col]]  <- factor(smeta[[first_factor_col]],  levels = first_levels)
  smeta[[second_factor_col]] <- factor(smeta[[second_factor_col]], levels = second_levels)
  f1 <- first_factor_col; f2 <- second_factor_col
  .msg("[meta] Samples: %d | %s levels: %s | %s levels: %s",
       ncol(counts), f1, paste(levels(smeta[[f1]]), collapse=", "),
       f2, paste(levels(smeta[[f2]]), collapse=", "))
  
  # ---- edgeR preprocessing ---------------------------------------------------
  .msg("[edgeR] DGEList, filterByExpr, TMM normalization...")
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = interaction(dge$samples[[f1]], dge$samples[[f2]], drop=TRUE))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  .msg("[edgeR] Kept %d/%d features after filtering.", nrow(dge), nrow(counts))
  
  # ---- A) ANOVA-like model ---------------------------------------------------
  .msg("[model] Fitting ANOVA-like model: ~ %s * %s ...", f1, f2)
  designA <- stats::model.matrix(stats::as.formula(paste0("~ ", f1, " * ", f2)), data = dge$samples)
  dgeA    <- edgeR::estimateDisp(dge, designA)
  fitA    <- edgeR::glmQLFit(dgeA, designA, robust = robust)
  
  .msg("[test] Omnibus QL F-test...")
  q_global <- edgeR::glmQLFTest(fitA, coef = 2:ncol(designA))
  tab_global <- .add_ids(edgeR::topTags(q_global, n=Inf)$table)
  
  .msg("[test] Main effect %s...", f1)
  q_f1 <- edgeR::glmQLFTest(fitA, coef = grep(paste0("^", f1), colnames(designA)))
  tab_f1 <- .add_ids(edgeR::topTags(q_f1, n=Inf)$table)
  
  .msg("[test] Main effect %s...", f2)
  q_f2 <- edgeR::glmQLFTest(fitA, coef = grep(paste0("^", f2), colnames(designA)))
  tab_f2 <- .add_ids(edgeR::topTags(q_f2, n=Inf)$table)
  
  .msg("[test] Interaction %s:%s...", f1, f2)
  q_int <- edgeR::glmQLFTest(fitA, coef = grep(":", colnames(designA)))
  tab_int <- .add_ids(edgeR::topTags(q_int, n=Inf)$table)
  if ("logFC" %in% colnames(tab_int)) tab_int$logFC <- NULL
  
  # ---- B) Simple effects model (posthoc) ------------------------------------
  .msg("[model] Fitting cell-means model for simple effects...")
  designM <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", f1, ":", f2)), data = dge$samples)
  dgeM    <- edgeR::estimateDisp(dge, designM)
  fitM    <- edgeR::glmQLFit(dgeM, designM, robust = robust)
  colsM   <- colnames(designM)
  
  cell_name  <- function(A, B) paste0(f1, A, ":", f2, B)
  which_cell <- function(A, B) match(cell_name(A, B), colsM)
  
  # build simple-effect contrasts
  contrasts_idx <- list()
  levA <- levels(smeta[[f1]]); levB <- levels(smeta[[f2]])
  for (a in levA){
    for (i in seq_len(length(levB)-1)){
      for (j in (i+1):length(levB)){
        contrasts_idx[[sprintf("SECOND_%s-%s|FIRST_%s", levB[j], levB[i], a)]] <- 
          c(which_cell(a, levB[j]), which_cell(a, levB[i]))
      }
    }
  }
  for (b in levB){
    for (i in seq_len(length(levA)-1)){
      for (j in (i+1):length(levA)){
        contrasts_idx[[sprintf("FIRST_%s-%s|SECOND_%s", levA[j], levA[i], b)]] <- 
          c(which_cell(levA[j], b), which_cell(levA[i], b))
      }
    }
  }
  
  run_contrast <- function(ix){
    contr <- numeric(ncol(designM))
    contr[ix[1]] <-  1
    contr[ix[2]] <- -1
    edgeR::glmQLFTest(fitM, contrast = contr)
  }
  
  .msg("[posthoc] Computing %d simple-effect contrasts...", length(contrasts_idx))
  posthoc_raw <- lapply(contrasts_idx, function(ix) edgeR::topTags(run_contrast(ix), n=Inf)$table)
  names(posthoc_raw) <- names(contrasts_idx)
  
  # ---- C) Gatekeeping --------------------------------------------------------
  posthoc_gate <- NULL
  if (isTRUE(gatekeep)) {
    .msg("[gatekeep] Applying thresholds: FDR <= %.3f ; |log2FC| >= %.2f ...",
         gatekeep_fdr, gatekeep_log2FC)
    gkeep <- rownames(tab_global)[tab_global$FDR <= gatekeep_fdr]
    if (gatekeep_log2FC > 0 && length(posthoc_raw)) {
      all_genes <- Reduce(union, lapply(posthoc_raw, rownames))
      max_abs <- setNames(rep(-Inf, length(all_genes)), all_genes)
      for (df in posthoc_raw) {
        if (!"logFC" %in% colnames(df)) next
        common <- intersect(names(max_abs), rownames(df))
        max_abs[common] <- pmax(max_abs[common], abs(df[common,"logFC"]), na.rm=TRUE)
      }
      big_eff <- names(max_abs)[max_abs >= gatekeep_log2FC]
      gkeep <- intersect(gkeep, big_eff)
    }
    posthoc_gate <- lapply(posthoc_raw, function(df){
      out <- df
      ok <- rownames(df) %in% gkeep
      newFDR <- rep(NA_real_, nrow(df))
      if (any(ok)) newFDR[ok] <- p.adjust(df$PValue[ok], method=p_adjust_method)
      out$FDR_gate <- newFDR
      out
    })
    names(posthoc_gate) <- names(posthoc_raw)
  }
  
  # ---- Wrap with IDs ---------------------------------------------------------
  posthoc     <- lapply(posthoc_raw, .add_ids)
  if (!is.null(posthoc_gate)) posthoc_gate <- lapply(posthoc_gate, .add_ids)
  
  # ---- D) Combined table -----------------------------------------------------
  .msg("[combine] Building combined table...")
  id_cols <- c("name_id","peak_id","gene_symbol")
  ANOVA_global <- .rename_block(tab_global,"ANOVA_global",id_cols)
  ANOVA_first  <- .rename_block(tab_f1,"ANOVA_first",id_cols)
  ANOVA_second <- .rename_block(tab_f2,"ANOVA_second",id_cols)
  ANOVA_inter  <- .rename_block(tab_int,"ANOVA_interaction",id_cols, drop_logFC=TRUE)
  ph_blocks    <- lapply(names(posthoc), function(nm) .rename_block(posthoc[[nm]], nm, id_cols))
  
  combined <- Reduce(function(a,b) merge(a,b, by=id_cols, all=TRUE),
                     c(list(ANOVA_global, ANOVA_first, ANOVA_second, ANOVA_inter), ph_blocks))
  
  .msg("[combine] Combined table: %d x %d.", nrow(combined), ncol(combined))
  .msg("[done] Two-factor ANOVA + post-hoc finished.")
  
  # ---- Output ---------------------------------------------------------------
  list(
    anova        = list(global=tab_global, main_first=tab_f1, main_second=tab_f2, interaction=tab_int),
    posthoc      = posthoc,
    posthoc_gate = posthoc_gate,
    combined     = combined,
    meta = list(
      first_factor  = list(name=f1, levels=levels(smeta[[f1]])),
      second_factor = list(name=f2, levels=levels(smeta[[f2]])),
      robust = robust,
      gatekeep = gatekeep,
      gatekeep_fdr = gatekeep_fdr,
      gatekeep_log2FC = gatekeep_log2FC,
      p_adjust_method = p_adjust_method
    )
  )
}


edger_twofactor_anova_and_posthoc <- function(
  counts,                    
  sample_meta,               
  sample_col         = "sample_id",
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  robust             = TRUE,
  verbose            = TRUE,
  gatekeep           = FALSE,          
  gatekeep_fdr       = 0.05,
  gatekeep_log2FC    = 0,
  p_adjust_method    = "BH"
){
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---- Helpers ---------------------------------------------------------------
  .parse_ids <- function(rn_vec){
    split_tokens <- strsplit(rn_vec, "[ -]+", perl = TRUE)
    peak_id <- vapply(split_tokens, function(tok){
      if (length(tok) >= 2) paste(tok[1:2], collapse = "-") else tok[1]
    }, character(1))
    gene_symbol <- vapply(split_tokens, function(tok){
      if (length(tok) >= 1) tok[length(tok)] else NA_character_
    }, character(1))
    data.frame(
      name_id     = rn_vec,
      peak_id     = peak_id,
      gene_symbol = gene_symbol,
      stringsAsFactors = FALSE
    )
  }
  .add_ids <- function(df){
    stopifnot(!is.null(rownames(df)))
    cbind(.parse_ids(rownames(df)), df, row.names = NULL)
  }
  .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
  
  .rename_block <- function(df, block, id_cols, drop_logFC = FALSE){
    measures <- c("logFC","logCPM","F","PValue","FDR")
    present  <- intersect(measures, colnames(df))
    if (drop_logFC) present <- setdiff(present, "logFC")
    x <- .pick(df, c(id_cols, present))
    new_meas_names <- setNames(paste0(block, "_", present), present)
    colnames(x) <- c(id_cols, unname(new_meas_names[present]))
    x
  }
  
  # ---- Metadata --------------------------------------------------------------
  smeta <- as.data.frame(sample_meta)
  idx <- match(colnames(counts), smeta[[sample_col]])
  if (anyNA(idx)) stop("Missing metadata for samples.")
  smeta <- smeta[idx, , drop = FALSE]
  
  # ---- Drop all-NA columns ---------------------------------------------------
  na_cols <- which(apply(counts, 2, function(x) all(is.na(x))))
  if (length(na_cols) > 0) {
    .msg("[preprocess] Dropping %d samples with all-NA counts: %s",
         length(na_cols), paste(colnames(counts)[na_cols], collapse=", "))
    counts <- counts[, -na_cols, drop=FALSE]
    smeta  <- smeta[-na_cols, , drop=FALSE]
  }
  if (ncol(counts) == 0) {
    warning("All columns dropped (all-NA). Returning empty result.")
    return(list(
      anova        = list(global=NULL, main_first=NULL, main_second=NULL, interaction=NULL),
      posthoc      = list(),
      posthoc_gate = list(),
      combined     = data.frame()
    ))
  }
  
  # ---- Factors ---------------------------------------------------------------
  if (is.null(first_levels))  first_levels  <- unique(as.character(smeta[[first_factor_col]]))
  if (is.null(second_levels)) second_levels <- unique(as.character(smeta[[second_factor_col]]))
  smeta[[first_factor_col]]  <- factor(smeta[[first_factor_col]],  levels = first_levels)
  smeta[[second_factor_col]] <- factor(smeta[[second_factor_col]], levels = second_levels)
  f1 <- first_factor_col; f2 <- second_factor_col
  .msg("[meta] Samples: %d | %s levels: %s | %s levels: %s",
       ncol(counts), f1, paste(levels(smeta[[f1]]), collapse=", "),
       f2, paste(levels(smeta[[f2]]), collapse=", "))
  
  # ---- edgeR preprocessing ---------------------------------------------------
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = interaction(dge$samples[[f1]], dge$samples[[f2]], drop=TRUE))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  
  .msg("[edgeR] Kept %d/%d features after filtering.", nrow(dge), nrow(counts))
  
  # ---- A) ANOVA-like model ---------------------------------------------------
  designA <- stats::model.matrix(stats::as.formula(paste0("~ ", f1, " * ", f2)), data = dge$samples)
  dgeA    <- edgeR::estimateDisp(dge, designA)
  fitA    <- edgeR::glmQLFit(dgeA, designA, robust = robust)
  
  # Global omnibus test
  q_global <- edgeR::glmQLFTest(fitA, coef = 2:ncol(designA))
  tab_global <- .add_ids(edgeR::topTags(q_global, n=Inf)$table)
  
  # Main effect: first factor
  f1_level2 <- levels(smeta[[f1]])[2]
  f1_coef   <- paste0(f1, f1_level2)
  f1_idx    <- which(colnames(designA) == f1_coef)
  q_f1      <- edgeR::glmQLFTest(fitA, coef = f1_idx)
  tab_f1    <- .add_ids(edgeR::topTags(q_f1, n=Inf)$table)
  
  # Main effect: second factor
  f2_level2 <- levels(smeta[[f2]])[2]
  f2_coef   <- paste0(f2, f2_level2)
  f2_idx    <- which(colnames(designA) == f2_coef)
  q_f2      <- edgeR::glmQLFTest(fitA, coef = f2_idx)
  tab_f2    <- .add_ids(edgeR::topTags(q_f2, n=Inf)$table)
  
  # Interaction
  inter_idx <- grep(":", colnames(designA))
  q_int     <- edgeR::glmQLFTest(fitA, coef = inter_idx)
  tab_int   <- .add_ids(edgeR::topTags(q_int, n=Inf)$table)
  if ("logFC" %in% colnames(tab_int)) tab_int$logFC <- NULL
  
  # ---- B) Simple effects (posthoc) ------------------------------------------
  designM <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", f1, ":", f2)), data = dge$samples)
  dgeM    <- edgeR::estimateDisp(dge, designM)
  fitM    <- edgeR::glmQLFit(dgeM, designM, robust = robust)
  colsM   <- colnames(designM)
  
  cell_name  <- function(A, B) paste0(f1, A, ":", f2, B)
  which_cell <- function(A, B) match(cell_name(A, B), colsM)
  
  contrasts_idx <- list()
  levA <- levels(smeta[[f1]]); levB <- levels(smeta[[f2]])
  for (a in levA){
    for (i in seq_len(length(levB)-1)){
      for (j in (i+1):length(levB)){
        contrasts_idx[[sprintf("SECOND_%s-%s|FIRST_%s", levB[j], levB[i], a)]] <- 
          c(which_cell(a, levB[j]), which_cell(a, levB[i]))
      }
    }
  }
  for (b in levB){
    for (i in seq_len(length(levA)-1)){
      for (j in (i+1):length(levA)){
        contrasts_idx[[sprintf("FIRST_%s-%s|SECOND_%s", levA[j], levA[i], b)]] <- 
          c(which_cell(levA[j], b), which_cell(levA[i], b))
      }
    }
  }
  
  run_contrast <- function(ix){
    contr <- numeric(ncol(designM))
    contr[ix[1]] <-  1
    contr[ix[2]] <- -1
    edgeR::glmQLFTest(fitM, contrast = contr)
  }
  
  posthoc_raw <- lapply(contrasts_idx, function(ix) edgeR::topTags(run_contrast(ix), n=Inf)$table)
  names(posthoc_raw) <- names(contrasts_idx)
  
  # ---- Gatekeeping -----------------------------------------------------------
  posthoc_gate <- NULL
  if (isTRUE(gatekeep)) {
    gkeep <- rownames(tab_global)[tab_global$FDR <= gatekeep_fdr]
    if (gatekeep_log2FC > 0 && length(posthoc_raw)) {
      all_genes <- Reduce(union, lapply(posthoc_raw, rownames))
      max_abs <- setNames(rep(-Inf, length(all_genes)), all_genes)
      for (df in posthoc_raw) {
        if (!"logFC" %in% colnames(df)) next
        common <- intersect(names(max_abs), rownames(df))
        max_abs[common] <- pmax(max_abs[common], abs(df[common,"logFC"]), na.rm=TRUE)
      }
      big_eff <- names(max_abs)[max_abs >= gatekeep_log2FC]
      gkeep <- intersect(gkeep, big_eff)
    }
    posthoc_gate <- lapply(posthoc_raw, function(df){
      out <- df
      ok <- rownames(df) %in% gkeep & !is.na(df$PValue)
      newFDR <- rep(NA_real_, nrow(df))
      if (any(ok)) newFDR[ok] <- p.adjust(df$PValue[ok], method=p_adjust_method)
      out$FDR_gate <- newFDR
      out
    })
    names(posthoc_gate) <- names(posthoc_raw)
  }
  
  # ---- Wrap IDs --------------------------------------------------------------
  posthoc     <- lapply(posthoc_raw, .add_ids)
  if (!is.null(posthoc_gate)) posthoc_gate <- lapply(posthoc_gate, .add_ids)
  
  # ---- Combined --------------------------------------------------------------
  id_cols <- c("name_id","peak_id","gene_symbol")
  ANOVA_global <- .rename_block(tab_global,"ANOVA_global",id_cols)
  ANOVA_first  <- .rename_block(tab_f1,"ANOVA_first",id_cols)
  ANOVA_second <- .rename_block(tab_f2,"ANOVA_second",id_cols)
  ANOVA_inter  <- .rename_block(tab_int,"ANOVA_interaction",id_cols, drop_logFC=TRUE)
  ph_blocks    <- lapply(names(posthoc), function(nm) .rename_block(posthoc[[nm]], nm, id_cols))
  
  combined <- Reduce(function(a,b) merge(a,b, by=id_cols, all=TRUE),
                     c(list(ANOVA_global, ANOVA_first, ANOVA_second, ANOVA_inter), ph_blocks))
  
  .msg("[done] Two-factor ANOVA + post-hoc finished.")
  
  list(
    anova        = list(global=tab_global, main_first=tab_f1, main_second=tab_f2, interaction=tab_int),
    posthoc      = posthoc,
    posthoc_gate = posthoc_gate,
    combined     = combined
  )
}

