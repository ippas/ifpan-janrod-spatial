edger_twofactor_anova_and_posthoc_v2 <- function(
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
  p_adjust_method    = "BH",
  main_effects       = c("baseline","marginal"),
  main_effects_weights = c("equal","proportional")
){
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  main_effects        <- match.arg(main_effects)
  main_effects_weights<- match.arg(main_effects_weights)
  
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
  levA <- levels(smeta[[f1]])
  levB <- levels(smeta[[f2]])
  .msg("[meta] Samples: %d | %s levels: %s | %s levels: %s",
       ncol(counts), f1, paste(levA, collapse=", "),
       f2, paste(levB, collapse=", "))
  
  # ---- edgeR preprocessing ---------------------------------------------------
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = interaction(dge$samples[[f1]], dge$samples[[f2]], drop=TRUE))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  .msg("[edgeR] Kept %d/%d features after filtering.", nrow(dge), nrow(counts))
  
  # ---- A) ANOVA-like model (global + interaction; and optionally baseline mains) ----
  designA <- stats::model.matrix(stats::as.formula(paste0("~ ", f1, " * ", f2)), data = dge$samples)
  dgeA    <- edgeR::estimateDisp(dge, designA)
  fitA    <- edgeR::glmQLFit(dgeA, designA, robust = robust)
  
  # Global omnibus test (wszystkie współczynniki poza interceptem)
  q_global   <- edgeR::glmQLFTest(fitA, coef = 2:ncol(designA))
  tab_global <- .add_ids(edgeR::topTags(q_global, n=Inf)$table)
  
  # Interaction (F-test wielu współczynników)
  inter_idx <- grep(":", colnames(designA))
  q_int     <- edgeR::glmQLFTest(fitA, coef = inter_idx)
  tab_int   <- .add_ids(edgeR::topTags(q_int, n=Inf)$table)
  if ("logFC" %in% colnames(tab_int)) tab_int$logFC <- NULL
  
  # (Opcjonalnie) baseline main effects z modelu z interakcją
  get_baseline_main <- function(which_factor = c("first","second")){
    which_factor <- match.arg(which_factor)
    if (which_factor == "first") {
      if (length(levA) < 2) stop("Need >=2 levels in first factor.")
      coef_name <- paste0(f1, levA[2])
    } else {
      if (length(levB) < 2) stop("Need >=2 levels in second factor.")
      coef_name <- paste0(f2, levB[2])
    }
    idx <- match(coef_name, colnames(designA))
    if (is.na(idx)) stop("Couldn't find coefficient: ", coef_name)
    .add_ids(edgeR::topTags(edgeR::glmQLFTest(fitA, coef = idx), n=Inf)$table)
  }
  
  # ---- B) Macierz komórkowa dla prostych efektów i kontrastów marginalnych ----
  designM <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", f1, ":", f2)), data = dge$samples)
  dgeM    <- edgeR::estimateDisp(dge, designM)
  fitM    <- edgeR::glmQLFit(dgeM, designM, robust = robust)
  colsM   <- colnames(designM)
  
  cell_name  <- function(A, B) paste0(f1, A, ":", f2, B)
  which_cell <- function(A, B) match(cell_name(A, B), colsM)
  
  # ---- Posthoc: wszystkie proste efekty (jak w Twojej wersji) ----------------
  contrasts_idx <- list()
  # SECOND: porównania poziomów f2 w ramach stałego f1
  for (a in levA){
    for (i in seq_len(length(levB)-1)){
      for (j in (i+1):length(levB)){
        contrasts_idx[[sprintf("SECOND_%s-%s|FIRST_%s", levB[j], levB[i], a)]] <- 
          c(which_cell(a, levB[j]), which_cell(a, levB[i]))
      }
    }
  }
  # FIRST: porównania poziomów f1 w ramach stałego f2
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
  
  # ---- Gatekeeping (opcjonalnie) --------------------------------------------
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
  
  # ---- Marginalne główne efekty (uśrednianie po drugim czynniku) -------------
  # Kontrasty dla 2-poziomowych czynników (drugi vs pierwszy), uśrednione po drugim czynniku.
  build_marginal_contrast_f1 <- function(){
    if (length(levA) != 2) {
      warning("main_effects='marginal' obecnie implementuje kontrast dla 2 poziomów pierwszego czynnika (levA[2] - levA[1]) uśredniony po wszystkich poziomach drugiego.")
    }
    # wagi po poziomach f2
    if (main_effects_weights == "equal") {
      wB <- setNames(rep(1/length(levB), length(levB)), levB)
    } else {
      # proporcjonalnie do liczby próbek w każdym poziomie f2
      tb <- table(smeta[[f2]])
      wB <- as.numeric(tb) / sum(tb)
      names(wB) <- names(tb)
    }
    contr <- numeric(ncol(designM))
    for (b in levB){
      contr[which_cell(levA[2], b)] <- contr[which_cell(levA[2], b)] + wB[b]
      contr[which_cell(levA[1], b)] <- contr[which_cell(levA[1], b)] - wB[b]
    }
    contr
  }
  build_marginal_contrast_f2 <- function(){
    if (length(levB) != 2) {
      warning("main_effects='marginal' obecnie implementuje kontrast dla 2 poziomów drugiego czynnika (levB[2] - levB[1]) uśredniony po wszystkich poziomach pierwszego.")
    }
    if (main_effects_weights == "equal") {
      wA <- setNames(rep(1/length(levA), length(levA)), levA)
    } else {
      tb <- table(smeta[[f1]])
      wA <- as.numeric(tb) / sum(tb)
      names(wA) <- names(tb)
    }
    contr <- numeric(ncol(designM))
    for (a in levA){
      contr[which_cell(a, levB[2])] <- contr[which_cell(a, levB[2])] + wA[a]
      contr[which_cell(a, levB[1])] <- contr[which_cell(a, levB[1])] - wA[a]
    }
    contr
  }
  
  if (main_effects == "baseline") {
    tab_f1 <- get_baseline_main("first")
    tab_f2 <- get_baseline_main("second")
  } else {
    # marginalne: używamy fitM i kontrastów z uśrednianiem
    q_f1m    <- edgeR::glmQLFTest(fitM, contrast = build_marginal_contrast_f1())
    tab_f1   <- .add_ids(edgeR::topTags(q_f1m, n=Inf)$table)
    q_f2m    <- edgeR::glmQLFTest(fitM, contrast = build_marginal_contrast_f2())
    tab_f2   <- .add_ids(edgeR::topTags(q_f2m, n=Inf)$table)
  }
  
  # ---- Wrap IDs --------------------------------------------------------------
  posthoc     <- lapply(posthoc_raw, .add_ids)
  if (!is.null(posthoc_gate)) posthoc_gate <- lapply(posthoc_gate, .add_ids)
  
  # ---- Combined --------------------------------------------------------------
  id_cols <- c("name_id","peak_id","gene_symbol")
  ANOVA_global <- .rename_block(tab_global,"ANOVA_global",id_cols)
  ANOVA_first  <- .rename_block(tab_f1, ifelse(main_effects=="baseline","ANOVA_first","ANOVA_first_marginal"), id_cols)
  ANOVA_second <- .rename_block(tab_f2, ifelse(main_effects=="baseline","ANOVA_second","ANOVA_second_marginal"), id_cols)
  ANOVA_inter  <- .rename_block(tab_int,"ANOVA_interaction",id_cols, drop_logFC=TRUE)
  ph_blocks    <- lapply(names(posthoc), function(nm) .rename_block(posthoc[[nm]], nm, id_cols))
  
  combined <- Reduce(function(a,b) merge(a,b, by=id_cols, all=TRUE),
                     c(list(ANOVA_global, ANOVA_first, ANOVA_second, ANOVA_inter), ph_blocks))
  
  .msg("[done] Two-factor ANOVA + post-hoc finished. Main effects: %s (%s weights).",
       main_effects, ifelse(main_effects_weights=="equal","equal","proportional"))
  
  list(
    anova        = list(
      global       = tab_global,
      main_first   = tab_f1,     # baseline lub marginal – zgodnie z argumentem
      main_second  = tab_f2,     # baseline lub marginal
      interaction  = tab_int
    ),
    posthoc      = posthoc,
    posthoc_gate = posthoc_gate,
    combined     = combined
  )
}
