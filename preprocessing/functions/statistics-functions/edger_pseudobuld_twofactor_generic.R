edger_twofactor_anova <- function(
  counts,
  sample_meta,
  sample_col         = "sample_id",
  first_factor_col   = "genotype",
  second_factor_col  = "treatment",
  first_levels       = NULL,
  second_levels      = NULL,
  robust             = TRUE,
  verbose            = TRUE
){
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # Basic checks
  stopifnot(is.matrix(counts) || inherits(counts, "Matrix"))
  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) stop("'counts' must have rownames = gene IDs.")
  if (is.null(colnames(counts))) stop("'counts' must have colnames = sample IDs.")
  stopifnot(all(c(sample_col, first_factor_col, second_factor_col) %in% colnames(sample_meta)))
  
  # Helpers
  .add_ids <- function(df){
    rn <- rownames(df)
    ensembl_id  <- sub("^([^_-]+)[_-].*$", "\\1", rn)
    gene_symbol <- ifelse(grepl("^[^_-]+[_-]", rn),
                          sub("^[^_-]+[_-](.*)$", "\\1", rn),
                          NA_character_)
    cbind(data.frame(name_id=rn, ensembl_id=ensembl_id, gene_symbol=gene_symbol, stringsAsFactors=FALSE),
          df, row.names=NULL)
  }
  .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
  .with_prefix <- function(df, prefix, id_cols){
    x <- .pick(df, c(id_cols, "logFC","logCPM","F","PValue","FDR"))
    nm <- colnames(x); nm <- ifelse(nm %in% id_cols, nm, paste0(prefix, nm))
    colnames(x) <- nm; x
  }
  
  # Metadata alignment
  .msg("[meta] Aligning metadata to count matrix...")
  smeta <- as.data.frame(sample_meta)
  if (any(duplicated(smeta[[sample_col]]))) {
    .msg("[meta] Duplicates in '%s' — keeping first occurrence.", sample_col)
    smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
  }
  idx <- match(colnames(counts), smeta[[sample_col]])
  if (anyNA(idx)) stop("Missing metadata for samples: ", paste(colnames(counts)[is.na(idx)], collapse=", "))
  smeta <- smeta[idx, , drop = FALSE]
  
  if (is.null(first_levels))  first_levels  <- unique(as.character(smeta[[first_factor_col]]))
  if (is.null(second_levels)) second_levels <- unique(as.character(smeta[[second_factor_col]]))
  smeta[[first_factor_col]]  <- factor(smeta[[first_factor_col]],  levels = first_levels)
  smeta[[second_factor_col]] <- factor(smeta[[second_factor_col]], levels = second_levels)
  .msg("[meta] Samples: %d | First '%s': %s | Second '%s': %s",
       ncol(counts), first_factor_col, paste(levels(smeta[[first_factor_col]]), collapse=", "),
       second_factor_col, paste(levels(smeta[[second_factor_col]]), collapse=", "))
  
  # edgeR core
  .msg("[edgeR] DGEList, filterByExpr, TMM normalization...")
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = interaction(dge$samples[[first_factor_col]],
                                                       dge$samples[[second_factor_col]], drop=TRUE))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  .msg("[edgeR] Kept %d/%d genes after filtering.", nrow(dge), nrow(counts))
  
  # A) ANOVA: ~ first * second
  f1 <- first_factor_col; f2 <- second_factor_col
  .msg("[model] Fitting ANOVA-like model: ~ %s * %s ...", f1, f2)
  design_anova <- stats::model.matrix(stats::as.formula(paste0("~ ", f1, " * ", f2)), data = dge$samples)
  anova_cols <- colnames(design_anova)
  dgeA <- edgeR::estimateDisp(dge, design_anova)
  fitA <- edgeR::glmQLFit(dgeA, design_anova, robust = robust)
  
  .msg("[test] Omnibus QL F-test (all non-intercept coefficients)...")
  q_global <- edgeR::glmQLFTest(fitA, coef = 2:ncol(design_anova))
  tab_global <- .add_ids(edgeR::topTags(q_global, n=Inf)$table)
  
  .msg("[test] Main effect '%s'...", f1)
  f1_idx <- grep(paste0("^", f1), anova_cols)
  q_f1   <- edgeR::glmQLFTest(fitA, coef = f1_idx)
  tab_f1 <- .add_ids(edgeR::topTags(q_f1, n=Inf)$table)
  
  .msg("[test] Main effect '%s'...", f2)
  f2_idx <- grep(paste0("^", f2), anova_cols)
  q_f2   <- edgeR::glmQLFTest(fitA, coef = f2_idx)
  tab_f2 <- .add_ids(edgeR::topTags(q_f2, n=Inf)$table)
  
  .msg("[test] Interaction '%s:%s'...", f1, f2)
  inter_idx <- grep(":", anova_cols)
  if (!length(inter_idx)) stop("No interaction columns found in design matrix.")
  q_int   <- edgeR::glmQLFTest(fitA, coef = inter_idx)
  tab_int <- .add_ids(edgeR::topTags(q_int, n=Inf)$table)
  
  # B) Post-hoc simple effects: ~ 0 + first:second
  .msg("[model] Fitting cell-means model for simple effects: ~ 0 + %s:%s ...", f1, f2)
  design_means <- stats::model.matrix(~ 0 + mouse_genotype:treatment, data = dge$samples)
  means_cols <- colnames(design_means)
  dgeM <- edgeR::estimateDisp(dge, design_means)      # <- NEW: estimate dispersions for this design
  fitM <- edgeR::glmQLFit(dgeM, design_means, robust = robust)
  cell_name  <- function(a, b) paste0(f1, a, ":", f2, b)
  which_cell <- function(a, b) match(cell_name(a, b), means_cols)
  
  posthoc <- list()
  if (length(second_levels) >= 2){
    .msg("[posthoc] SECOND within each FIRST level...")
    for (a in first_levels){
      for (i in seq_len(length(second_levels)-1)){
        for (j in (i+1):length(second_levels)){
          b1 <- second_levels[i]; b2 <- second_levels[j]
          i1 <- which_cell(a,b1); i2 <- which_cell(a,b2)
          if (is.na(i1) || is.na(i2)) next
          contr <- numeric(ncol(design_means)); contr[i2] <- 1; contr[i1] <- -1
          q <- edgeR::glmQLFTest(fitM, contrast = contr)
          posthoc[[paste0("SECOND_", b2, "-", b1, "|FIRST_", a)]] <- .add_ids(edgeR::topTags(q, n=Inf)$table)
        }
      }
    }
  }
  if (length(first_levels) >= 2){
    .msg("[posthoc] FIRST within each SECOND level...")
    for (b in second_levels){
      for (i in seq_len(length(first_levels)-1)){
        for (j in (i+1):length(first_levels)){
          a1 <- first_levels[i]; a2 <- first_levels[j]
          i1 <- which_cell(a1,b); i2 <- which_cell(a2,b)
          if (is.na(i1) || is.na(i2)) next
          contr <- numeric(ncol(design_means)); contr[i2] <- 1; contr[i1] <- -1
          q <- edgeR::glmQLFTest(fitM, contrast = contr)
          posthoc[[paste0("FIRST_", a2, "-", a1, "|SECOND_", b)]] <- .add_ids(edgeR::topTags(q, n=Inf)$table)
        }
      }
    }
  }
  
  # C) Combined table
  .msg("[combine] Building combined table...")
  id_cols <- c("name_id","ensembl_id","gene_symbol")
  g <- .pick(tab_global, c(id_cols, "logCPM","F","PValue","FDR"))
  colnames(g) <- sub("^logCPM$","ANOVA_global_logCPM", colnames(g))
  colnames(g) <- sub("^F$","ANOVA_global_F", colnames(g))
  colnames(g) <- sub("^PValue$","ANOVA_global_PValue", colnames(g))
  colnames(g) <- sub("^FDR$","ANOVA_global_FDR", colnames(g))
  
  main_f1 <- .pick(tab_f1, c(id_cols, "F","PValue","FDR"))
  colnames(main_f1)[-(1:3)] <- paste0("ANOVA_first_", colnames(main_f1)[-(1:3)])
  
  main_f2 <- .pick(tab_f2, c(id_cols, "F","PValue","FDR"))
  colnames(main_f2)[-(1:3)] <- paste0("ANOVA_second_", colnames(main_f2)[-(1:3)])
  
  inter <- .pick(tab_int, c(id_cols, "F","PValue","FDR"))
  colnames(inter)[-(1:3)] <- paste0("ANOVA_interaction_", colnames(inter)[-(1:3)])
  
  posthoc_prefixed <- lapply(names(posthoc), function(nm){
    .with_prefix(posthoc[[nm]], paste0(nm, "_"), id_cols)
  })
  names(posthoc_prefixed) <- names(posthoc)
  
  combined <- Reduce(function(a,b) merge(a,b,by=id_cols,all=TRUE),
                     c(list(g, main_f1, main_f2, inter), unname(posthoc_prefixed)))
  rownames(combined) <- NULL
  .msg("[combine] Combined table: %d genes x %d columns.", nrow(combined), ncol(combined))
  
  .msg("[done] Two-factor edgeR ANOVA finished.")
  list(
    anova = list(
      global        = tab_global,
      main_first    = tab_f1,
      main_second   = tab_f2,
      interaction   = tab_int
    ),
    posthoc = posthoc,
    combined = combined,
    meta = list(
      first_factor  = list(name = first_factor_col,  levels = first_levels),
      second_factor = list(name = second_factor_col, levels = second_levels),
      robust = robust
    )
  )
}
