# ============================================================
# edgeR two-factor pseudobulk analysis per spatial cluster
#
# Input:
#   - cluster_summary_data:
#       output of compute_cluster_sum_mean_v3()
#       each cluster must contain:
#         $sum       raw pseudobulk counts: features x samples
#         $peak      original feature identifiers
#         $gene      gene symbols
#         $sample_qc spot-level QC table
#
#   - sample_meta:
#       metadata with sample ID, genotype, and treatment columns
#
# Output:
#   $results$cluster_X$full_results
#       nine full edgeR result tables:
#         interaction_genotype_x_treatment
#         genotype_main_wt_vs_del
#         treatment_main_sal_vs_ris
#         salWt_vs_salDel
#         risWt_vs_risDel
#         salWt_vs_risWt
#         salDel_vs_risDel
#         salWt_vs_risDel
#         salDel_vs_risWt
#
#   Each result table contains:
#     name_id, peak, gene, log2FC, logCPM, F, p_value, fdr
#
#   Per-cluster filtering output:
#     $preprocessing_summary
#     $features_removed_filterByExpr
#     $features_retained_filterByExpr
#     $genes_fully_removed_filterByExpr
#     $genes_partially_removed_filterByExpr
#     $genes_retained_filterByExpr
#     $removed_annotation_filterByExpr
#     $retained_annotation_filterByExpr
#     $fully_removed_annotation_filterByExpr
#     $removed_features_from_partially_retained_genes
#     $filterByExpr_parameters
# ============================================================

library(dplyr)
library(edgeR)

run_edger_twoFactors_perClusters <- function(
  cluster_summary_data,
  sample_meta,
  sample_col = "sample_ID",
  genotype_col = "mouse_genotype",
  treatment_col = "treatment",
  genotype_levels = c("wtwt", "wtdel"),
  treatment_levels = c("saline", "risperidone"),
  cell_labels = c("salWt", "salDel", "risWt", "risDel"),
  count_slot = "sum",
  robust = TRUE,
  verbose = TRUE
) {
  
  # ============================================================
  # Helper: print messages only when verbose = TRUE
  # ============================================================
  
  .msg <- function(...) {
    if (isTRUE(verbose)) {
      message(sprintf(...))
    }
  }
  
  # ============================================================
  # Helper: remove only the terminal "-GENE" suffix from a peak ID
  #
  # Example:
  #   peak-341618-Slc6a6  +  Slc6a6
  #   -> peak-341618
  # ============================================================
  
  clean_peak_ids <- function(peak_values, gene_values) {
    
    mapply(
      function(current_peak, current_gene) {
        
        if (is.na(current_peak)) {
          return(NA_character_)
        }
        
        current_peak <- as.character(current_peak)
        
        if (
          is.na(current_gene) ||
          !nzchar(current_gene)
        ) {
          return(current_peak)
        }
        
        suffix <- paste0("-", as.character(current_gene))
        
        if (endsWith(current_peak, suffix)) {
          return(
            substr(
              current_peak,
              1,
              nchar(current_peak) - nchar(suffix)
            )
          )
        }
        
        return(current_peak)
      },
      peak_values,
      gene_values,
      USE.NAMES = FALSE
    )
  }
  
  # ============================================================
  # Basic input checks
  # ============================================================
  
  if (!is.list(cluster_summary_data)) {
    stop("cluster_summary_data must be a named list.")
  }
  
  if (!is.data.frame(sample_meta)) {
    sample_meta <- as.data.frame(sample_meta)
  }
  
  if (length(genotype_levels) != 2) {
    stop("genotype_levels must contain exactly two values.")
  }
  
  if (length(treatment_levels) != 2) {
    stop("treatment_levels must contain exactly two values.")
  }
  
  if (length(cell_labels) != 4) {
    stop(
      "cell_labels must contain exactly four labels in this order: ",
      "salWt, salDel, risWt, risDel."
    )
  }
  
  required_meta_columns <- c(
    sample_col,
    genotype_col,
    treatment_col
  )
  
  missing_meta_columns <- setdiff(
    required_meta_columns,
    colnames(sample_meta)
  )
  
  if (length(missing_meta_columns) > 0) {
    stop(
      "Missing columns in sample_meta: ",
      paste(missing_meta_columns, collapse = ", ")
    )
  }
  
  all_clusters <- names(cluster_summary_data)
  
  if (
    is.null(all_clusters) ||
    length(all_clusters) == 0
  ) {
    stop("cluster_summary_data has no named clusters.")
  }
  
  # ============================================================
  # Prepare sample metadata
  # ============================================================
  
  genotype_levels <- trimws(as.character(genotype_levels))
  treatment_levels <- trimws(as.character(treatment_levels))
  
  sample_meta <- sample_meta %>%
    mutate(
      .sample_id = as.character(.data[[sample_col]]),
      .genotype = trimws(as.character(.data[[genotype_col]])),
      .treatment = trimws(as.character(.data[[treatment_col]]))
    )
  
  if (anyDuplicated(sample_meta$.sample_id)) {
    duplicated_ids <- sample_meta$.sample_id[
      duplicated(sample_meta$.sample_id)
    ]
    
    stop(
      "Duplicated sample IDs in sample_meta: ",
      paste(unique(duplicated_ids), collapse = ", ")
    )
  }
  
  observed_genotypes <- sort(unique(sample_meta$.genotype))
  observed_treatments <- sort(unique(sample_meta$.treatment))
  
  if (!all(genotype_levels %in% observed_genotypes)) {
    stop(
      "genotype_levels do not match metadata. Observed values: ",
      paste(observed_genotypes, collapse = ", ")
    )
  }
  
  if (!all(treatment_levels %in% observed_treatments)) {
    stop(
      "treatment_levels do not match metadata. Observed values: ",
      paste(observed_treatments, collapse = ", ")
    )
  }
  
  # ============================================================
  # Define four factorial groups
  #
  # salWt  = saline + wtwt
  # salDel = saline + wtdel
  # risWt  = risperidone + wtwt
  # risDel = risperidone + wtdel
  # ============================================================
  
  cell_map <- data.frame(
    .genotype = rep(genotype_levels, times = 2),
    .treatment = rep(treatment_levels, each = 2),
    cell_group = cell_labels,
    stringsAsFactors = FALSE
  )
  
  sample_meta <- sample_meta %>%
    left_join(
      cell_map,
      by = c(".genotype", ".treatment")
    ) %>%
    mutate(
      cell_group = factor(
        cell_group,
        levels = cell_labels
      )
    )
  
  if (any(is.na(sample_meta$cell_group))) {
    stop(
      "Could not map all samples to the four factorial groups."
    )
  }
  
  # ============================================================
  # Define contrasts
  # Model columns are always:
  # salWt | salDel | risWt | risDel
  # ============================================================
  
  sal_wt <- cell_labels[1]
  sal_del <- cell_labels[2]
  ris_wt <- cell_labels[3]
  ris_del <- cell_labels[4]
  
  make_contrast <- function(values) {
    stats::setNames(values, cell_labels)
  }
  
  contrast_vectors <- list(
    
    interaction_genotype_x_treatment = make_contrast(
      c(1, -1, -1, 1)
    ),
    
    genotype_main_wt_vs_del = make_contrast(
      c(0.5, -0.5, 0.5, -0.5)
    ),
    
    treatment_main_sal_vs_ris = make_contrast(
      c(0.5, 0.5, -0.5, -0.5)
    ),
    
    salWt_vs_salDel = make_contrast(
      c(1, -1, 0, 0)
    ),
    
    risWt_vs_risDel = make_contrast(
      c(0, 0, 1, -1)
    ),
    
    salWt_vs_risWt = make_contrast(
      c(1, 0, -1, 0)
    ),
    
    salDel_vs_risDel = make_contrast(
      c(0, 1, 0, -1)
    ),
    
    salWt_vs_risDel = make_contrast(
      c(1, 0, 0, -1)
    ),
    
    salDel_vs_risWt = make_contrast(
      c(0, 1, -1, 0)
    )
  )
  
  contrast_definitions <- data.frame(
    test = names(contrast_vectors),
    test_type = c(
      "interaction",
      "main_effect_genotype",
      "main_effect_treatment",
      rep("pairwise", 6)
    ),
    contrast = c(
      paste0(
        "(", ris_del, " - ", ris_wt, ") - ",
        "(", sal_del, " - ", sal_wt, ")"
      ),
      paste0(
        "0.5*(", sal_wt, " - ", sal_del, ") + ",
        "0.5*(", ris_wt, " - ", ris_del, ")"
      ),
      paste0(
        "0.5*(", sal_wt, " + ", sal_del, ") - ",
        "0.5*(", ris_wt, " + ", ris_del, ")"
      ),
      paste0(sal_wt, " - ", sal_del),
      paste0(ris_wt, " - ", ris_del),
      paste0(sal_wt, " - ", ris_wt),
      paste0(sal_del, " - ", ris_del),
      paste0(sal_wt, " - ", ris_del),
      paste0(sal_del, " - ", ris_wt)
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================
  # Run edgeR for one cluster
  # ============================================================
  
  process_one_cluster <- function(cluster_name) {
    
    tryCatch({
      
      .msg("\n============================================================")
      .msg("[cluster %s] Starting.", cluster_name)
      
      cluster_data <- cluster_summary_data[[cluster_name]]
      
      if (inherits(cluster_data, "cluster_error")) {
        stop(
          "Upstream summary error: ",
          cluster_data$error_message
        )
      }
      
      if (is.null(cluster_data[[count_slot]])) {
        stop(
          "Slot '", count_slot,
          "' is missing in ", cluster_name
        )
      }
      
      counts_input <- cluster_data[[count_slot]]
      
      if (!is.matrix(counts_input)) {
        counts_input <- as.matrix(counts_input)
      }
      
      if (!is.numeric(counts_input)) {
        stop("Count matrix must be numeric.")
      }
      
      if (is.null(rownames(counts_input))) {
        stop("Count matrix has no feature rownames.")
      }
      
      if (is.null(colnames(counts_input))) {
        stop("Count matrix has no sample column names.")
      }
      
      if (anyDuplicated(rownames(counts_input))) {
        stop("Duplicated feature IDs in count matrix.")
      }
      
      if (anyDuplicated(colnames(counts_input))) {
        stop("Duplicated sample IDs in count matrix.")
      }
      
      if (any(counts_input < 0, na.rm = TRUE)) {
        stop("Negative values detected in raw counts.")
      }
      
      .msg(
        "[cluster %s] Input: %d features × %d pseudobulks.",
        cluster_name,
        nrow(counts_input),
        ncol(counts_input)
      )
      
      # ----------------------------------------------------------
      # Feature annotation
      # ----------------------------------------------------------
      
      peak_input <- if (
        !is.null(cluster_data$peak) &&
        length(cluster_data$peak) == nrow(counts_input)
      ) {
        as.character(cluster_data$peak)
      } else {
        rownames(counts_input)
      }
      
      gene_input <- if (
        !is.null(cluster_data$gene) &&
        length(cluster_data$gene) == nrow(counts_input)
      ) {
        as.character(cluster_data$gene)
      } else {
        rownames(counts_input)
      }
      
      feature_annotation <- data.frame(
        name_id = rownames(counts_input),
        peak = clean_peak_ids(
          peak_values = peak_input,
          gene_values = gene_input
        ),
        gene = gene_input,
        stringsAsFactors = FALSE
      )
      
      # ----------------------------------------------------------
      # Align metadata to count matrix columns
      # ----------------------------------------------------------
      
      metadata_index <- match(
        colnames(counts_input),
        sample_meta$.sample_id
      )
      
      if (anyNA(metadata_index)) {
        missing_metadata <- colnames(counts_input)[
          is.na(metadata_index)
        ]
        
        stop(
          "Missing metadata for samples: ",
          paste(missing_metadata, collapse = ", ")
        )
      }
      
      smeta <- sample_meta[
        metadata_index,
        ,
        drop = FALSE
      ]
      
      # ----------------------------------------------------------
      # Sample-level QC inherited from cluster-summary object
      # ----------------------------------------------------------
      
      sample_qc <- smeta %>%
        transmute(
          sample_id = .sample_id,
          genotype = .genotype,
          treatment = .treatment,
          cell_group = as.character(cell_group)
        )
      
      if (
        !is.null(cluster_data$sample_qc) &&
        is.data.frame(cluster_data$sample_qc) &&
        "sample" %in% colnames(cluster_data$sample_qc)
      ) {
        
        cluster_sample_qc <- cluster_data$sample_qc %>%
          transmute(
            sample_id = as.character(sample),
            n_spots = as.integer(n_spots),
            pass_min_spots = as.logical(pass_min_spots)
          )
        
        sample_qc <- sample_qc %>%
          left_join(
            cluster_sample_qc,
            by = "sample_id"
          )
        
      } else {
        
        sample_qc <- sample_qc %>%
          mutate(
            n_spots = NA_integer_,
            pass_min_spots = NA
          )
      }
      
      # ----------------------------------------------------------
      # Remove pseudobulks that are unusable for edgeR
      #
      # All-NA pseudobulks occur when sample × cluster has fewer
      # spots than min_number_spots in the upstream summarization.
      # ----------------------------------------------------------
      
      all_na_columns <- colSums(!is.na(counts_input)) == 0
      
      partial_na_columns <- !all_na_columns &
        colSums(is.na(counts_input)) > 0
      
      if (any(partial_na_columns)) {
        stop(
          "Partial NA values found in samples: ",
          paste(
            colnames(counts_input)[partial_na_columns],
            collapse = ", "
          )
        )
      }
      
      input_library_size <- colSums(
        counts_input,
        na.rm = TRUE
      )
      
      zero_library_columns <- !all_na_columns &
        input_library_size == 0
      
      keep_sample_columns <- !all_na_columns &
        !zero_library_columns
      
      sample_qc <- sample_qc %>%
        mutate(
          input_library_size = as.numeric(input_library_size),
          count_status = case_when(
            all_na_columns ~ "dropped_all_na",
            zero_library_columns ~ "dropped_zero_library",
            TRUE ~ "used_in_edgeR"
          )
        )
      
      .msg(
        "[cluster %s] Pseudobulks removed: all-NA = %d; zero-library = %d.",
        cluster_name,
        sum(all_na_columns),
        sum(zero_library_columns)
      )
      
      if (!any(keep_sample_columns)) {
        stop("No pseudobulks remain after sample filtering.")
      }
      
      counts <- counts_input[
        ,
        keep_sample_columns,
        drop = FALSE
      ]
      
      smeta_used <- smeta[
        keep_sample_columns,
        ,
        drop = FALSE
      ]
      
      smeta_used$cell_group <- factor(
        smeta_used$cell_group,
        levels = cell_labels
      )
      
      # ----------------------------------------------------------
      # Group counts after pseudobulk filtering
      # ----------------------------------------------------------
      
      group_counts <- data.frame(
        cell_group = cell_labels,
        n_samples = as.integer(
          table(
            factor(
              smeta_used$cell_group,
              levels = cell_labels
            )
          )
        ),
        stringsAsFactors = FALSE
      )
      
      if (any(group_counts$n_samples == 0)) {
        missing_groups <- group_counts %>%
          filter(n_samples == 0) %>%
          pull(cell_group)
        
        stop(
          "No usable pseudobulk in group(s): ",
          paste(missing_groups, collapse = ", ")
        )
      }
      
      .msg(
        "[cluster %s] Group sizes: %s",
        cluster_name,
        paste(
          paste0(
            group_counts$cell_group,
            "=n",
            group_counts$n_samples
          ),
          collapse = "; "
        )
      )
      
      # ----------------------------------------------------------
      # Four-group edgeR design
      # ----------------------------------------------------------
      
      design <- model.matrix(
        ~ 0 + cell_group,
        data = smeta_used
      )
      
      colnames(design) <- cell_labels
      
      if (qr(design)$rank < ncol(design)) {
        stop("Design matrix is not full rank.")
      }
      
      if (nrow(design) <= qr(design)$rank) {
        stop(
          "No residual degrees of freedom remain after fitting the model."
        )
      }
      
      # ----------------------------------------------------------
      # edgeR preprocessing and QL model fitting
      # ----------------------------------------------------------
      
      .msg(
        "[cluster %s] Running filterByExpr() on raw pseudobulk counts.",
        cluster_name
      )
      
      y <- edgeR::DGEList(
        counts = counts
      )
      
      # ----------------------------------------------------------
      # filterByExpr() parameters and cluster-specific threshold
      #
      # The filtering rule is not a single fixed raw-count cutoff.
      # edgeR calculates a CPM cutoff from the median raw library
      # size of the current cluster and requires expression in a
      # design-dependent minimum number of pseudobulks.
      # ----------------------------------------------------------
      
      filter_min_count <- 10
      filter_min_total_count <- 15
      filter_large_n <- 10
      filter_min_prop <- 0.7
      
      raw_library_sizes_before_filter <- y$samples$lib.size
      median_library_size_before_filter <- median(
        raw_library_sizes_before_filter
      )
      
      design_q <- qr.Q(
        qr(design)
      )
      
      hat_values <- rowSums(
        design_q ^ 2
      )
      
      min_samples_required_raw <- 1 / max(hat_values)
      
      min_samples_required <- min_samples_required_raw
      
      if (min_samples_required > filter_large_n) {
        min_samples_required <- filter_large_n +
          (min_samples_required - filter_large_n) *
          filter_min_prop
      }
      
      cpm_threshold <- filter_min_count /
        median_library_size_before_filter *
        1e6
      
      filterByExpr_parameters <- data.frame(
        cluster = cluster_name,
        min_count = filter_min_count,
        min_total_count = filter_min_total_count,
        large_n = filter_large_n,
        min_prop = filter_min_prop,
        median_library_size_before_filter =
          median_library_size_before_filter,
        min_library_size_before_filter =
          min(raw_library_sizes_before_filter),
        max_library_size_before_filter =
          max(raw_library_sizes_before_filter),
        cpm_threshold = cpm_threshold,
        min_samples_required_raw =
          min_samples_required_raw,
        min_samples_required =
          min_samples_required,
        min_samples_required_ceiling =
          ceiling(min_samples_required),
        stringsAsFactors = FALSE
      )
      
      keep_genes <- edgeR::filterByExpr(
        y,
        design = design,
        min.count = filter_min_count,
        min.total.count = filter_min_total_count,
        large.n = filter_large_n,
        min.prop = filter_min_prop
      )
      
      if (sum(keep_genes) == 0) {
        stop("filterByExpr() removed all features.")
      }
      
      .msg(
        paste0(
          "[cluster %s] filterByExpr(): CPM >= %.4f in >= %.2f ",
          "pseudobulks (rounded up: %d); total counts >= %d."
        ),
        cluster_name,
        cpm_threshold,
        min_samples_required,
        ceiling(min_samples_required),
        filter_min_total_count
      )
      
      .msg(
        "[cluster %s] filterByExpr(): retained %d / %d features.",
        cluster_name,
        sum(keep_genes),
        length(keep_genes)
      )
      
      y <- y[
        keep_genes,
        ,
        keep.lib.sizes = FALSE
      ]
      
      .msg(
        "[cluster %s] Calculating TMM normalization factors.",
        cluster_name
      )
      
      y <- edgeR::calcNormFactors(
        y,
        method = "TMM"
      )
      
      .msg(
        "[cluster %s] Estimating dispersions and fitting QL model.",
        cluster_name
      )
      
      y <- edgeR::estimateDisp(
        y,
        design = design,
        robust = robust
      )
      
      fit <- edgeR::glmQLFit(
        y,
        design = design,
        robust = robust
      )
      
      # ----------------------------------------------------------
      # Add TMM / library information to sample QC
      # ----------------------------------------------------------
      
      normalization_qc <- data.frame(
        sample_id = colnames(y$counts),
        library_size_after_filterByExpr = y$samples$lib.size,
        tmm_norm_factor = y$samples$norm.factors,
        effective_library_size =
          y$samples$lib.size * y$samples$norm.factors,
        stringsAsFactors = FALSE
      )
      
      sample_qc <- sample_qc %>%
        left_join(
          normalization_qc,
          by = "sample_id"
        )
      
      # ----------------------------------------------------------
      # Features retained and removed by filterByExpr()
      # ----------------------------------------------------------
      
      annotation_retained <- feature_annotation[
        keep_genes,
        ,
        drop = FALSE
      ]
      
      annotation_removed <- feature_annotation[
        !keep_genes,
        ,
        drop = FALSE
      ]
      
      # ----------------------------------------------------------
      # Gene-level interpretation of filterByExpr()
      #
      # A gene is fully removed only when none of its features /
      # peaks passed filterByExpr().
      #
      # A gene is partially removed when at least one feature was
      # removed but at least one other feature for the same gene
      # was retained.
      # ----------------------------------------------------------
      
      genes_input <- unique(feature_annotation$gene)
      genes_retained <- unique(annotation_retained$gene)
      genes_removed_any_peak <- unique(annotation_removed$gene)
      
      genes_fully_removed <- setdiff(
        genes_input,
        genes_retained
      )
      
      genes_partially_removed <- intersect(
        genes_removed_any_peak,
        genes_retained
      )
      
      fully_removed_annotation <- annotation_removed %>%
        filter(gene %in% genes_fully_removed)
      
      removed_features_from_partially_retained_genes <- annotation_removed %>%
        filter(gene %in% genes_partially_removed)
      
      # Reorder retained annotation to match the edgeR count matrix
      annotation_tested <- annotation_retained[
        match(
          rownames(y$counts),
          annotation_retained$name_id
        ),
        ,
        drop = FALSE
      ]
      
      # ----------------------------------------------------------
      # Prepare one concise full result table
      # ----------------------------------------------------------
      
      make_result_table <- function(qlf_test) {
        
        result_table <- as.data.frame(qlf_test$table)
        
        data.frame(
          name_id = annotation_tested$name_id,
          peak = annotation_tested$peak,
          gene = annotation_tested$gene,
          log2FC = result_table$logFC,
          logCPM = result_table$logCPM,
          F = result_table$F,
          p_value = result_table$PValue,
          fdr = p.adjust(
            result_table$PValue,
            method = "BH"
          ),
          stringsAsFactors = FALSE
        ) %>%
          arrange(
            fdr,
            p_value
          )
      }
      
      # ----------------------------------------------------------
      # All nine tests use the same filtered and fitted model
      # ----------------------------------------------------------
      
      .msg(
        "[cluster %s] Calculating interaction, main effects and six contrasts.",
        cluster_name
      )
      
      full_results <- lapply(
        contrast_vectors,
        function(current_contrast) {
          
          qlf_test <- edgeR::glmQLFTest(
            fit,
            contrast = current_contrast
          )
          
          make_result_table(qlf_test)
        }
      )
      
      # ----------------------------------------------------------
      # Compact preprocessing summary
      # ----------------------------------------------------------
      
      preprocessing_summary <- data.frame(
        cluster = cluster_name,
        status = "completed",
        n_pseudobulks_input = ncol(counts_input),
        n_pseudobulks_below_min_spots = sum(
          !sample_qc$pass_min_spots,
          na.rm = TRUE
        ),
        n_pseudobulks_removed_all_na = sum(all_na_columns),
        n_pseudobulks_removed_zero_library = sum(zero_library_columns),
        n_pseudobulks_used_edgeR = ncol(counts),
        n_features_input = nrow(counts_input),
        n_features_removed_filterByExpr = nrow(annotation_removed),
        n_features_retained_filterByExpr = nrow(annotation_retained),
        filter_min_count = filter_min_count,
        filter_min_total_count = filter_min_total_count,
        filter_median_library_size = median_library_size_before_filter,
        filter_cpm_threshold = cpm_threshold,
        filter_min_samples_required = min_samples_required,
        filter_min_samples_required_ceiling =
          ceiling(min_samples_required),
        n_unique_genes_input = length(genes_input),
        n_unique_genes_fully_removed_filterByExpr =
          length(genes_fully_removed),
        n_unique_genes_partially_removed_filterByExpr =
          length(genes_partially_removed),
        n_unique_genes_retained_filterByExpr =
          length(genes_retained),
        stringsAsFactors = FALSE
      )
      
      .msg("[cluster %s] Finished successfully.", cluster_name)
      
      list(
        status = "completed",
        full_results = full_results,
        preprocessing_summary = preprocessing_summary,
        # Feature-level filtering output
        features_removed_filterByExpr = annotation_removed$name_id,
        features_retained_filterByExpr = annotation_retained$name_id,
        removed_annotation_filterByExpr = annotation_removed,
        retained_annotation_filterByExpr = annotation_retained,
        
        # Gene-level filtering output
        # fully removed = no feature / peak for this gene was retained
        genes_fully_removed_filterByExpr = genes_fully_removed,
        
        # partially removed = some features removed, at least one retained
        genes_partially_removed_filterByExpr = genes_partially_removed,
        
        # retained = at least one feature / peak for this gene was retained
        genes_retained_filterByExpr = genes_retained,
        
        # Annotation tables supporting gene-level interpretation
        fully_removed_annotation_filterByExpr =
          fully_removed_annotation,
        
        removed_features_from_partially_retained_genes =
          removed_features_from_partially_retained_genes,
        
        filterByExpr_parameters = filterByExpr_parameters,
        sample_qc = sample_qc,
        group_counts = group_counts,
        design = design
      )
      
    }, error = function(e) {
      
      .msg(
        "[cluster %s] FAILED: %s",
        cluster_name,
        e$message
      )
      
      list(
        status = "failed",
        error_message = e$message,
        preprocessing_summary = data.frame(
          cluster = cluster_name,
          status = "failed",
          error_message = e$message,
          stringsAsFactors = FALSE
        )
      )
    })
  }
  
  # ============================================================
  # Run sequentially for all clusters
  # ============================================================
  
  .msg("============================================================")
  .msg("Starting edgeR two-factor analysis.")
  .msg("Clusters: %d", length(all_clusters))
  .msg("Count slot: %s", count_slot)
  .msg("Groups: %s", paste(cell_labels, collapse = ", "))
  .msg("============================================================")
  
  start_time <- Sys.time()
  
  results <- lapply(
    all_clusters,
    process_one_cluster
  )
  
  names(results) <- all_clusters
  
  completed_clusters <- names(results)[
    vapply(
      results,
      function(x) x$status == "completed",
      logical(1)
    )
  ]
  
  failed_clusters <- names(results)[
    vapply(
      results,
      function(x) x$status == "failed",
      logical(1)
    )
  ]
  
  runtime_seconds <- round(
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )
    ),
    1
  )
  
  .msg("============================================================")
  .msg("Finished in %s seconds.", runtime_seconds)
  .msg("Completed clusters: %d", length(completed_clusters))
  .msg("Failed clusters: %d", length(failed_clusters))
  
  if (length(failed_clusters) > 0) {
    .msg(
      "Failed: %s",
      paste(failed_clusters, collapse = ", ")
    )
  }
  
  .msg("============================================================")
  
  list(
    results = results,
    contrast_definitions = contrast_definitions,
    contrast_vectors = contrast_vectors,
    meta = list(
      sample_col = sample_col,
      genotype_col = genotype_col,
      treatment_col = treatment_col,
      genotype_levels = genotype_levels,
      treatment_levels = treatment_levels,
      cell_labels = cell_labels,
      count_slot = count_slot,
      robust = robust,
      completed_clusters = completed_clusters,
      failed_clusters = failed_clusters,
      runtime_seconds = runtime_seconds
    )
  )
}

# ============================================================
# Example call
# ============================================================

# ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors <-
#   run_edger_twoFactors_perClusters(
#     cluster_summary_data = ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
#
#     sample_meta = metadata_ris3q29 %>%
#       filter(sample_ID != "S13839Nr3"),
#
#     sample_col = "sample_ID",
#     genotype_col = "mouse_genotype",
#     treatment_col = "treatment",
#
#     genotype_levels = c("wtwt", "wtdel"),
#     treatment_levels = c("saline", "risperidone"),
#
#     cell_labels = c(
#       "salWt",
#       "salDel",
#       "risWt",
#       "risDel"
#     ),
#
#     count_slot = "sum",
#     robust = TRUE,
#     verbose = TRUE
#   )
