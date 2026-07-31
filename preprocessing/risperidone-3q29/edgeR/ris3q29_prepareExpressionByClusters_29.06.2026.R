# ##############################################################################
# Pseudobulk SUM counts by cluster
#
# Input:
#   wtDel_summary_statistics$raw_data$resolution_0.4
#   metadata_ris3q29
#
# Important:
#   Only sample-cluster combinations with >=20 spots are available in the
#   original SUM matrices.
#   Missing sample-cluster combinations remain complete NA columns.
#
# Final files:
#   n25Samples_min20SpotsPerCluster_rawPseudobulkSumCounts_allClusters_2026-06-29.xlsx
#   n25Samples_min20SpotsPerCluster_log2QuantileNormalizedPseudobulkSumCounts_allClusters_2026-06-29.xlsx
#
# Final objects retained:
#   raw_sum_all_clusters
#   log2_qn_sum_all_clusters
# ##############################################################################


# ==============================================================================
# 1. Libraries and paths
# ==============================================================================

library(dplyr)

output_dir <- "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics"


# ==============================================================================
# 2. Helper: parse peak and gene IDs
# ==============================================================================

parse_ids <- function(feature_ids) {
  
  tokens <- strsplit(feature_ids, "[ -]+", perl = TRUE)
  
  data.frame(
    name_id = feature_ids,
    
    peak_id = vapply(tokens, function(x) {
      if (length(x) >= 2) {
        paste(x[1:2], collapse = "-")
      } else {
        x[1]
      }
    }, character(1)),
    
    gene_symbol = vapply(tokens, function(x) {
      x[length(x)]
    }, character(1)),
    
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. Create raw pseudobulk SUM tables: one data frame per cluster
# ==============================================================================

make_raw_sum_by_cluster <- function(
  raw_data_by_resolution,
  groups = c("control", "experiment"),
  count_slot = "sum"
) {
  
  result <- vector(
    mode = "list",
    length = length(raw_data_by_resolution)
  )
  
  names(result) <- names(raw_data_by_resolution)
  
  for (cluster_name in names(raw_data_by_resolution)) {
    
    message("[cluster] ", cluster_name)
    
    cluster_object <- raw_data_by_resolution[[cluster_name]]
    
    missing_groups <- setdiff(
      groups,
      names(cluster_object)
    )
    
    if (length(missing_groups) > 0) {
      stop(
        "Missing groups in cluster ",
        cluster_name,
        ": ",
        paste(missing_groups, collapse = ", ")
      )
    }
    
    count_list <- lapply(groups, function(group_name) {
      
      counts <- cluster_object[[group_name]][[count_slot]]
      
      if (is.null(counts)) {
        stop(
          "Missing ",
          count_slot,
          " counts in cluster ",
          cluster_name,
          " / group ",
          group_name
        )
      }
      
      counts <- as.matrix(counts)
      
      if (is.null(rownames(counts)) || is.null(colnames(counts))) {
        stop(
          "Counts matrix lacks rownames or colnames in cluster ",
          cluster_name,
          " / group ",
          group_name
        )
      }
      
      counts
    })
    
    reference_features <- rownames(count_list[[1]])
    
    count_list <- lapply(count_list, function(counts) {
      
      if (!setequal(reference_features, rownames(counts))) {
        stop(
          "Different feature sets between groups in cluster: ",
          cluster_name
        )
      }
      
      counts[
        reference_features,
        ,
        drop = FALSE
      ]
    })
    
    counts_combined <- do.call(
      cbind,
      count_list
    )
    
    if (anyDuplicated(colnames(counts_combined))) {
      stop(
        "Duplicated sample IDs in cluster: ",
        cluster_name
      )
    }
    
    result[[cluster_name]] <- cbind(
      parse_ids(rownames(counts_combined)),
      as.data.frame(
        counts_combined,
        check.names = FALSE
      )
    )
  }
  
  result
}


# ==============================================================================
# 4. Create raw SUM count tables per cluster
# ==============================================================================

raw_sum_by_cluster <- make_raw_sum_by_cluster(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  groups = c("control", "experiment"),
  count_slot = "sum"
)


# ==============================================================================
# 5. Identify the real 25 samples present in raw data
# ==============================================================================

annotation_cols <- c(
  "name_id",
  "peak_id",
  "gene_symbol"
)

sample_ids_present_in_raw <- unique(
  unlist(
    lapply(raw_sum_by_cluster, function(df) {
      setdiff(
        colnames(df),
        annotation_cols
      )
    })
  )
)

metadata_sample_ids <- as.character(
  metadata_ris3q29$sample_ID
)

unknown_raw_samples <- setdiff(
  sample_ids_present_in_raw,
  metadata_sample_ids
)

if (length(unknown_raw_samples) > 0) {
  stop(
    "Raw data contain sample IDs absent from metadata: ",
    paste(unknown_raw_samples, collapse = ", ")
  )
}

# Keep metadata order, but only samples actually present in the raw matrices.
all_sample_ids <- metadata_sample_ids[
  metadata_sample_ids %in% sample_ids_present_in_raw
]

if (length(all_sample_ids) != 25L) {
  stop(
    "Expected 25 real samples in raw matrices, found: ",
    length(all_sample_ids)
  )
}

if (anyDuplicated(all_sample_ids)) {
  stop("Duplicated sample IDs detected.")
}


# ==============================================================================
# 6. Ensure every cluster has all 25 sample columns
# Missing sample-cluster combinations become complete NA columns
# ==============================================================================

missing_samples_by_cluster <- lapply(
  raw_sum_by_cluster,
  function(df) {
    setdiff(
      all_sample_ids,
      colnames(df)
    )
  }
)

if (any(lengths(missing_samples_by_cluster) > 0)) {
  
  message("\nSamples added as complete NA columns:")
  
  print(
    missing_samples_by_cluster[
      lengths(missing_samples_by_cluster) > 0
    ]
  )
}

raw_sum_by_cluster <- lapply(raw_sum_by_cluster, function(df) {
  
  missing_samples <- setdiff(
    all_sample_ids,
    colnames(df)
  )
  
  for (sample_id in missing_samples) {
    df[[sample_id]] <- NA_real_
  }
  
  df[
    ,
    c(
      annotation_cols,
      all_sample_ids
    ),
    drop = FALSE
  ]
})


# ==============================================================================
# 7. Log2(count + 1) + quantile normalization per cluster
# Complete NA sample columns remain NA
# ==============================================================================

log2_qn_one_cluster <- function(df) {
  
  count_matrix <- as.matrix(
    df[
      ,
      all_sample_ids,
      drop = FALSE
    ]
  )
  
  storage.mode(count_matrix) <- "numeric"
  
  # TRUE only for samples having any valid value in this cluster.
  valid_samples <- colSums(
    !is.na(count_matrix)
  ) > 0
  
  # A valid sample should not contain partial NA values.
  partial_na_samples <- (
    colSums(is.na(count_matrix)) > 0
  ) & valid_samples
  
  if (any(partial_na_samples)) {
    stop(
      "Unexpected partial NA values in samples: ",
      paste(
        colnames(count_matrix)[partial_na_samples],
        collapse = ", "
      )
    )
  }
  
  normalized_matrix <- matrix(
    NA_real_,
    nrow = nrow(count_matrix),
    ncol = ncol(count_matrix),
    dimnames = list(
      rownames(count_matrix),
      colnames(count_matrix)
    )
  )
  
  # Quantile normalization only for samples with data.
  # Entirely missing sample columns stay NA.
  if (any(valid_samples)) {
    
    normalized_valid_matrix <- preprocessCore::normalize.quantiles(
      log2(
        count_matrix[
          ,
          valid_samples,
          drop = FALSE
        ] + 1
      )
    )
    
    colnames(normalized_valid_matrix) <- colnames(
      count_matrix[
        ,
        valid_samples,
        drop = FALSE
      ]
    )
    
    normalized_matrix[
      ,
      valid_samples
    ] <- normalized_valid_matrix
  }
  
  cbind(
    df[
      ,
      annotation_cols,
      drop = FALSE
    ],
    
    as.data.frame(
      normalized_matrix,
      check.names = FALSE
    )
  )
}


# ==============================================================================
# 8. Create normalized tables per cluster
# ==============================================================================

log2_qn_sum_by_cluster <- lapply(
  raw_sum_by_cluster,
  log2_qn_one_cluster
)


# ==============================================================================
# 9. Combine all clusters into two final data frames
# ==============================================================================

raw_sum_all_clusters <- dplyr::bind_rows(
  raw_sum_by_cluster,
  .id = "cluster"
)

log2_qn_sum_all_clusters <- dplyr::bind_rows(
  log2_qn_sum_by_cluster,
  .id = "cluster"
)

rownames(raw_sum_all_clusters) <- NULL
rownames(log2_qn_sum_all_clusters) <- NULL


# ==============================================================================
# 10. Check that both final tables have 25 sample columns
# ==============================================================================

technical_cols <- c(
  "cluster",
  annotation_cols
)

raw_sample_cols <- setdiff(
  colnames(raw_sum_all_clusters),
  technical_cols
)

normalized_sample_cols <- setdiff(
  colnames(log2_qn_sum_all_clusters),
  technical_cols
)

if (length(raw_sample_cols) != 25L) {
  stop(
    "Raw table has ",
    length(raw_sample_cols),
    " sample columns instead of 25."
  )
}

if (length(normalized_sample_cols) != 25L) {
  stop(
    "Normalized table has ",
    length(normalized_sample_cols),
    " sample columns instead of 25."
  )
}

if (!identical(raw_sample_cols, all_sample_ids)) {
  stop("Raw sample column order differs from metadata order.")
}

if (!identical(normalized_sample_cols, all_sample_ids)) {
  stop("Normalized sample column order differs from metadata order.")
}


# ==============================================================================
# 11. Create readable sample names
# sal_wt_ID / sal_del_ID / ris_wt_ID / ris_del_ID
# ==============================================================================

sample_rename <- metadata_ris3q29 %>%
  dplyr::filter(
    sample_ID %in% all_sample_ids
  ) %>%
  dplyr::transmute(
    sample_ID = as.character(sample_ID),
    
    treatment_short = dplyr::recode(
      treatment,
      "saline" = "sal",
      "risperidone" = "ris",
      .default = NA_character_
    ),
    
    genotype_short = dplyr::recode(
      mouse_genotype,
      "wtwt" = "wt",
      "wtdel" = "del",
      .default = NA_character_
    ),
    
    new_name = paste(
      treatment_short,
      genotype_short,
      sample_ID,
      sep = "_"
    )
  )

if (anyNA(sample_rename$new_name)) {
  stop(
    "Unknown treatment or mouse_genotype in metadata."
  )
}

if (anyDuplicated(sample_rename$new_name)) {
  stop("Duplicated renamed sample names were created.")
}

rename_lookup <- stats::setNames(
  sample_rename$new_name,
  sample_rename$sample_ID
)


# ==============================================================================
# 12. Rename sample columns in both final tables
# ==============================================================================

rename_sample_columns <- function(df, rename_lookup) {
  
  new_colnames <- colnames(df)
  
  matched <- match(
    new_colnames,
    names(rename_lookup)
  )
  
  to_rename <- !is.na(matched)
  
  new_colnames[to_rename] <- unname(
    rename_lookup[
      matched[to_rename]
    ]
  )
  
  colnames(df) <- new_colnames
  
  df
}

raw_sum_all_clusters <- rename_sample_columns(
  raw_sum_all_clusters,
  rename_lookup
)

log2_qn_sum_all_clusters <- rename_sample_columns(
  log2_qn_sum_all_clusters,
  rename_lookup
)


# ==============================================================================
# 13. Save final Excel files
# ==============================================================================

raw_output_file <- file.path(
  output_dir,
  "n25Samples_min20SpotsPerCluster_rawPseudobulkSumCounts_allClusters_2026-06-29.xlsx"
)

log2_qn_output_file <- file.path(
  output_dir,
  "n25Samples_min20SpotsPerCluster_log2QuantileNormalizedPseudobulkSumCounts_allClusters_2026-06-29.xlsx"
)

openxlsx::write.xlsx(
  raw_sum_all_clusters,
  file = raw_output_file,
  overwrite = TRUE
)

openxlsx::write.xlsx(
  log2_qn_sum_all_clusters,
  file = log2_qn_output_file,
  overwrite = TRUE
)

message("\nSaved raw counts:")
message(raw_output_file)

message("\nSaved log2 + quantile-normalized counts:")
message(log2_qn_output_file)


# ==============================================================================
# 14. Remove temporary objects
# ==============================================================================

rm(
  list = intersect(
    c(
      "output_dir",
      "annotation_cols",
      "sample_ids_present_in_raw",
      "metadata_sample_ids",
      "unknown_raw_samples",
      "all_sample_ids",
      "missing_samples_by_cluster",
      "raw_sum_by_cluster",
      "log2_qn_sum_by_cluster",
      "technical_cols",
      "raw_sample_cols",
      "normalized_sample_cols",
      "sample_rename",
      "rename_lookup",
      "raw_output_file",
      "log2_qn_output_file",
      "parse_ids",
      "make_raw_sum_by_cluster",
      "log2_qn_one_cluster",
      "rename_sample_columns"
    ),
    ls()
  )
)