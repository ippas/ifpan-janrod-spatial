# ==============================================================================
# Libraries
# ==============================================================================

library(dplyr)
library(purrr)
library(openxlsx)
library(preprocessCore)


# ==============================================================================
# 1. Define samples included in the analysis
# ==============================================================================

samples_ris3q29 <- metadata_ris3q29 %>%
  filter(sample_ID != "S13839Nr3") %>%
  pull(sample_ID)


# ==============================================================================
# 2. Calculate pseudobulk raw sum and mean expression per cluster
# ==============================================================================

ris3q29_n25Samples_rawRes04_min15Spots_sumMean <- compute_cluster_sum_mean_v3(
  spatial_data = ris3q29_st_data,
  resolution = 0.4,
  samples = samples_ris3q29,
  data_type = "raw_data",
  min_number_spots = 15,
  num_cores = 24,
  verbose = TRUE
)


# Inspect raw summed counts for cluster 0
ris3q29_n25Samples_rawRes04_min15Spots_sumMean$cluster_0$sum


# Inspect the complete object structure
ris3q29_n25Samples_rawRes04_min15Spots_sumMean %>%
  str()


# ==============================================================================
# 3. Run edgeR two-factor analysis separately for each cluster
# ==============================================================================

ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors <- 
  run_edger_twoFactors_perClusters(
    cluster_summary_data = ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
    
    sample_meta = metadata_ris3q29 %>%
      filter(sample_ID != "S13839Nr3"),
    
    sample_col = "sample_ID",
    genotype_col = "mouse_genotype",
    treatment_col = "treatment",
    
    genotype_levels = c(
      "wtwt",
      "wtdel"
    ),
    
    treatment_levels = c(
      "saline",
      "risperidone"
    ),
    
    cell_labels = c(
      "salWt",
      "salDel",
      "risWt",
      "risDel"
    ),
    
    count_slot = "sum",
    robust = TRUE,
    verbose = TRUE
  )


# ==============================================================================
# 4. Prepare sample names with experimental-group prefixes
# ==============================================================================

sample_names <- metadata_ris3q29 %>%
  filter(sample_ID != "S13839Nr3") %>%
  transmute(
    sample_ID,
    
    prefixed_sample_id = paste0(
      case_when(
        mouse_genotype == "wtwt"  & treatment == "saline"      ~ "salWt",
        mouse_genotype == "wtdel" & treatment == "saline"      ~ "salDel",
        mouse_genotype == "wtwt"  & treatment == "risperidone" ~ "risWt",
        mouse_genotype == "wtdel" & treatment == "risperidone" ~ "risDel"
      ),
      "_",
      sample_ID
    )
  )


# ==============================================================================
# 5. Export raw summed counts per cluster
# ==============================================================================

ris3q29_n25Samples_rawRes04_min15Spots_allClusters_sum <- 
  imap_dfr(
    ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
    
    ~ {
      counts <- as.data.frame(
        .x$sum,
        check.names = FALSE
      )
      
      colnames(counts) <- sample_names$prefixed_sample_id[
        match(
          colnames(counts),
          sample_names$sample_ID
        )
      ]
      
      cbind(
        data.frame(
          name_id = rownames(.x$sum),
          peak = sub("-[^-]+$", "", rownames(.x$sum)),
          gene = .x$gene,
          cluster = .y,
          check.names = FALSE
        ),
        counts
      )
    }
  )

rownames(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_sum
) <- NULL


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_sum,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_rawSumCountsPerCluster.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 6. Export raw mean expression per spot and cluster
# ==============================================================================

ris3q29_n25Samples_rawRes04_min15Spots_allClusters_mean <- 
  imap_dfr(
    ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
    
    ~ {
      mean_expression <- as.data.frame(
        .x$mean,
        check.names = FALSE
      )
      
      colnames(mean_expression) <- sample_names$prefixed_sample_id[
        match(
          colnames(mean_expression),
          sample_names$sample_ID
        )
      ]
      
      cbind(
        data.frame(
          name_id = rownames(.x$mean),
          peak = sub("-[^-]+$", "", rownames(.x$mean)),
          gene = .x$gene,
          cluster = .y,
          check.names = FALSE
        ),
        mean_expression
      )
    }
  )

rownames(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_mean
) <- NULL


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_mean,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_rawMeanCountsPerSpot.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 7. Log2-transform and quantile-normalize raw summed counts within each cluster
# ==============================================================================

ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum <- 
  imap_dfr(
    ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
    
    ~ {
      current_values <- .x$sum
      
      valid_samples <- colSums(!is.na(current_values)) > 0
      
      normalized_values <- matrix(
        NA_real_,
        nrow = nrow(current_values),
        ncol = ncol(current_values),
        dimnames = dimnames(current_values)
      )
      
      if (sum(valid_samples) > 0) {
        
        log2_values <- log2(
          current_values[, valid_samples, drop = FALSE] + 1
        )
        
        normalized_values[, valid_samples] <- 
          preprocessCore::normalize.quantiles(log2_values)
      }
      
      colnames(normalized_values) <- sample_names$prefixed_sample_id[
        match(
          colnames(normalized_values),
          sample_names$sample_ID
        )
      ]
      
      cbind(
        data.frame(
          name_id = rownames(normalized_values),
          peak = sub("-[^-]+$", "", rownames(normalized_values)),
          gene = .x$gene,
          cluster = .y,
          check.names = FALSE
        ),
        as.data.frame(
          normalized_values,
          check.names = FALSE
        )
      )
    }
  )

rownames(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum
) <- NULL


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_log2QuantileNormalizedRawSumCountsPerCluster.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 8. Log2-transform and quantile-normalize raw mean counts within each cluster
# ==============================================================================

ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean <- 
  imap_dfr(
    ris3q29_n25Samples_rawRes04_min15Spots_sumMean,
    
    ~ {
      current_values <- .x$mean
      
      valid_samples <- colSums(!is.na(current_values)) > 0
      
      normalized_values <- matrix(
        NA_real_,
        nrow = nrow(current_values),
        ncol = ncol(current_values),
        dimnames = dimnames(current_values)
      )
      
      if (sum(valid_samples) > 0) {
        
        log2_values <- log2(
          current_values[, valid_samples, drop = FALSE] + 1
        )
        
        normalized_values[, valid_samples] <- 
          preprocessCore::normalize.quantiles(log2_values)
      }
      
      colnames(normalized_values) <- sample_names$prefixed_sample_id[
        match(
          colnames(normalized_values),
          sample_names$sample_ID
        )
      ]
      
      cbind(
        data.frame(
          name_id = rownames(normalized_values),
          peak = sub("-[^-]+$", "", rownames(normalized_values)),
          gene = .x$gene,
          cluster = .y,
          check.names = FALSE
        ),
        as.data.frame(
          normalized_values,
          check.names = FALSE
        )
      )
    }
  )

rownames(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean
) <- NULL


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_log2QuantileNormalizedRawMeanCountsPerSpot.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 9. Calculate group means for log2 quantile-normalized summed counts
# ==============================================================================

salWt_columns <- grep(
  "^salWt_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum),
  value = TRUE
)

salDel_columns <- grep(
  "^salDel_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum),
  value = TRUE
)

risWt_columns <- grep(
  "^risWt_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum),
  value = TRUE
)

risDel_columns <- grep(
  "^risDel_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum),
  value = TRUE
)


ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum_groupMeans <-
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum %>%
  transmute(
    name_id,
    peak,
    gene,
    cluster,
    
    mean_salWt = rowMeans(
      as.matrix(pick(all_of(salWt_columns))),
      na.rm = TRUE
    ),
    
    mean_salDel = rowMeans(
      as.matrix(pick(all_of(salDel_columns))),
      na.rm = TRUE
    ),
    
    mean_risWt = rowMeans(
      as.matrix(pick(all_of(risWt_columns))),
      na.rm = TRUE
    ),
    
    mean_risDel = rowMeans(
      as.matrix(pick(all_of(risDel_columns))),
      na.rm = TRUE
    )
  ) %>%
  mutate(
    across(
      starts_with("mean_"),
      ~ replace(.x, is.nan(.x), NA_real_)
    )
  )


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNsum_groupMeans,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_log2QuantileNormalizedRawSumCountsPerCluster_groupMeans.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 10. Calculate group means for log2 quantile-normalized mean counts per spot
# ==============================================================================

salWt_columns <- grep(
  "^salWt_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean),
  value = TRUE
)

salDel_columns <- grep(
  "^salDel_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean),
  value = TRUE
)

risWt_columns <- grep(
  "^risWt_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean),
  value = TRUE
)

risDel_columns <- grep(
  "^risDel_",
  colnames(ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean),
  value = TRUE
)


ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean_groupMeans <-
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean %>%
  transmute(
    name_id,
    peak,
    gene,
    cluster,
    
    mean_salWt = rowMeans(
      as.matrix(pick(all_of(salWt_columns))),
      na.rm = TRUE
    ),
    
    mean_salDel = rowMeans(
      as.matrix(pick(all_of(salDel_columns))),
      na.rm = TRUE
    ),
    
    mean_risWt = rowMeans(
      as.matrix(pick(all_of(risWt_columns))),
      na.rm = TRUE
    ),
    
    mean_risDel = rowMeans(
      as.matrix(pick(all_of(risDel_columns))),
      na.rm = TRUE
    )
  ) %>%
  mutate(
    across(
      starts_with("mean_"),
      ~ replace(.x, is.nan(.x), NA_real_)
    )
  )


openxlsx::write.xlsx(
  ris3q29_n25Samples_rawRes04_min15Spots_allClusters_log2QNmean_groupMeans,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_log2QuantileNormalizedRawMeanCountsPerSpot_groupMeans.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 11. Inspect edgeR results after corrected contrast names
# ==============================================================================

names(
  ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors$results$cluster_0$full_results
)

ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors$contrast_definitions

ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors$results$cluster_0$full_results$salDel_vs_salWt %>%
  head()


# ==============================================================================
# 12. Combine edgeR full results from all clusters
# ==============================================================================

comparison_prefixes <- c(
  interaction_genotype_x_treatment = "interaction",
  genotype_main_del_vs_wt = "genotype_main_del_vs_wt",
  treatment_main_ris_vs_sal = "treatment_main_ris_vs_sal",
  salDel_vs_salWt = "salDel_vs_salWt",
  risDel_vs_risWt = "risDel_vs_risWt",
  risWt_vs_salWt = "risWt_vs_salWt",
  risDel_vs_salDel = "risDel_vs_salDel",
  risDel_vs_salWt = "risDel_vs_salWt",
  risWt_vs_salDel = "risWt_vs_salDel"
)

statistics_columns <- c(
  "log2FC",
  "logCPM",
  "F",
  "p_value",
  "fdr"
)


ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors_combinedResults_df <- 
  ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors$results %>%
  
  purrr::imap_dfr(function(cluster_result, cluster_name) {
    
    full_results <- cluster_result$full_results
    
    if (is.null(full_results) || length(full_results) == 0) {
      return(NULL)
    }
    
    full_results <- full_results %>%
      
      purrr::imap(function(result_table, comparison_name) {
        
        if (is.null(result_table) || nrow(result_table) == 0) {
          return(NULL)
        }
        
        prefix <- comparison_prefixes[[comparison_name]]
        
        if (is.null(prefix)) {
          stop(
            "Missing prefix for comparison: ",
            comparison_name
          )
        }
        
        result_table <- as.data.frame(
          result_table,
          check.names = FALSE
        )
        
        names(result_table)[
          names(result_table) %in% statistics_columns
        ] <- paste0(
          prefix,
          "_",
          names(result_table)[
            names(result_table) %in% statistics_columns
          ]
        )
        
        result_table
      }) %>%
      
      purrr::compact()
    
    if (length(full_results) == 0) {
      return(NULL)
    }
    
    full_results %>%
      
      purrr::reduce(
        dplyr::full_join,
        by = c("name_id", "peak", "gene")
      ) %>%
      
      dplyr::mutate(
        cluster = cluster_name,
        .before = name_id
      )
  })


ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors_combinedResults_df %>%
  head()


# ==============================================================================
# 13. Export combined edgeR full results from all clusters
# ==============================================================================

openxlsx::write.xlsx(
  ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors_combinedResults_df,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "allClusters_edgeRTwoFactorsFullResults.xlsx"
  ),
  
  rowNames = FALSE,
  overwrite = TRUE
)


# ==============================================================================
# 14. Save complete edgeR two-factor analysis object as RData
# ==============================================================================

save(
  ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors,
  
  file = paste0(
    "/home/mateusz/projects/ifpan-janrod-spatial/results/",
    "risperidone-3q29/edger_statistics/",
    "pseudobulkRawSum_edgerTwoFactors_n25SamplesRes0.4min15Spots/",
    "ris3q29_n25SamplesRes0.4min15Spots_",
    "edgeRTwoFactorsFullStatistics.RData"
  )
)


# ==============================================================================
# 15. Quick final check of corrected columns
# ==============================================================================

ris3q29_n25Samples_rawSumRes0.4min15spots_edgerTwoFactors_combinedResults_df %>%
  select(
    cluster,
    name_id,
    peak,
    gene,
    salDel_vs_salWt_log2FC,
    risDel_vs_risWt_log2FC,
    risWt_vs_salWt_log2FC,
    risDel_vs_salDel_log2FC,
    risDel_vs_salWt_log2FC,
    risWt_vs_salDel_log2FC
  ) %>%
  head()