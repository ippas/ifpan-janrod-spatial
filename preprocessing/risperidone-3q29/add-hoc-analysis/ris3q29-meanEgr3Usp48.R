load(file = "results/risperidone-3q29/ris3q29_st_data.RData")

load(file = "results/risperidone-3q29/wtDel_summary_statistics.RData")

ris3q29_st_data$quantile_normalize_resolution_0.4$annotate$gene_name %>% unique %>% length()



imap_dfr(
  wtDel_summary_statistics$quantile_normalize$resolution_0.4,
  function(cluster_data, cluster_name) {
    cbind(
      cluster_data$control$mean,
      cluster_data$experiment$mean
    ) %>%
      as.data.frame() %>%
      rownames_to_column(var = "peak") %>%
      mutate(
        gene = sub(".*-", "", peak),
        cluster = cluster_name
      ) %>%
      select(cluster, peak, gene, everything())
  }
) %>% 
  filter(peak %in% c("peak-340087-Egr4",
                     "peak-294663-Usp48")) %>% 
  filter(cluster %in% c("cluster_0", "cluster_15")) %>% head 
  write.table(file = "results/risperidone-3q29/ad-hoc-results/Egr4Usp48-meanPerSample-quantileRes0.4-wtDel.csv", 
              quote = F,
              sep = ",",
              col.names = T,
              row.names = F)

  
  