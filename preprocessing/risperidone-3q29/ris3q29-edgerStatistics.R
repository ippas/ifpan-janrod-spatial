
# ##############################################################################
# ---- prepare data ----
# ##############################################################################

wtDel_summary_statistics$raw_data$resolution_0.1$cluster_1$control$sum 
wtDel_summary_statistics$raw_data$resolution_0.1$cluster_1$experiment$sum

cbind(
  wtDel_summary_statistics$raw_data$resolution_0.4$cluster_$control$sum,
  wtDel_summary_statistics$raw_data$resolution_0.4$cluste$experiment$sum
) -> data_sum_per_sample

metadata_ris3q29 %>% head

metadata_ris3q29$mouse_genotype %>% unique()

  metadata_ris3q29$treatment %>% unique()


# ##############################################################################
# ---- statistics ----
# ##############################################################################

edger_twofactor_anova(
  counts = data_sum_per_sample,                    # matrix: genes x samples (rownames=genes, colnames=samples)
  sample_meta = metadata_ris3q29,               # data.frame with columns: sample_col, first_factor_col, second_factor_col
  sample_col = "sample_ID",
  first_factor_col   = "mouse_genotype",
  second_factor_col  = "treatment",
  first_levels       = c("wtwt", "wtdel"),          # e.g., c("WT","MUT"); if NULL -> inferred from data (keeps order of appearance)
  second_levels      = c("saline", "risperidone"),          # e.g., c("CTRL","DEX","CORT")
  robust             = TRUE,
  verbose            = TRUE
)  -> tmp
  
  
tmp$combined %>% head
tmp$anova$global %>% head
  

res <- edger_twofactor_anova_and_posthoc(
  counts            = data_sum_per_sample,
  sample_meta       = metadata_ris3q29,
  sample_col        = "sample_ID",
  first_factor_col  = "mouse_genotype",           # wtwt / wtdel
  second_factor_col = "treatment",                # saline / risperidone
  first_levels      = c("wtwt","wtdel"),
  second_levels     = c("saline","risperidone"),
  robust            = TRUE,
  verbose           = TRUE,
  gatekeep          = T,
  gatekeep_fdr      = 0.05,
  gatekeep_log2FC   = 1,        # np. wymagaj co najmniej 0.5 w którymś z prostych efektów
  p_adjust_method   = "BH"
)


res_edger_MarginalProportional <- run_edger_twofactor_for_clusters_v2(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            = metadata_ris3q29,
  sample_col             = "sample_ID",
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  first_levels           = c("wtwt","wtdel"),
  second_levels          = c("saline","risperidone"),
  main_effects           = "marginal",
  main_effects_weights   = "proportional",   # <<< kluczowe
  clusters               = NULL
)

res$anova$global %>% 
  filter(FDR < 0.1)

res$anova$main_first

res$anova$interaction

res$anova$global %>% 
  filter(FDR < 0.05)



res_all <- run_edger_twofactor_for_clusters(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            =  metadata_ris3q29,
  sample_col             = "sample_ID",
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  first_levels           = c("wtwt","wtdel"),
  second_levels          = c("saline","risperidone"),
  groups                 = c("control","experiment"),
  count_slot             = "sum",
  clusters               = NULL,   # all clusters
  robust                 = TRUE,
  verbose                = TRUE,
  return_combined        = TRUE
)
res_all$results$cluster_19$anova$global
res_all$results %>% 
  lapply(., function(x){
    x$anova$global %>% filter(FDR < 0.1)
  })


res_all$results %>% 
  lapply(., function(x){
    x$anova$interaction %>% filter(PValue < 0.05)
  }) %>% write_cluster_list_to_xlsx(file_path = "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics/edger_interaction_res0.4_pval0.05.xlsx")

res_all$results %>% 
  lapply(., function(x){
    x$anova$global %>% filter(FDR < 0.05)
  }) %>% write_cluster_list_to_xlsx(file_path = "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics/edger_global_res0.4_fdr0.05.xlsx")

res_all$results %>% 
  imap_dfr(~ {
    .x$anova$global %>% 
      filter(PValue < 0.05) %>% 
      select(name_id, peak_id, gene_symbol) %>% 
      mutate(cluster = .y)
  }) %>% 
  group_by(name_id, gene_symbol, peak_id) %>% 
  summarise(
    clusters = paste(unique(cluster), collapse = ", "),
    n_clusters = n_distinct(cluster),
    .groups = "drop"
  ) %>% 
  arrange(desc(n_clusters)) %>% 
  select(!c(clusters)) 
    
library(openxlsx)

res_all$results %>% 
  imap_dfr(~ {
    .x$anova$global %>% 
      filter(FDR < 0.01) %>% 
      select(name_id, peak_id, gene_symbol) %>% 
      mutate(cluster = .y)
  }) %>% 
  group_by(name_id, gene_symbol, peak_id) %>% 
  summarise(
    clusters   = paste(unique(cluster), collapse = ", "),
    n_clusters = n_distinct(cluster),
    .groups = "drop"
  ) %>% 
  arrange(desc(n_clusters)) %>% 
  select(!clusters) %>% 
  { 
    list(
      df           = .,
      unique_genes = tibble(gene_symbol = unique(.$gene_symbol))
    )
  } %>% 
  { 
    wb <- createWorkbook()
    walk2(., names(.), ~ {
      addWorksheet(wb, .y)
      writeData(wb, .y, .x)
      freezePane(wb, sheet = .y, firstRow = TRUE)  # blokada pierwszego wiersza
    })
    saveWorkbook(wb, "results/risperidone-3q29/edger_statistics/combined_genesGlobalResultsAnovaEdger_fdr0.01.xlsx", overwrite = TRUE)
  }



res_all$results %>% 
  lapply(., function(x){
    x$anova$global %>% filter(PValue < 0.01) %>% nrow
  }) %>% unlist()

res_all$results %>% 
  imap_dfr(~ {
    .x$anova$global %>% 
      filter(FDR < 0.05) %>% 
      select(name_id, peak_id, gene_symbol) %>% 
      mutate(cluster = .y)
  }) %>% 
  group_by(name_id, gene_symbol, peak_id) %>% 
  summarise(
    clusters   = paste(unique(cluster), collapse = ", "),
    n_clusters = n_distinct(cluster),
    .groups = "drop"
  ) %>% 
  arrange(desc(n_clusters)) %>% 
  select(!clusters) %>% 
  { 
    list(
      df           = .,
      unique_genes = tibble(gene_symbol = unique(.$gene_symbol))
    )
  } %>% 
  { 
    wb <- createWorkbook()
    walk2(., names(.), ~ {
      addWorksheet(wb, .y)
      writeData(wb, .y, .x)
      freezePane(wb, sheet = .y, firstRow = TRUE)  # blokada pierwszego wiersza
    })
    saveWorkbook(wb, "results/risperidone-3q29/edger_statistics/combined_genesGlobalResultsAnovaEdger_fdr0.05.xlsx", overwrite = TRUE)
  }



res_all$results$cluster_0$anova$main_second

# wtSal vs wtRis 
# detSal vs delRis
# allSal vs allRis

res_all$results$cluster_0$anova$global %>%head

# ##############################################################################
# ----  ----
# ##############################################################################

df1 <- res_all$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% 
      filter(abs(logFC) > 0.8)
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_wtSalWtRis",
                 "logCPM_wtSalWtRis", "F_wtSalWtRis", "PValue_wtSalWtRis", "FDR_wtSalWtRis"))

df1 %>% filter(gene_symbol == "Arc")

df2 <- res_all$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` %>%
      filter(abs(logFC) < 0.2)
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_delSalDelRis",
                 "logCPM_delSalDelRis", "F_delSalDelRis", "PValue_delSalDelRis", "FDR_delSalDelRis"))

df2 %>% filter(gene_symbol == "Arc")


df3 <- res_all$results %>%
  lapply(., function(x){
    x$anova$main_second %>% filter(PValue < 0.01)
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_allSalAllRis",
                 "logCPM_allSalAllRis", "F_allSalAllRis", "PValue_allSalAllRis", "FDR_allSalAllRis"))

# --- joinowanie ---
df_joined <- df1 %>%
  inner_join(df2, by = c("cluster", "name_id", "peak_id", "gene_symbol")) %>%
  inner_join(df3, by = c("cluster", "name_id", "peak_id", "gene_symbol"))

df_joined %>% filter(gene_symbol == "Arc")


df_joined %>% 
  mutate(check_logFC = ifelse(logFC_wtSalWtRis == logFC_allSalAllRis, T, F)) %>% 
  .$check_logFC %>% table


write.xlsx(
  df_joined,
  file = "results/risperidone-3q29/edger_statistics/edger_wtSalWtRisLog2ratioGreater0.8_delSalDelRisLog2ratioLess0.2_allSalAllDelPvalue0.01.xlsx"
)
# ##############################################################################
# ----  ----
# ##############################################################################
df1 <- res_all$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_wtSalWtRis",
                 "logCPM_wtSalWtRis", "F_wtSalWtRis", "PValue_wtSalWtRis", "FDR_wtSalWtRis"))

df1 %>% filter(gene_symbol == "Arc")

df2 <- res_all$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_delSalDelRis",
                 "logCPM_delSalDelRis", "F_delSalDelRis", "PValue_delSalDelRis", "FDR_delSalDelRis"))

df2 %>% filter(gene_symbol == "Arc")


df3 <- res_all$results %>%
  lapply(., function(x){
    x$anova$main_second 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_allSalAllRis",
                 "logCPM_allSalAllRis", "F_allSalAllRis", "PValue_allSalAllRis", "FDR_allSalAllRis"))

# --- joinowanie ---
df_joined <- df1 %>%
  inner_join(df2, by = c("cluster", "name_id", "peak_id", "gene_symbol")) %>%
  inner_join(df3, by = c("cluster", "name_id", "peak_id", "gene_symbol"))

df_joined %>% 
  filter(PValue_allSalAllRis < 0.01) %>% 
  head
  



df_joined %>% 
  filter(PValue_allSalAllRis < 0.01) %>%
  ggplot(aes(x = logFC_delSalDelRis, y = logFC_wtSalWtRis)) +
  geom_point(alpha = 0.6, size = 0.8) +
  # linie pomocnicze
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = c(-1.5, 1.5), linetype = "dashed", color = "red", linewidth = 0.6) +
  # etykiety genów (|x|>1.2 lub |y|>1.8)
  geom_text(
    data = df_joined %>%
      filter(PValue_allSalAllRis < 0.01) %>%
      filter(abs(logFC_delSalDelRis) > 1.2 | abs(logFC_wtSalWtRis) > 1.8),
    aes(label = gene_symbol),
    vjust = -0.4, size = 2.8, color = "black"
  )
# jeśli wolisz bez nakładania etykiet:
# ggrepel::geom_text_repel(aes(label = gene_symbol), size = 2.8, max.overlaps = Inf, min.segment.length = 0)
+
  labs(
    x = "logFC delSalDelRis",
    y = "logFC wtSalWtRis",
    title = "Porównanie logFC między warunkami (facety po klastrach)"
  ) +
  xlim(-3, 3.5) +
  ylim(-3, 3.5) +
  facet_wrap(~ cluster, ncol = 4, scales = "fixed") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

write.xlsx(
  df_joined,
  file = "results/risperidone-3q29/edger_statistics/edger_wtSalWtRisLog2ratioLess0.2_delSalDelRisLog2ratioGreather0.8_allSalAllDelPvalue0.01.xlsx"
)



library(dplyr)
library(ggplot2)

df_points <- df_joined %>%
  filter(PValue_allSalAllRis < 0.01) %>%
  mutate(
    cluster = factor(cluster, levels = paste0("cluster_", 0:19)) # kolejność numeryczna
  )

df_labels <- df_points %>%
  filter(abs(logFC_delSalDelRis) > 1.5 | abs(logFC_wtSalWtRis) > 1.5)

ggplot(df_points, aes(x = logFC_delSalDelRis, y = logFC_wtSalWtRis)) +
  geom_point(alpha = 0.6, size = 0.8) +
  # szare linie w 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # czerwone linie na ±0.5
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "red", linewidth = 0.6) +
  # etykiety genów
  geom_text(
    data = df_labels,
    aes(x = logFC_delSalDelRis, y = logFC_wtSalWtRis, label = gene_symbol),
    vjust = -0.4, size = 2.8, color = "black", inherit.aes = FALSE
  ) +
  labs(
    x = "logFC delSalDelRis",
    y = "logFC wtSalWtRis",
    title = "Porównanie logFC między warunkami (klastry)"
  ) +
  # xlim(-3, 3.5) +
  # ylim(-3, 3.5) +
  facet_wrap(~ cluster, ncol = 5, scales = "fixed") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


df_joined %>% 
  filter(abs(logFC_wtSalWtRis) > 0.8 & PValue_wtSalWtRis < 0.01) %>% nrow

df_joined %>% 
  filter(abs(logFC_delSalDelRis) > 0.8 & PValue_delSalDelRis < 0.01) %>% nrow
