res_edger_MarginalProportional

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


logcpm_res0.4 <- prepare_logCPM_per_sample_for_clusters(
  raw_data_by_resolution  = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta             = metadata_ris3q29,
  sample_col              = "sample_ID",
  groups                  = c("control","experiment"),  # zostaw jak w Twojej strukturze
  count_slot              = "sum",
  clusters                = NULL,        # wszystkie
  filter_by_expr          = TRUE,        # odsiej nisko-ekspresyjne (jak w edgeR)
  zscore_rows             = FALSE,       # na surowym logCPM (bez Z-score)
  return_format           = "list",
  verbose                 = TRUE
)

logcpm_res0.4$logCPM_by_cluster$cluster_0

# ##############################################################################
# ---- select interaction ----
# ##############################################################################
res_edger_MarginalProportional$results %>% 
  lapply(., function(x){
    x$anova$interaction %>% filter(PValue < 0.05)
  }) 

res_with_all <- res_edger_MarginalProportional$results %>% 
  lapply(function(x){
    x$anova$interaction %>% filter(PValue < 0.05)
  })


res_with_all$all_clusters %>% 
  filter(FDR < 0.1)


# dodajemy wersję scaloną
res_with_all$all_clusters <- bind_rows(res_with_all, .id = "cluster")
res_with_all %>% write_cluster_list_to_xlsx(file_path = "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/edger_statistics/edgerMarginalProportion_interaction_res0.4_pval0.05.xlsx")


# ##############################################################################
# ---- combine data ----
# ##############################################################################


df1 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_wtSalWtRis",
                 "logCPM_wtSalWtRis", "F_wtSalWtRis", "PValue_wtSalWtRis", "FDR_wtSalWtRis"))

df1 %>% filter(gene_symbol == "Arc")

df2 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_delSalDelRis",
                 "logCPM_delSalDelRis", "F_delSalDelRis", "PValue_delSalDelRis", "FDR_delSalDelRis"))

df2 %>% filter(gene_symbol == "Arc")


df3 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$anova$main_second 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_allSalAllRis",
                 "logCPM_allSalAllRis", "F_allSalAllRis", "PValue_allSalAllRis", "FDR_allSalAllRis"))

# --- joinowanie ---
df_joined <- df1 %>%
  inner_join(df2, by = c("cluster", "name_id", "peak_id", "gene_symbol")) %>%
  inner_join(df3, by = c("cluster", "name_id", "peak_id", "gene_symbol"))


# ##############################################################################
# ---- LABEL ----
# ##############################################################################

df_joined %>% 
  filter(abs(logFC_delSalDelRis) < 0.2 & abs(logFC_wtSalWtRis) > 0.8 & PValue_allSalAllRis < 0.01) %>% 
  dplyr::rename(logCPM = "logCPM_wtSalWtRis") %>% 
  select(!c(logCPM_delSalDelRis, logCPM_allSalAllRis)) %>%
  write.xlsx(
    file = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_wtSalWtRisLog2ratioGreater0.8_delSalDelRisLog2ratioLess0.2_allSalAllDelPvalue0.01.xlsx"
  )
  
df_joined %>% 
  filter(abs(logFC_delSalDelRis) > 0.8 & abs(logFC_wtSalWtRis) < 0.2 & PValue_allSalAllRis < 0.01) %>% 
  dplyr::rename(logCPM = "logCPM_wtSalWtRis") %>% 
  select(!c(logCPM_delSalDelRis, logCPM_allSalAllRis)) %>% 
  write.xlsx(
    file = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_wtSalWtRisLog2ratioLess0.2_delSalDelRisLog2ratioGreater0.8_allSalAllDelPvalue0.01.xlsx"
  )

# ##############################################################################
# ---- visualization ----
# ##############################################################################
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
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "red", linewidth = 0.6) +
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


# ##############################################################################
# ---- second combined data ----
# ##############################################################################


df1 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_wtSalWtRis",
                 "logCPM_wtSalWtRis", "F_wtSalWtRis", "PValue_wtSalWtRis", "FDR_wtSalWtRis"))

df2 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logFC_delSalDelRis",
                 "logCPM_delSalDelRis", "F_delSalDelRis", "PValue_delSalDelRis", "FDR_delSalDelRis"))

df3 <- res_edger_MarginalProportional$results %>%
  lapply(., function(x){
    x$anova$interaction 
  }) %>% bind_rows(.id = "cluster") %>%
  set_colnames(c("cluster", "name_id", "peak_id", "gene_symbol", "logCPM", "F_interaction", "PValue_interaction", "FDR_interaction"))

# --- joinowanie ---
df_joined <- df1 %>%
  inner_join(df2, by = c("cluster", "name_id", "peak_id", "gene_symbol")) %>%
  inner_join(df3, by = c("cluster", "name_id", "peak_id", "gene_symbol"))


df_joined %>% 
  select(-c(logCPM_delSalDelRis, logCPM_wtSalWtRis)) %>% 
  filter(PValue_interaction < 0.01) %>% 
  filter(abs(logFC_delSalDelRis) < 0.2) %>% 
  filter(abs(logFC_wtSalWtRis) > 0.8) %>% 
  write.xlsx(
    file = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_wtSalWtRisLog2ratioGreater0.8_delSalDelRisLog2ratioLess0.2_interactionPvalue0.01.xlsx"
  )


df_joined %>% 
  select(-c(logCPM_delSalDelRis, logCPM_wtSalWtRis)) %>% 
  filter(PValue_interaction < 0.01) %>% 
  filter(abs(logFC_delSalDelRis) > 0.8) %>% 
  filter(abs(logFC_wtSalWtRis) < 0.2) %>% 
  write.xlsx(
    file = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_wtSalWtRisLog2ratioLessr0.2_delSalDelRisLog2ratioGreater0.8_interactionPvalue0.01.xlsx"
  )



df_joined %>% 
  select(-c(logCPM_delSalDelRis, logCPM_wtSalWtRis)) %>% 
  filter(PValue_interaction < 0.01) %>% 
  filter(abs(logFC_delSalDelRis) > 0.8) %>%
  filter(abs(logFC_wtSalWtRis) < 0.2) %>% 
  write.xlsx(
    file = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_wtSalWtRisLog2ratio_delSalDelRisLog2ratio_interactionPvalue0.01.xlsx"
  )
  
  
  # Wtsal vs wtRis
  # DelSal vs DelRis
  # WtSal vs DelSal
  # WtRis vs DelRis
  # AllWt vs AllDel
  # i interakcja
  # [12:31]
  # w .csv i  .xlsx
  
res_edger_MarginalProportional$results$cluster_0$posthoc$`SECOND_risperidone-saline|FIRST_wtwt`

combine_with_clusters <- function(results, 
                                  slot = c("anova", "posthoc"), 
                                  element = NULL,
                                  filter_expr = NULL) {
  slot <- match.arg(slot)
  
  tmp <- lapply(results, function(x){
    # sięgnij do results[[cluster]][[slot]][[element]]
    obj <- x[[slot]]
    if (!is.null(element)) {
      if (!element %in% names(obj)) {
        stop(paste("Element", element, "not found in slot", slot))
      }
      obj <- obj[[element]]
    }
    
    # opcjonalne filtrowanie
    if (!is.null(filter_expr)) {
      obj <- dplyr::filter(obj, !!rlang::parse_expr(filter_expr))
    }
    
    obj
  })
  
  combined_data <- dplyr::bind_rows(tmp, .id = "cluster")
  
  # do listy wyników dokładamy combined
  tmp$all_clusters <- combined_data
  
  tmp
}

combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "anova", element = "interaction") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_interactionPValue0.05.R")

combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "anova", element = "main_first") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_allWtAllDelPValue0.05.R")

combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "posthoc", element = "SECOND_risperidone-saline|FIRST_wtwt") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_salWtRisWt_PValue0.05.R")

combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "posthoc", element = "SECOND_risperidone-saline|FIRST_wtdel") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_salDelRisDel_PValue0.05.R")


combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "posthoc", element = "FIRST_wtdel-wtwt|SECOND_risperidone") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_risWtRisDel_PValue0.05.R")

combine_with_clusters(res_edger_MarginalProportional$results,
                      slot = "posthoc", element = "FIRST_wtdel-wtwt|SECOND_saline") %>% 
  write_cluster_list_to_xlsx(file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_salWtSisDel_PValue0.05.xlsx")



# ##############################################################################
# ---- merge data to one table ----
# ##############################################################################


# --- klucze wspólne
keys <- c("cluster", "name_id", "peak_id", "gene_symbol")

t_global <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "anova", element = "global") %>% 
  .$all_cluster %>%
  select(-logCPM) %>% 
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",  "logFC.mouse_genotypewtdel_global", "logFC.treatmentrisperidone_global",
                 "logFC.mouse_genotypewtdel.treatmentrisperidone_global",
                 "F_global","P_global","FDR_global"))

t_interaction <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "anova", element = "interaction") %>% 
  .$all_cluster %>% 
  set_colnames(c("cluster","name_id","peak_id","gene_symbol","logCPM",
                 "F_interaction","P_interaction","FDR_interaction"))

t_main_first <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "anova", element = "main_first") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>% 
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_allWtAllDel","F_allWtAllDel","P_allWtAllDel","FDR_allWtAllDel"))

t_main_second <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "anova", element = "main_second") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>% 
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_allSalAllRis","F_allSalAllRis","P_allSalAllRis","FDR_allSalAllRis"))

t_posthoc_salRisWt <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "posthoc", element = "SECOND_risperidone-saline|FIRST_wtwt") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>%  
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_salRisWt","F_salRisWt","P_salRisWt","FDR_salRisWt"))

t_posthoc_salRisDel <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "posthoc", element = "SECOND_risperidone-saline|FIRST_wtdel") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>%  
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_salRisDel","F_salRisDel","P_salRisDel","FDR_salRisDel"))

t_posthoc_risWtDel <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "posthoc", element = "FIRST_wtdel-wtwt|SECOND_risperidone") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>%  
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_risWtDel","F_risWtDel","P_risWtDel","FDR_risWtDel"))

t_posthoc_salWtDel <-
  combine_with_clusters(res_edger_MarginalProportional$results,
                        slot = "posthoc", element = "FIRST_wtdel-wtwt|SECOND_saline") %>% 
  .$all_cluster %>% 
  select(-logCPM) %>%  
  # poprawka literówki: FDR_selWtDel -> FDR_salWtDel
  set_colnames(c("cluster","name_id","peak_id","gene_symbol",
                 "logFC_salWtDel","F_salWtDel","P_salWtDel","FDR_salWtDel"))

# --- sklejanie w jedną dużą tabelę
big_tbl <- reduce(
  list(t_global, t_interaction, t_main_first, t_main_second, t_posthoc_salRisWt,
       t_posthoc_salRisDel, t_posthoc_risWtDel, t_posthoc_salWtDel),
  full_join, by = keys
) %>% as.data.frame() %>% 
  # filter(P_interaction < 0.05) %>% 
  arrange(cluster, name_id, peak_id, gene_symbol) 

list(combined_data = big_tbl) %>% write_cluster_list_to_xlsx(., file_path = "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_09.09.2025.xlsx")

rm(
  t_interaction,
  t_main_first,
  t_posthoc_salRisWt,
  t_posthoc_salRisDel,
  t_posthoc_risWtDel,
  t_posthoc_salWtDel
)

sample(c(1:nrow(t_posthoc_salRisWt)), 10)

t_main_second %>% 
  mutate(cluster = factor(cluster, 
                          levels = paste0("cluster_", 0:19))) %>% 
  ggplot(aes(x = cluster, y = logFC_allSalAllRis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size = 0.1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of log2ratio per cluster; comparision salWt vs risWt",
    x = "Cluster",
    y = "log2ratio (salRisWt)"
  )

  

anova_res <- aov(logFC_salRisWt ~ cluster, data = t_posthoc_salRisWt %>% 
                   filter(cluster %in% c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", 
                                         "cluster_6", "cluster_7", "cluster_8", "cluster_9", "cluster_10", "cluster_11")) %>% sample(1000, replace = T))
summary(anova_res)

summary_stats <- t_posthoc_salRisWt %>%
  group_by(cluster) %>%
  summarise(
    n       = n(),
    min     = min(logFC_salRisWt, na.rm = TRUE),
    Q1      = quantile(logFC_salRisWt, 0.25, na.rm = TRUE),
    median  = median(logFC_salRisWt, na.rm = TRUE),
    Q3      = quantile(logFC_salRisWt, 0.75, na.rm = TRUE),
    max     = max(logFC_salRisWt, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

logCPM_ris3q29 <- prepare_logCPM_per_sample_for_clusters(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            = metadata_ris3q29,
  sample_col             = "sample_ID",
  groups                 = c("control","experiment"),  # <<< ważne
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  clusters               = NULL,
  filter_by_expr         = TRUE,
  prior.count            = 2,
  return_format          = "long",
  verbose                = TRUE
)

#6.6 cluster_16
logCPM_ris3q29 %>%
  filter(gene_symbol == "Arc") %>% 
  filter(cluster == "cluster_16") %>% .$logCPM %>% mean


res_edger_v2 <- run_edger_twofactor_for_clusters_v2(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            = metadata_ris3q29,
  sample_col             = "sample_ID",
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  first_levels           = c("wtwt","wtdel"),
  second_levels          = c("saline","risperidone"),
  groups                 = c("control","experiment"),
  main_effects           = "marginal",
  main_effects_weights   = "proportional",
  return_preprocessing_data = T    # <<< nowa opcja
)


ris3q29_bundle_data_to_visualization <- res_edger_v2$results %>%
  # 1) surowe macierze
  map(~ as.matrix(.x$counts)) %>%
  # 2) zbuduj listę z raw + qn + metadata
  { list(
    raw_counts_by_cluster = .,
    qn_counts_by_cluster  = map(., ~ {
      m <- preprocessCore::normalize.quantiles(as.matrix(.x))
      dimnames(m) <- dimnames(.x)
      m
    }),
    metadata = metadata_ris3q29
  )
  }

ris3q29_bundle_data_to_visualization


plot_gene_boxplot_from_bundle <- function(
  bundle,
  gene_symbol,
  cluster,
  data_type = c("raw","qn"),
  log2_transform = FALSE,
  pseudocount = 1,
  group_cols = c("treatment","mouse_genotype"),
  title = NULL,
  y_limits = NULL,
  levels = NULL,     # kolejność grup
  palette = NULL     # named vector z kolorami
){
  
  #' Boxplot dla wybranego genu z bundla ris3q29
  #'
  #' Funkcja tworzy wykres typu boxplot z opcjonalnym log2-transform, 
  #' z wyborem danych raw lub quantile-normalized, dla podanego genu i klastra.
  #' Obsługuje sytuację, gdy ten sam gene_symbol ma kilka rekordów (wykres facet_wrap).
  #'
  #' @param bundle list; obiekt np. `ris3q29_bundle_data_to_visualization`, zawierający 
  #'        `raw_counts_by_cluster`, `qn_counts_by_cluster` oraz `metadata`.
  #' @param gene_symbol character; np. `"Arc"`, symbol genu do wizualizacji.
  #' @param cluster character; np. `"cluster_16"`, nazwa klastra.
  #' @param data_type character; `"raw"` lub `"qn"` (raw counts albo quantile normalized).
  #' @param log2_transform logical; czy zastosować log2-transformację wartości (domyślnie `FALSE`).
  #' @param pseudocount numeric; dodatek do wartości przed log2 (domyślnie `1`).
  #' @param group_cols character vector; kolumny z metadanych do grupowania (domyślnie `c("treatment","mouse_genotype")`).
  #' @param title character; tytuł wykresu (domyślnie generowany automatycznie).
  #' @param y_limits numeric vector długości 2; ograniczenia osi Y (opcjonalnie).
  #' @param levels character vector; kolejność wyświetlania grup (opcjonalnie).
  #' @param palette named character vector; własna paleta kolorów (np. `c("CTRL_WT"="#1f77b4", "DEX_WT"="#ff7f0e")`).
  #'
  #' @return Obiekt klasy `ggplot`.
  #' @examples
  #' plot_gene_boxplot_from_bundle(
  #'   bundle      = ris3q29_bundle_data_to_visualization,
  #'   gene_symbol = "Arc",
  #'   cluster     = "cluster_16",
  #'   data_type   = "qn",
  #'   log2_transform = TRUE,
  #'   levels      = c("CTRL_WT", "DEX_WT", "CORT_WT"),
  #'   palette     = c("CTRL_WT"="#1f77b4", "DEX_WT"="#ff7f0e", "CORT_WT"="#2ca02c")
  #' )
  stopifnot(is.list(bundle), "metadata" %in% names(bundle))
  data_type <- match.arg(data_type)
  
  # wybór macierzy
  mat_list_name <- switch(
    data_type,
    raw = "raw_counts_by_cluster",
    qn  = "qn_counts_by_cluster"
  )
  if (!(mat_list_name %in% names(bundle))) {
    stop(sprintf("W bundlu brakuje elementu '%s'.", mat_list_name))
  }
  if (!(cluster %in% names(bundle[[mat_list_name]]))) {
    stop(sprintf("Klaster '%s' nie istnieje w '%s'.", cluster, mat_list_name))
  }
  
  mat <- bundle[[mat_list_name]][[cluster]]
  if (is.null(dim(mat))) stop("Wybrany element nie jest macierzą.")
  
  # do długiego formatu i filtr po genie
  df_long <- mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name_id") %>%
    tidyr::separate("name_id", into = c("p1","p2","gene_symbol"), sep = "-", remove = FALSE) %>%
    dplyr::filter(gene_symbol == !!gene_symbol) %>%
    tidyr::pivot_longer(
      cols = -c(name_id, p1, p2, gene_symbol),
      names_to = "sample_ID",
      values_to = "value"
    ) %>%
    dplyr::left_join(bundle$metadata, by = "sample_ID")
  
  n_genes <- length(unique(df_long$name_id))
  if (n_genes == 0) {
    stop(sprintf("Nie znaleziono w klastrze '%s' wiersza dla gene_symbol == '%s'.", cluster, gene_symbol))
  } else if (n_genes > 1) {
    message(sprintf("⚠️ Uwaga: znaleziono %d rekordów dla gene_symbol '%s' w klastrze '%s'. Wszystkie zostaną pokazane (facet_wrap).",
                    n_genes, gene_symbol, cluster))
  }
  
  # log2 transform (opcjonalnie)
  if (isTRUE(log2_transform)) {
    df_long <- df_long %>%
      dplyr::mutate(value = log2(value + pseudocount))
    ylab_txt <- sprintf("log2(%s counts + %g)", data_type, pseudocount)
  } else {
    ylab_txt <- sprintf("%s counts", data_type)
  }
  
  # kolumna 'group'
  if (length(group_cols) > 0) {
    missing_cols <- setdiff(group_cols, colnames(df_long))
    if (length(missing_cols)) {
      stop(sprintf("Brakuje kolumn w metadata: %s", paste(missing_cols, collapse = ", ")))
    }
    df_long <- df_long %>%
      dplyr::mutate(group = interaction(dplyr::across(dplyr::all_of(group_cols)), sep = "_", drop = TRUE))
  } else {
    df_long <- df_long %>% dplyr::mutate(group = "all")
  }
  
  # ustawienie kolejności grup
  if (!is.null(levels)) {
    df_long$group <- factor(df_long$group, levels = levels)
  }
  
  if (is.null(title)) {
    title <- sprintf("%s | %s | %s", gene_symbol, cluster, toupper(data_type))
  }
  
  # wykres
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = title,
      x = "group",
      y = ylab_txt
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
  
  if (!is.null(y_limits) && length(y_limits) == 2) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }
  
  # własna paleta kolorów
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  
  # facet jeśli wiele rekordów
  if (n_genes > 1) {
    p <- p + ggplot2::facet_wrap(~ name_id, scales = "free_y")
  }
  
  return(p)http://localhost:9777/graphics/plot_zoom_png?width=815&height=900
}


plot_gene_boxplot_from_bundle(
  bundle      = ris3q29_bundle_data_to_visualization,
  gene_symbol = "Gfap",
  cluster     = "cluster_16",
  data_type   = "qn"
)


save(
  ris3q29_bundle_data_to_visualization,
  plot_gene_boxplot_from_bundle,
  file = "results/risperidone-3q29/Rdata/ris3q29_bundle_data_to_visualization_simpleVisMagda.RData"
)
