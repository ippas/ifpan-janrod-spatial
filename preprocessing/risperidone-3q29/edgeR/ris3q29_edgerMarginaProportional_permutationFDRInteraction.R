
# ##############################################################################
# ---- functions ----
# ##############################################################################
# ============================================================
# 🧮 Funkcja wysokiego poziomu — łączenie wyników edgeR (2 czynniki)
# ============================================================

combine_edger_results_full <- function(res_edger_object,
                                       label = "MarginalProportional",
                                       verbose = TRUE) {
  start_time <- Sys.time()
  if (verbose) message("🚀 Rozpoczynam scalanie wyników edgeR: ", label)
  
  # Klucze wspólne
  keys <- c("cluster", "name_id", "peak_id", "gene_symbol")
  
  # Pomocnicza funkcja lokalna
  combine_with_clusters <- function(results, 
                                    slot = c("anova", "posthoc"), 
                                    element = NULL,
                                    filter_expr = NULL) {
    slot <- match.arg(slot)
    
    tmp <- lapply(results, function(x){
      obj <- x[[slot]]
      if (!is.null(element)) {
        if (!element %in% names(obj)) {
          stop(paste("Element", element, "not found in slot", slot))
        }
        obj <- obj[[element]]
      }
      if (!is.null(filter_expr)) {
        obj <- dplyr::filter(obj, !!rlang::parse_expr(filter_expr))
      }
      obj
    })
    
    combined_data <- dplyr::bind_rows(tmp, .id = "cluster")
    tmp$all_cluster <- combined_data
    tmp
  }
  
  # ============================================================
  # 📊 Tworzenie tabel cząstkowych
  # ============================================================
  results <- res_edger_object$results
  
  t_global <- combine_with_clusters(results, slot = "anova", element = "global")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC.mouse_genotypewtdel_global",
                       "logFC.treatmentrisperidone_global",
                       "logFC.mouse_genotypewtdel.treatmentrisperidone_global",
                       "F_global","P_global","FDR_global"))
  
  t_interaction <- combine_with_clusters(results, slot = "anova", element = "interaction")$all_cluster %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol","logCPM",
                       "F_interaction","P_interaction","FDR_interaction"))
  
  t_main_first <- combine_with_clusters(results, slot = "anova", element = "main_first")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_allWtAllDel","F_allWtAllDel","P_allWtAllDel","FDR_allWtAllDel"))
  
  t_main_second <- combine_with_clusters(results, slot = "anova", element = "main_second")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_allSalAllRis","F_allSalAllRis","P_allSalAllRis","FDR_allSalAllRis"))
  
  t_posthoc_salRisWt <- combine_with_clusters(results, slot = "posthoc",
                                              element = "SECOND_risperidone-saline|FIRST_wtwt")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_salRisWt","F_salRisWt","P_salRisWt","FDR_salRisWt"))
  
  t_posthoc_salRisDel <- combine_with_clusters(results, slot = "posthoc",
                                               element = "SECOND_risperidone-saline|FIRST_wtdel")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_salRisDel","F_salRisDel","P_salRisDel","FDR_salRisDel"))
  
  t_posthoc_risWtDel <- combine_with_clusters(results, slot = "posthoc",
                                              element = "FIRST_wtdel-wtwt|SECOND_risperidone")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_risWtDel","F_risWtDel","P_risWtDel","FDR_risWtDel"))
  
  t_posthoc_salWtDel <- combine_with_clusters(results, slot = "posthoc",
                                              element = "FIRST_wtdel-wtwt|SECOND_saline")$all_cluster %>%
    dplyr::select(-logCPM) %>%
    rlang::set_names(c("cluster","name_id","peak_id","gene_symbol",
                       "logFC_salWtDel","F_salWtDel","P_salWtDel","FDR_salWtDel"))
  
  # ============================================================
  # 🔗 Łączenie wszystkiego
  # ============================================================
  full_results_df <- purrr::reduce(
    list(
      t_global,
      t_interaction,
      t_main_first,
      t_main_second,
      t_posthoc_salRisWt,
      t_posthoc_salRisDel,
      t_posthoc_risWtDel,
      t_posthoc_salWtDel
    ),
    dplyr::full_join, by = keys
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(cluster, name_id, peak_id, gene_symbol)
  
  end_time <- Sys.time()
  if (verbose) {
    message("✅ Scalanie zakończone: ", label)
    message("⏱️ Czas wykonania: ", round(end_time - start_time, 2), " ", attr(end_time - start_time, "units"))
  }
  
  return(full_results_df)
}


run_permutation_edger_FDR <- function(
  raw_data_by_resolution,
  sample_meta,
  n_permutations = 100,
  seed = 777,
  memory_limit_gb = 150
) {
  stopifnot(is.data.frame(sample_meta))
  stopifnot(n_permutations > 0)
  
  total_start <- Sys.time()
  set.seed(seed)
  seeds <- sample(1e6:9e6, n_permutations)
  
  times <- numeric(n_permutations)
  all_results <- vector("list", n_permutations)
  
  message("🚀 Start analizy FDR z ", n_permutations, " permutacjami...\n")
  
  for (i in seq_len(n_permutations)) {
    iter_start <- Sys.time()
    
    # 🔹 Pomiar zużycia RAM
    mem_gb <- as.numeric(pryr::mem_used()) / (1024^3)
    if (mem_gb > memory_limit_gb) {
      stop(
        sprintf("❌ Przekroczono limit pamięci (%.1f GB > %.0f GB)", 
                mem_gb, memory_limit_gb)
      )
    }
    
    # 🔹 Permutacja metadanych
    set.seed(seeds[i])
    perm_meta <- sample_meta %>%
      dplyr::mutate(sample_ID = sample(sample_ID))
    
    # 🚀 Uruchom analizę edgeR
    res <- run_edger_twofactor_for_clusters_v2(
      raw_data_by_resolution = raw_data_by_resolution,
      sample_meta            = perm_meta,
      sample_col             = "sample_ID",
      first_factor_col       = "mouse_genotype",
      second_factor_col      = "treatment",
      first_levels           = c("wtwt", "wtdel"),
      second_levels          = c("saline", "risperidone"),
      main_effects           = "marginal",
      main_effects_weights   = "proportional",
      clusters               = NULL
    )
    
    # 🔗 Połącz wyniki
    full_res <- combine_edger_results_full(
      res_edger_object = res,
      label = paste0("perm_", i)
    )
    
    # 📊 Filtrowanie i agregacja wyników po klastrze
    perm_summary <- full_res %>%
      filter(P_interaction < 0.01) %>%
      filter(
        (abs(logFC_salRisDel) > 0.8 & abs(logFC_salRisWt) < 0.2) |
          (abs(logFC_salRisDel) < 0.2 & abs(logFC_salRisWt) > 0.8)
      ) %>%
      split(.$cluster) %>%
      map_int(nrow) %>%
      enframe(name = "cluster", value = "n_results") %>%
      right_join(
        tibble(cluster = paste0("cluster_", 0:19)),
        by = "cluster"
      ) %>%
      mutate(n_results = replace_na(n_results, 0)) %>%
      mutate(cluster = factor(cluster, levels = paste0("cluster_", 0:19))) %>%
      arrange(cluster) %>%
      bind_rows(
        tibble(
          cluster = "total_results",
          n_results = sum(.$n_results, na.rm = TRUE)
        )
      )
    
    all_results[[i]] <- perm_summary
    
    # ⏱️ Pomiar czasu iteracji
    iter_end <- Sys.time()
    times[i] <- as.numeric(iter_end - iter_start, units = "mins")
    
    # 📊 Postęp i ETA (w minutach)
    avg_time <- mean(times[1:i])
    eta <- avg_time * (n_permutations - i)
    message(sprintf(
      "🔁 Iteracja %d/%d ukończona (%.2f min, ETA ~ %.1f min, RAM: %.1f GB)",
      i, n_permutations, times[i], eta, mem_gb
    ))
  }
  
  # 🧩 Scalanie wszystkich wyników permutacji
  merged_results <- all_results %>%
    set_names(paste0("perm_", seq_along(.))) %>%
    imap(~ rename(.x, !!paste0("nResults_", .y) := n_results)) %>%
    reduce(full_join, by = "cluster") %>%
    arrange(cluster) %>%
    mutate(
      mean_nResults = rowMeans(select(., starts_with("nResults_")), na.rm = TRUE),
      sd_nResults = apply(select(., starts_with("nResults_")), 1, sd, na.rm = TRUE)
    )
  
  # 🕒 Podsumowanie czasu
  total_end <- Sys.time()
  total_elapsed <- as.numeric(total_end - total_start, units = "mins")
  time_summary <- tibble(
    total_time_min = total_elapsed,
    mean_time_min = mean(times),
    sd_time_min = sd(times)
  )
  
  message("\n✅ Wszystkie permutacje zakończone.")
  message("⏱️ Całkowity czas: ", round(time_summary$total_time_min, 1), " min")
  message("📊 Średni czas: ", round(time_summary$mean_time_min, 2), " ± ", round(time_summary$sd_time_min, 2), " min")
  
  return(list(
    merged_results = merged_results,
    permutation_results = all_results,
    time_summary = time_summary
  ))
}



# ##############################################################################
# ---- run analysis ----
# ##############################################################################
resultsFDR_perm100_interactionEdger <- run_permutation_edger_FDR(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta = metadata_ris3q29 %>%  filter(sample_ID != "S13839Nr3"),
  n_permutations = 100,
  seed = 777,
  memory_limit_gb = 100
)

resultsFDR_perm100_interactionEdger$merged_results %>% 
  select(-c(cluster, mean_nResults, sd_nResults)) 

c(cluster_0 = 3, 
  cluster_1 = 2,
  cluster_2 = 1,
  cluster_3 = 2,
  cluster_4 = 3,
  cluster_5 = 5,
  cluster_6 = 1,
  cluster_7 = 3,
  cluster_8 = 4,
  cluster_9 = 2,
  cluster_10 = 1,
  cluster_11 = 1,
  cluster_12 = 1,
  cluster_13 = 1,
  cluster_14 = 7,
  cluster_15 = 1,
  cluster_16 = 9,
  cluster_17 = 5,
  cluster_18 = 10,
  cluster_19 = 2,
  total_results = 64
)


resultsFDR_perm100_interactionEdger$merged_results %>%
  pivot_longer(
    cols = starts_with("nResults_perm_"),
    names_to = "perm_id",
    values_to = "nResults"
  ) %>%
  left_join(
    tibble(
      cluster = names(c(
        cluster_0 = 3, cluster_1 = 2, cluster_2 = 1, cluster_3 = 2,
        cluster_4 = 3, cluster_5 = 5, cluster_6 = 1, cluster_7 = 3,
        cluster_8 = 4, cluster_9 = 2, cluster_10 = 1, cluster_11 = 1,
        cluster_12 = 1, cluster_13 = 1, cluster_14 = 7, cluster_15 = 1,
        cluster_16 = 9, cluster_17 = 5, cluster_18 = 10, cluster_19 = 2,
        total_results = 64
      )),
      nResults_original = as.numeric(c(
        cluster_0 = 3, cluster_1 = 2, cluster_2 = 1, cluster_3 = 2,
        cluster_4 = 3, cluster_5 = 5, cluster_6 = 1, cluster_7 = 3,
        cluster_8 = 4, cluster_9 = 2, cluster_10 = 1, cluster_11 = 1,
        cluster_12 = 1, cluster_13 = 1, cluster_14 = 7, cluster_15 = 1,
        cluster_16 = 9, cluster_17 = 5, cluster_18 = 10, cluster_19 = 2,
        total_results = 64
      ))
    ),
    by = "cluster"
  ) %>%
  group_by(cluster) %>%
  summarise(
    n_perm = n(),
    n_Results_original = unique(nResults_original),
    n_betterRandom = sum(nResults > nResults_original[1]),
    n_betterOrEqualRandom = sum(nResults >= nResults_original[1]),
    permFDR_betterRanodm = (n_betterRandom) / (n_perm),
    permFDR_betterOrEqualRanodm = (n_betterOrEqualRandom) / (n_perm),
    .groups = "drop"
  ) %>%
  arrange(factor(cluster, levels = c(paste0("cluster_", 0:19), "total_results"))) %>% 
  as.data.frame()


save(resFDR_perm100_interactionEdger,
     file = "results/risperidone-3q29/Rdata/permutationFDR_edgerInteraction/resFDR_perm100_interactionEdger_23.10.2025.RData")

# ##############################################################################
# ---- run 1000 permutation ----
# ##############################################################################
resultsFDR_perm1000_interactionEdger <- run_permutation_edger_FDR(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta = metadata_ris3q29 %>%  filter(sample_ID != "S13839Nr3"),
  n_permutations = 1000,
  seed = 777,
  memory_limit_gb = 100
)

save(resFDR_perm1000_interactionEdger,
     file = "results/risperidone-3q29/Rdata/permutationFDR_edgerInteraction/resFDR_perm1000_interactionEdger_23.10.2025.RData")


resultsFDR_perm100_interactionEdger$merged_results %>% 
  select(-c(cluster, mean_nResults, sd_nResults)) 

c(cluster_0 = 3, 
  cluster_1 = 2,
  cluster_2 = 1,
  cluster_3 = 2,
  cluster_4 = 3,
  cluster_5 = 5,
  cluster_6 = 1,
  cluster_7 = 3,
  cluster_8 = 4,
  cluster_9 = 2,
  cluster_10 = 1,
  cluster_11 = 1,
  cluster_12 = 1,
  cluster_13 = 1,
  cluster_14 = 7,
  cluster_15 = 1,
  cluster_16 = 9,
  cluster_17 = 5,
  cluster_18 = 10,
  cluster_19 = 2,
  total_results = 64
)


resultsFDR_perm1000_interactionEdger$merged_results %>%
  pivot_longer(
    cols = starts_with("nResults_perm_"),
    names_to = "perm_id",
    values_to = "nResults"
  ) %>%
  left_join(
    tibble(
      cluster = names(c(
        cluster_0 = 3, cluster_1 = 2, cluster_2 = 1, cluster_3 = 2,
        cluster_4 = 3, cluster_5 = 5, cluster_6 = 1, cluster_7 = 3,
        cluster_8 = 4, cluster_9 = 2, cluster_10 = 1, cluster_11 = 1,
        cluster_12 = 1, cluster_13 = 1, cluster_14 = 7, cluster_15 = 1,
        cluster_16 = 9, cluster_17 = 5, cluster_18 = 10, cluster_19 = 2,
        total_results = 64
      )),
      nResults_original = as.numeric(c(
        cluster_0 = 3, cluster_1 = 2, cluster_2 = 1, cluster_3 = 2,
        cluster_4 = 3, cluster_5 = 5, cluster_6 = 1, cluster_7 = 3,
        cluster_8 = 4, cluster_9 = 2, cluster_10 = 1, cluster_11 = 1,
        cluster_12 = 1, cluster_13 = 1, cluster_14 = 7, cluster_15 = 1,
        cluster_16 = 9, cluster_17 = 5, cluster_18 = 10, cluster_19 = 2,
        total_results = 64
      ))
    ),
    by = "cluster"
  ) %>%
  group_by(cluster) %>%
  summarise(
    n_perm = n(),
    n_Results_original = unique(nResults_original),
    n_betterRandom = sum(nResults > nResults_original[1]),
    n_betterOrEqualRandom = sum(nResults >= nResults_original[1]),
    permFDR_betterRanodm = (n_betterRandom) / (n_perm),
    permFDR_betterOrEqualRanodm = (n_betterOrEqualRandom) / (n_perm),
    .groups = "drop"
  ) %>%
  arrange(factor(cluster, levels = c(paste0("cluster_", 0:19), "total_results"))) %>% 
  as.data.frame()

resultsFDR_perm1000_interactionEdger$merged_results$mean_nResults
resultsFDR_perm1000_interactionEdger$merged_results %>% 
  filter(cluster == "total_results") %>% 
  select(-c(sd_nResults,mean_nResults, cluster)) %>% 
  tail(1) 


resultsFDR_perm1000_interactionEdger$merged_results %>%
  filter(cluster == "total_results") %>%
  select(nResults) %>%                  # wybieramy kolumnę z liczbą wyników
  arrange(nResults) %>%                 # sortujemy rosnąco
  mutate(rank = row_number()) %>%       # dodajemy rangę
  mutate(FDR_perm = rank / n()) %>%     # wyliczamy procentowo pozycję (np. 43/1000)
  mutate(is_target = nResults == 64) %>%# zaznaczamy gdzie jest 64
  filter(is_target | row_number() %in% c(1, n())) 


resultsFDR_perm1000_interactionEdger$merged_results %>% 
  filter(cluster == "cluster_19") %>% 
  select(-c(sd_nResults,mean_nResults, cluster)) %>% 
  tail(1) %>% as.numeric() %>% as.data.frame() %>% 
  set_colnames("nResults") %>% 
  arrange(desc(nResults)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(FDR_perm = rank / n()) %>% 
  mutate(is_target = nResults == 2) %>% 
  filter(is_target | row_number() %in% c(1, n())) 
  

