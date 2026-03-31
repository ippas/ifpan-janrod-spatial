
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
  data,
  metadata,
  n_permutations = 100,
  seed = 777
) {
  # 📋 Kontrola wejść
  stopifnot(is.data.frame(metadata))
  stopifnot(n_permutations > 0)
  
  # 🕒 Start pomiaru czasu całkowitego
  total_start <- Sys.time()
  
  # 🔢 Ustaw wektor seedów
  set.seed(seed)
  seeds <- sample(1e6:9e6, n_permutations)
  
  # 📊 Ramka do zapisu wyników
  all_results <- vector("list", n_permutations)
  times <- numeric(n_permutations)
  
  message("🚀 Start analizy FDR z ", n_permutations, " permutacjami...\n")
  
  for (i in seq_len(n_permutations)) {
    iter_start <- Sys.time()
    
    # 🔹 Permutuj metadane
    perm_meta <- metadata %>%
      dplyr::mutate(sample_ID = sample(sample_ID))
    
    # 🔹 Uruchom edgeR
    res <- run_edger_twofactor_for_clusters_v2(
      raw_data_by_resolution = data,
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
    
    # 📉 Filtrowanie wyników
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
    times[i] <- as.numeric(iter_end - iter_start, units = "secs")
    
    # 📊 Postęp i ETA
    avg_time <- mean(times[1:i])
    eta <- avg_time * (n_permutations - i)
    message(sprintf(
      "🔁 Iteracja %d/%d ukończona (%.1fs, ETA ~ %.1f s)",
      i, n_permutations, times[i], eta
    ))
  }
  
  # 🕒 Czas całkowity
  total_end <- Sys.time()
  total_elapsed <- total_end - total_start
  
  # 📈 Podsumowanie czasu
  time_summary <- tibble(
    total_time_min = as.numeric(total_elapsed, units = "mins"),
    mean_time_s = mean(times),
    sd_time_s = sd(times)
  )
  
  message("\n✅ Wszystkie permutacje zakończone.")
  message("⏱️ Całkowity czas: ", round(time_summary$total_time_min, 1), " min")
  message("📊 Średni czas: ", round(time_summary$mean_time_s, 2), " s ± ", round(time_summary$sd_time_s, 2))
  
  return(list(
    permutation_results = all_results,
    time_summary = time_summary
  ))
}



# ##############################################################################
# ---- run analysis ----
# ##############################################################################
# 🕒 Start pomiaru czasu
start_time <- Sys.time()

# 🚀 Uruchom analizę edgeR (efekty marginalne, wagi proporcjonalne)
res_edger_MarginalProportional_23.10.2025 <- run_edger_twofactor_for_clusters_v2(
  raw_data_by_resolution = wtDel_summary_statistics$raw_data$resolution_0.4,
  sample_meta            =   filter(sample_ID != "S13839Nr3") %>% 
    dplyr::mutate(sample_ID = sample(sample_ID)),
  sample_col             = "sample_ID",
  first_factor_col       = "mouse_genotype",
  second_factor_col      = "treatment",
  first_levels           = c("wtwt", "wtdel"),
  second_levels          = c("saline", "risperidone"),
  main_effects           = "marginal",          # ← efekty główne (marginal)
  main_effects_weights   = "proportional",      # ← wagi proporcjonalne dla niebalansowanych grup
  clusters               = NULL                 # ← NULL = przetwarza wszystkie klastry
)

# 🕒 Koniec pomiaru czasu
end_time <- Sys.time()

# ⏱️ Oblicz czas trwania
elapsed_time <- end_time - start_time

# 📊 Podsumowanie
message("✅ Analiza zakończona pomyślnie.")
message("⏱️ Czas wykonania: ", round(elapsed_time, 2), " ", attr(elapsed_time, "units"))


# 🔗 Łączenie wszystkich wyników edgeR (two-factor, marginal-proportional, ris3q29)
edger_twoFactorMarginalProportional_fullResults_ris3q29 <-  combine_edger_results_full(
  res_edger_object = res_edger_MarginalProportional_23.10.2025,
  label = "ris3q29_MarginalProportional"
)

edger_twoFactorMarginalProportional_fullResults_ris3q29 %>% 
  filter(P_interaction < 0.01) %>%
  filter((abs(logFC_salRisDel) > 0.8 & abs(logFC_salRisWt) < 0.2) | (abs(logFC_salRisDel) < 0.2 & abs(logFC_salRisWt) > 0.8)) %>% nrow

edger_twoFactorMarginalProportional_fullResults_ris3q29 %>%
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
  ) %>% as.data.frame()


