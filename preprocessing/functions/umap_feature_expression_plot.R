umap_feature_expression_plot <- function(ldopa_integrate,
                                         spatial_data,
                                         type_data,
                                         peak_id,
                                         samples,
                                         plot_title = "UMAP Visualization with Feature Expression",  # Nowy argument dla tytułu
                                         min_percentile = 0.01,
                                         max_percentile = 0.99,
                                         normalization = TRUE,
                                         low_color = "gray99",
                                         high_color = "red",
                                         na_color = "grey99",
                                         point_size = 1,
                                         alpha = 1,
                                         save_to_png = FALSE,        # Nowy argument: czy zapisać do PNG
                                         file_path = NULL,           # Nowy argument: ścieżka do zapisu
                                         width = 8,                  # Opcjonalny argument: szerokość wykresu
                                         height = 6,                 # Opcjonalny argument: wysokość wykresu
                                         dpi = 300) {                # Opcjonalny argument: rozdzielczość
  
  # Wyciągnij koordynaty UMAP oraz barcodes
  umap_coordinates <- ldopa_integrate@reductions$umap@cell.embeddings
  umap_table <- data.frame(
    barcode = rownames(umap_coordinates),
    X = umap_coordinates[, 1], # Oś X
    Y = umap_coordinates[, 2]  # Oś Y
  )
  
  # Przygotowanie tabeli ekspresji
  expression_table <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    dplyr::rename(expression = ".") %>%
    rownames_to_column("sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample")) %>%
    filter(sample %in% samples) %>%
    mutate(
      max_perc = as.numeric(quantile(expression, probs = c(max_percentile), na.rm = TRUE)),
      min_perc = as.numeric(quantile(expression, probs = c(min_percentile), na.rm = TRUE)),
      expression = ifelse(expression > max_perc, max_perc, expression),
      expression = ifelse(expression < min_perc, min_perc, expression),
      expression = if (normalization == TRUE) {
        (expression - min(expression, na.rm = TRUE)) / (max(expression, na.rm = TRUE) - min(expression, na.rm = TRUE))
      } else {
        expression
      },
      barcode = paste(sample, barcode, sep = "_")
    ) %>%
    dplyr::select(barcode, expression)
  
  # Merge dwóch tabel
  merged_table <- dplyr::left_join(umap_table, expression_table, by = "barcode")
  
  # Sprawdzenie poprawności mergu
  all_matched <- nrow(merged_table) == nrow(umap_table)
  missing_matches <- sum(is.na(merged_table$expression))
  if (!all_matched || missing_matches > 0) {
    cat("Uwaga: Niektóre wiersze nie zostały dopasowane. Niedopasowane wiersze:", missing_matches, "\n")
  }
  
  # Tworzenie wykresu
  plot <- ggplot(merged_table, aes(x = X, y = Y, fill = expression)) +
    geom_point(shape = 21, size = point_size, alpha = alpha, stroke = 0) +  # Punkty bez obramówki
    scale_fill_gradient(low = low_color, high = high_color, na.value = na_color) +  # Skala kolorów
    labs(
      title = plot_title,  # Użycie nowego argumentu dla tytułu
      x = "UMAP 1",
      y = "UMAP 2",
      fill = "Expression"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
  
  # Opcjonalny zapis do pliku
  if (save_to_png) {
    # Ustawienie domyślnej ścieżki, jeśli nie została podana
    if (is.null(file_path)) {
      file_path <- getwd()
    }
    
    # Upewnienie się, że ścieżka istnieje
    if (!dir.exists(file_path)) {
      dir.create(file_path, recursive = TRUE)
      cat("Utworzono katalog:", file_path, "\n")
    }
    
    # Przygotowanie bezpiecznej nazwy pliku
    safe_title <- stringr::str_replace_all(plot_title, "[: ]", "_")
    safe_title <- stringr::str_replace_all(safe_title, "_+", "_")
    file_name <- paste0(safe_title, ".png")
    full_path <- file.path(file_path, file_name)
    
    # Dostosowanie motywu do białego tła
    plot_to_save <- plot + theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
    
    # Zapis wykresu
    ggsave(filename = full_path, plot = plot_to_save, width = width, height = height, dpi = dpi)
    cat("Wykres zapisano do pliku:", full_path, "\n")
  }
  
  return(plot)
}


# 
# # Wywołanie funkcji z Twoimi danymi
# plot <- umap_feature_expression_plot(
#   ldopa_integrate = ldopa_integrate,
#   spatial_data = ldopa_st_data,
#   type_data = "quantile_normalize_resolution_1",
#   peak_id = "ldopa-peak-17758",
#   samples = c(samples_saline, samples_ldopa),
#   min_percentile = 0.00,
#   max_percentile = 1,
#   normalization = TRUE,
#   low_color = "gray99",
#   high_color = "red",
#   na_color = "grey99",
#   point_size = 1,
#   alpha = 1
# )
# 
# # Wyświetlenie wykresu
# print(plot)


umap_gene_expression_plot <- function(spatial_data,
                                      type_data,
                                      ncol,
                                      gene,
                                      ldopa_integrate,
                                      samples,
                                      min_percentile = 0.01,
                                      max_percentile = 0.99,
                                      normalization = TRUE,
                                      low_color = "gray95",
                                      high_color = "red",
                                      na_color = "grey99",
                                      point_size = 1,
                                      alpha = 1,
                                      save_to_png = FALSE,        # Nowy argument: czy zapisać do PNG
                                      file_path = NULL,           # Nowy argument: ścieżka do zapisu
                                      width = 8,                  # Opcjonalny argument: szerokość wykresu
                                      height = 6,                 # Opcjonalny argument: wysokość wykresu
                                      dpi = 300,                  # Opcjonalny argument: rozdzielczość
                                      ...) {
  
  # Wyciągnięcie wektora peak_id dla danego genu
  vector_peak <- spatial_data[[type_data]]$annotate %>%
    dplyr::filter(gene_name == gene) %>%
    dplyr::select(peak_id) %>% 
    pull(peak_id)
  
  print(paste("Peaki dla genu", gene, ":", paste(vector_peak, collapse = ", ")))
  
  # Lista do przechowywania wykresów
  plot_list <- list()
  
  # Iteracja po każdym peak_id
  for (peak in vector_peak) {
    # Tworzenie tytułu wykresu
    plot_title <- paste(gene, peak, sep = ": ")
    
    # Generowanie wykresu UMAP
    plot <- umap_feature_expression_plot(
      ldopa_integrate = ldopa_integrate,
      spatial_data = spatial_data,
      type_data = type_data,
      peak_id = peak,
      samples = samples,
      plot_title = plot_title,
      min_percentile = min_percentile,
      max_percentile = max_percentile,
      normalization = normalization,
      low_color = low_color,
      high_color = high_color,
      na_color = na_color,
      point_size = point_size,
      alpha = alpha,
      save_to_png = save_to_png,    # Przekazanie argumentu
      file_path = file_path,        # Przekazanie argumentu
      width = width,                # Przekazanie argumentu
      height = height,              # Przekazanie argumentu
      dpi = dpi                    # Przekazanie argumentu
    )
    
    # Dodanie wykresu do listy
    plot_list[[peak]] <- plot
    
    # Opcjonalne wyświetlenie wykresu
    print(plot_title)
    print(plot)
  }
  
  # Organizacja wykresów w siatkę
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol, guides = "collect") +
    patchwork::plot_annotation(title = paste("UMAP Plots dla Genu:", gene))
  
  # Wyświetlenie połączonego wykresu
  print(combined_plot)
  
  # Opcjonalny zapis połączonego wykresu
  if (save_to_png) {
    # Ustawienie domyślnej ścieżki, jeśli nie została podana
    if (is.null(file_path)) {
      file_path <- getwd()
    }
    
    # Upewnienie się, że ścieżka istnieje
    if (!dir.exists(file_path)) {
      dir.create(file_path, recursive = TRUE)
      cat("Utworzono katalog:", file_path, "\n")
    }
    
    # Przygotowanie bezpiecznej nazwy pliku
    safe_title <- gsub("[/:*?\"<>|\\s]", "_", paste("Combined_", gene, sep = ""))
    file_name <- paste0(safe_title, ".png")
    full_path <- file.path(file_path, file_name)
    
    # # Zapis połączonego wykresu
    # ggsave(filename = full_path, plot = combined_plot, width = width, height = height, dpi = dpi)
    # cat("Połączony wykres zapisano do pliku:", full_path, "\n")
  }
  
  # Zwrócenie listy wykresów oraz połączonego wykresu
  return(list(individual_plots = plot_list, combined_plot = combined_plot))
}



umap_gene_expression_plot(
  spatial_data = risperidone_st_data_half,
  type_data = "quantile_normalize_resolution_0.4",
  ncol = 2,
  gene = "Gad2",
  ldopa_integrate = risperidone_integrate_half,
  samples = c(samples_saline, samples_risperidone),
  min_percentile = 0.01,
  max_percentile = 0.99,
  normalization = TRUE,
  low_color = "gray95",
  high_color = "red",
  na_color = "grey99",
  point_size = 1.5,
  alpha = 1,
  save_to_png = TRUE,                        # Włączenie zapisu do PNG
  file_path = "./results/risperidone/umap-gene-expression",                # Podanie ścieżki do katalogu
  width = 10,                                # Szerokość wykresu w calach
  height = 8,                                # Wysokość wykresu w calach
  dpi = 300                                  # Rozdzielczość
)


# Definiowanie wektora z nazwami genów
genes <- c("Drd1", "Drd2", "Adora2a", "Ppp1r1b", "Ppp1r2", "Gad2", "Ecel1", "Gfra1")




