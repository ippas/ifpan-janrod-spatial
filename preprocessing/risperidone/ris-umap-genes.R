spatial_gene_plot(spatial_data = risperidone_st_data_half,
                  type_data = "raw_data",
                  gene = "Gad2",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)

# Wywołanie funkcji z Twoimi danymi
plot <- umap_feature_expression_plot(
  ldopa_integrate = risperidone_integrate_half,
  spatial_data = risperidone_st_data_half,
  type_data = "quantile_normalize_resolution_0.4",
  peak_id = "risperidone-peak-125061",
  samples = c(samples_saline, samples_risperidone),
  min_percentile = 0.01,
  max_percentile = 0.99,
  normalization = TRUE,
  low_color = "gray95",
  high_color = "red",
  na_color = "grey99",
  point_size = 1,
  alpha = 1,
  plot_title = paste("Gad2:", "risperidone-peak-125061")
)

# Wyświetlenie wykresu
print(plot)




# Definiowanie wektora z nazwami genów
genes <- c("Drd1", "Drd2", "Adora2a", "Ppp1r1b", "Ppp1r2", "Gad2", "Ecel1", "Gfra1")

# Pętla iterująca po każdym genie w wektorze
for (gene in genes) {
  # Wywołanie funkcji umap_gene_expression_plot dla aktualnego genu
  umap_gene_expression_plot(
    spatial_data = risperidone_st_data_half,
    type_data = "quantile_normalize_resolution_0.4",
    ncol = 2,
    gene = gene,
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
  
  # Opcjonalnie: Informacja o postępie w konsoli
  message("Wykres dla genu ", gene, " został zapisany.")
}


