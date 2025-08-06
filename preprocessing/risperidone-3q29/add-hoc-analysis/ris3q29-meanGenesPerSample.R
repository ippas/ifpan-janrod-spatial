ris3q29_st_data$raw_data$data[1:10, 1:10]


mat <- ris3q29_st_data$raw_data$data
genes <- sub('.*-(.*)$', '\\1', rownames(mat))

# Dla każdej kolumny
gene_counts_per_spot <- apply(mat, 2, function(col) {
  # bierzemy tylko niezerowe wiersze
  present_genes <- genes[col != 0]
  length(unique(present_genes))
})

# Z wektora numerycznego do tibble
df <- enframe(gene_counts_per_spot, name = "sample_barcode", value = "n_genes")

# Rozdziel kolumnę "sample_barcode" na "sample" i "barcode"
df <- df %>%
  separate(sample_barcode, into = c("sample", "barcode"), sep = "_")

df$n_genes %>% median

# Podsumowanie grupami
df %>%
  group_by(sample) %>%
  summarise(
    mean_genes = mean(n_genes, na.rm = TRUE),
    sd_genes = sd(n_genes, na.rm = TRUE),
    n_spots = n()
  ) %>% as.data.frame() %>% 
  filter(sample != "S13839Nr3") -> summary_df

write.table(
  summary_df,
  file = "results/risperidone-3q29/ad-hoc-results/ris3q29-meanGenesPerSample.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_df,
  file = "results/risperidone-3q29/ad-hoc-results/ris3q29-meanGenesPerSample.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)


summary_df$mean_genes %>% mean
