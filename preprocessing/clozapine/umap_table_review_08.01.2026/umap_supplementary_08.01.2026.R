DimPlot(
  clozapine_integrate_half,
  reduction = "seuart_clusters"
)


svg(filename = "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/umapClozapineRes0.45_08.01.2026.svg", 
    width = 10,
    height = 8)
Embeddings(
  object    = clozapine_integrate_half,
  reduction = "umap"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_barcode") %>% 
  left_join(., 
            clozapine_st_data_half$clusters %>% 
              mutate(sample_barcode = paste0(sample, "_", barcode)),
            by = "sample_barcode") %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster_resolution_0.45)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = palette_allen) +
  ggtitle("cluster_resolution_0.45")
dev.off()


svg(filename = "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/umapClozapineRes0.4_08.01.2026.svg", 
    width = 10,
    height = 8)
Embeddings(
  object    = clozapine_integrate_half,
  reduction = "umap"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_barcode") %>% 
  left_join(., 
            clozapine_st_data_half$clusters %>% 
              mutate(sample_barcode = paste0(sample, "_", barcode)),
            by = "sample_barcode") %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster_resolution_0.4)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = palette_allen) +
  ggtitle("cluster_resolution_0.4")

clozapine_st_data_half$clusters %>% 
  mutate(sample_barcode = paste0(sample, "_", barcode))

dev.off()
