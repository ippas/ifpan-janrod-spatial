
DimPlot(risperidone_integrate_half, reduction = "umap",
        split.by = "sample", ncol = 4)

DimPlot(ris3q29_integrate, reduction = "umap")


spatial_gene_plot(spatial_data = ris3q29_st_data,
                  type_data = "raw_data",
                  gene = "Sgk1",
                  samples =  c(samples_salWt, samples_risWt),
                  # samples = sample_names,
                  min_percentile = 0.00,
                  max_percentile = 0.99,
                  size = 0.7,
                  ncol = 7,
                  tif_image = T,
                  normalization = F)

spatial_cluster_select(
  spatial_data = ris3q29_st_data,
  resolution = 0.4,
  samples = c(samples_salWt, samples_risWt),
  palette = palette_allen,
  size = 1,
  select_clusters = c(0),
  tif_image = T,
  ncol = 4
)


spatial_feature_aggregated_cluster_plot(
  spatial_data = ris3q29_st_data,
  summary_statistics = risWtSalWt_summary_statistics,
  peak_id = "peak-35029-Sgk1",
  type_data = "raw_data",
  samples = c(samples_salWt, samples_risWt),
  cluster_resolution = "cluster_resolution_0.4",
  summary_metric = "mean",
  normalization = TRUE,
  tif_image = TRUE,
  show_legend = TRUE,
  return_list = FALSE
)

spatial_feature_aggregated_cluster_plot(
  spatial_data = ris3q29_st_data,
  summary_statistics = risWtSalWt_summary_statistics,
  type_data = "quantile_normalize_resolution_0.4",
  samples = c(samples_salWt, samples_risWt),
  peak_id = "peak-35029-Sgk1",
  cluster_resolution = "cluster_resolution_0.4",        # lub "cluster_resolution_0.4" – obie zadziałają
  summary_metric = "mean",
  clusters = c(0, 3, 4),                     # tylko te klastry będą wyświetlane
  normalization = TRUE,
  tif_image = TRUE,
  show_legend = TRUE
)

peak_id <- "peak-35029-Sgk1"
ris3q29_st_data$raw_data$data[peak_id, ] %>% str()

# spatial_cluster()

sample_names

tmp <- FindClusters(ris3q29_integrate, resolution = 0.25)

DimPlot(tmp, reduction = "umap")

DimPlot(risperidone_integrate_half, reduction = "umap", label = TRUE)


risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_16$control$mean %>% head


#
ris3q29_integrate@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "barcode") 


svg("results/risperidone-3q29/figures/umap/umap-res0.4-26samples-risWtDelsalWtDel.svg", 
    width = 13, height = 12)

ris3q29_st_data$clusters %>% 
  dplyr::select(c(sample, barcode, cluster_resolution_0.4)) %>% 
  mutate(barcode = paste0(sample, "_", barcode)) %>% 
  dplyr::select(c(barcode, cluster_resolution_0.4)) %>% 
  left_join(
    {
      ris3q29_integrate@reductions$umap@cell.embeddings %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "barcode")
    },
    by = "barcode"
  ) %>% 
  mutate(cluster = as.factor(cluster_resolution_0.4)) %>% 
  {
    p <- ggplot(., aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
      geom_point(size = 0.01) +
      scale_color_manual(
        values = palette_allen[levels(.$cluster)],
        labels = as.character(as.integer(levels(.$cluster)) + 1)
      ) +
      guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_classic() +
      theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14)
      )
    print(p)
  }

dev.off()

  
palette_allen
  
