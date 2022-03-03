#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################



##### testing code to refresh how visualize data
resolution = 0.2
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)


data_cluster <- FindClusters(integrated_analysis, resolution = resolution) %>% 
  .$seurat_clusters %>%
  as.data.frame() %>% 
  rename(cluster = ".") %>% 
  rownames_to_column(var = "sample_barcode") %>%
  separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  left_join(., bcs_merge, by = c("barcode", "sample"))

plot_clusters(data_cluster)

colfilt_anno %>%
  filter(gene_name == "Fos") %>%
  select(gene_name, peak_id)


plot_feature(data_cluster = data_cluster,
             peak_id = "merged-samples-peak-94135",
             size = 1)





spatial_feature_plot(spatial_data = spatial_transcriptomic_data,
                     type_modification = "range_normalize",
                     peak_id = "merged-samples-peak-173070",
                     samples = samples_name,
                     min_percentile = 0.05,
                     max_percentile = 0.99,
                     size = 1,
                     normalization = TRUE)

