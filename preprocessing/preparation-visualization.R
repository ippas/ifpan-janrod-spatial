#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

# visualization of UMAP
resolution = 0.2
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)

DimPlot(data_cluster, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(seurat_object = integrated_analysis,
                spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                samples = samples_name,
                palette = palette_cluster, 
                size= 1)

# visualize features 

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_modification = "raw_data",
                  gene = "Junb",
                  samples = samples_name,
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 1,
                  filt_score_int = 0,
                  normalization = F) 


genetoplot = "Arc"

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_modification = "raw_data",
                  gene = genetoplot,
                  samples = samples_name,
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  filt_score_int = 0,
                  normalization = T) 



spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_modification = "range_normalize",
                  gene = genetoplot,
                  samples = samples_name,
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  filt_score_int = 0,
                  normalization = T)



spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_modification = "seurat",
                  gene = genetoplot,
                  samples = samples_name,
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  filt_score_int = 0,
                  normalization = T)

# visualize genes using seurat
spatial_gene_plot_seurat(data = spatial_transcriptomic_data$raw_data$annotate, 
                         gene = "",
                         filt_score_int = 0)



# old code

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

