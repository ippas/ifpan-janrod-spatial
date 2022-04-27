#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

# visualization of UMAP
resolution = 1
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)

DimPlot(data_cluster, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                samples = samples_name,
                palette = palette_cluster, 
                size= 1, 
                ncol = 4)

# visualize features 
spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = "Junb",
                  samples = samples_name,
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 1,
                  ncol = 4,
                  normalization = F)


genetoplot = "Egr1"

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  normalization = T,
                  ncol = 4) 

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "range_normalize",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  normalization = T,
                  ncol = 3)


spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "seurat",
                  gene = genetoplot,
                  samples = samples_name,
                  min_percentile = 0.0,
                  max_percentile = 1,
                  normalization = T,
                  ncol = 4)

spatial_gene_plot_cluster(spatial_data = spatial_transcriptomic_data,
                          type_data = "raw_data",
                          gene = "Junb",
                          samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                          min_percentile = 0.0,
                          max_percentile = 1,
                          size = 1,
                          normalization = T,
                          resolution = 1,
                          clusters = 0,
                          ncol = 3) 



spatial_interest_cluster(cluster = 0,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                         size= 1, 
                         ncol = 3)






# old code
data_cluster <- FindClusters(integrated_analysis, resolution = resolution) %>% 
  .$seurat_clusters %>%
  as.data.frame() %>% 
  rename(cluster = ".") %>% 
  rownames_to_column(var = "sample_barcode") %>%
  separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  left_join(., bcs_merge, by = c("barcode", "sample"))
