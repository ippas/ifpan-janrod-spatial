#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

# visualization of UMAP
resolution = 1
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 3)

DimPlot(data_cluster, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                samples = samples_name,
                palette = palette_allen, 
                size= 1.2, 
                ncol = 3)

# visualize features 
spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = "Camk4",
                  samples = samples_name,
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 1,
                  ncol = 3,
                  normalization = T)


genetoplot = "Egr1"

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  tif_image = T,
                  normalization = T,
                  alpha = 1,
                  ncol = 3) 

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

tmp <- FindMarkers(data_cluster, ident.1 = 4, min.pct = 0.25)

tmp %>% rownames_to_column(., "peak_id") %>%
  left_join(., {mutate(info_peaks, peak_id = str_replace_all(peak_id, "_", "-")) %>%
      select(peak_id, gene_name)}, by = "peak_id") %>% head(20)


  
  info_peaks %>% 
  mutate(info_peaks, peak_id = str_replace_all(peak_id, "_", "-")) %>%
  select(peak_id, gene_name) %>% 
  filter(peak_id == "merged-samples-peak-199960")

  

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = "Cck",
                  samples = samples_name[c(1,5)],
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 1.6,
                  ncol = 2,
                  normalization = T)



