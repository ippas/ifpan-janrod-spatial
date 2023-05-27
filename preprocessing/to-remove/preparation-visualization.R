#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################
# example to run function
samples_vector1 <- meta_data %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]


samples_vector2 <- meta_data %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]


# visualization of UMAP
resolution = 0.5
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)

DimPlot(data_cluster, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(spatial_data = risperidone_st_data,
                resolution = 0.5,
                samples = c(samples_vector1, samples_vector2),
                palette = palette_allen, 
                size= 1.0, 
                ncol = 4)


spatial_interest_cluster(cluster = 12,
                         # seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = c(samples_vector1, samples_vector2),
                         size= 1,
                         ncol = 4)


# visualize features 
spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "seurat",
                  gene = "Polr3e",
                  samples =  c(samples_vector1, samples_vector2),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
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



