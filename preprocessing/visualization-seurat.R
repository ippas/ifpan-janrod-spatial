#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################


integrated_analysis <- FindClusters(integrated_analysis, resolution = 0.5)


# adding images to integrated_analysis
for(sample_name in samples_name) {
  dir_image <- paste("data/spaceranger-corrected/",
                     sample_name,
                     "/outs/spatial/",
                     sep = "")
   
  sample_image <- Read10X_Image(dir_image)
  
  integrated_analysis@images[[sample_name]] <- sample_image
  
  rm(sample_image,
     sample_name)
}

# changing coordinates inside integrated_analysis
for(sample_name in samples_name) {
  new_rownames <- integrated_analysis@images[[sample_name]]@coordinates %>%
    rownames() %>%
    paste(sample_name, ., sep = "_") 
  
  rownames(integrated_analysis@images[[sample_name]]@coordinates) <- new_rownames
  
  rm(sample_name,
     new_rownames)
}


#####################################################
# prepare function to visualize using Seurat object #
#####################################################

spatial_feature_plot_seurat <- function(..., samples){
  plot_list <- list()
  # function create plot for all slide inside seurat object
  
  for(sample_index in 1:length(samples)){
    
    sample_plot <-SpatialFeaturePlot(...)[[sample_index]] +
      xlim(0,557)
    
    plot_list[[samples[sample_index]]] <- sample_plot
    
    rm(sample_plot)
  }
  
  wrap_plot <- wrap_plots(plot_list) + 
    plot_layout(ncol = 2, guides = "collect")
  
  rm(plot_list)
  
  return(wrap_plot)
}



spatial_feature_plot_seurat(object = integrated_analysis,
                            features = c("merged-samples-peak-173070"),
                            samples = samples_name)



spatial_gene_plot_seurat <- function(data, gene, filt_score_int = 1000){
  # function create plots for all peaks describe to gene
  # gene name gene for which create plots
  # filt_score_int threshold to filter peaks using metrics from MACS3
  
  # find peaks for interest gene, and using filter threshold for peaks
  data %>%
    filter(gene_name == gene) %>%
    filter(score_int..10.log10pvalue. > filt_score_int) %>%
    select(peak_id) %>% .[,1] -> vector_peak
  
  # create plot for find peaks
  for(peak in vector_peak){
    spatial_feature_plot_seurat(object = integrated_analysis,
                         features = c(peak),
                         samples = samples_name) +
      plot_annotation(title = paste(gene, peak, sep = ": ")) -> peak_plot 
      print(peak_plot)
      
   rm(peak_plot)   
  }
}






spatial_transcriptomic_data$raw_data$annotate %>%
filter(gene_name == "Mag") %>%
  filter(score_int..10.log10pvalue. > 1000) %>%
  select(peak_id) %>% .[,1]












 




