#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################
install.packages("roxygen2")
require(roxygen2)
install.packages("devtools")
require(devtools)

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
    spatial_feature_plot(object = integrated_analysis,
                         features = c(peak),
                         samples = samples_name) %>% 
      print()
  }
}


spatial_gene_plot_seurat(data = sample_anno_all, gene = "Mag")




spatial_feature_plot <- function(data,
                                 data_cluster, 
                                 peak_id,
                                 samples,
                                 barcode_data,
                                 images_tibble,
                                 min_percentile = 0.01,
                                 max_percentile = 0.99,
                                 size = 1.2,
                                 normalization = TRUE){
  
  # normalize <- function(x) {
  #   return((x - min(x)) / (max(x) - min(x)))
  # }
  
  # prepare data to create plot
  preprocessing_data <- data[peak_id, ] %>%
    as.data.frame() %>%
    rename(value = ".") %>% 
    rownames_to_column("sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    # left_join(., data_cluster, by = c("barcode", "sample")) %>%
    left_join(., bcs_merge, by = c("barcode", "sample")) %>%
    mutate(max_perc = as.numeric(quantile(value, probs = c(max_percentile))),
           min_perc = as.numeric(quantile(value, probs = c(min_percentile))),
           # trimming the extreme values to the adopted percentiles
           value = ifelse(value < max_perc, value, max_perc),
           value = ifelse(value > min_perc, value, min_perc),
           # zero one normalization if normalization = TRUE
           value = if (normalization == TRUE){
             (value - min(value)) / (max(value) - min(value))
           } else {
             value
           })

  plot_list <- list()

  # creating plot
  for (i in 1:length(samples)){
    sample_plot <- filter(preprocessing_data, sample == samples[i]) %>%
      # main part of creating plot
      {ggplot(., aes(x = imagecol, y = imagerow, fill = value)) +
          geom_spatial(data = images_tibble[images_tibble$sample == samples[i],],
                       aes(grob = grob),
                       x = 0.5, y = 0.5) +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0.25) +
          coord_cartesian(expand = FALSE) +
          # chose correct scale fill gradientn during condition 
          {if(normalization == TRUE){
            scale_fill_gradientn(colours = myPalette(100), limits = c(0,1))
          } else {
            scale_fill_gradientn(colours = myPalette(100))
          }} +
          xlim(0, max(barcode_data %>%
                        filter(sample == samples[i]) %>%
                        select(width)))+
          ylim(max(barcode_data %>%
                     filter(sample == samples[i]) %>%
                     select(height)), 0) +
          # aesthetics plot
          xlab("") +
          ylab("") +
          ggtitle(samples[i])+
          theme_set(theme_bw(base_size = 10))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank())}


    plot_list[[samples[i]]] <- sample_plot

  }
 
  wrap_plot <- wrap_plots(plot_list) +
    plot_layout(ncol = 3, guides = "collect")

  rm(plot_list)

  return(wrap_plot)

}



spatial_feature_plot(data = colfilt_norm_data,
                     data_cluster = data_cluster,
                     samples = samples_name,
                     peak_id = "merged-samples-peak-173070",
                     barcode_data = bcs_merge,
                     images_tibble = images_tibble,
                     max_percentile = 0.99,
                     min_percentile = 0.05,
                     normalization = TRUE)
 












 




