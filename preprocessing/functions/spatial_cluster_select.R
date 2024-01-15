spatial_cluster_select <- function(spatial_data,
                            resolution,
                            # seurat_object,
                            samples,
                            palette, 
                            select_clusters = NULL,
                            size = 1.2,
                            alpha = 1,
                            tif_image = TRUE,
                            ncol){
  
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  # prepare data for plots
  data_cluster <- spatial_data$clusters %>% 
    dplyr::select(sample, barcode, name_column_resolution) %>%
    dplyr::rename(cluster = name_column_resolution) %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
  if (is.null(select_clusters)) {
    data_cluster
  } else {
    data_cluster %>% 
      filter(cluster %in% select_clusters) -> data_cluster
  }
  
  plot_list <- list()
  
  # prepare list of plots
  for (sample_id in samples){
    sample_plot <- filter(data_cluster, sample == sample_id) %>%
      {ggplot(., aes(x = imagecol, y = imagerow, fill = factor(cluster))) +
          {if(tif_image == TRUE){
            geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
                         aes(grob = grob),
                         x = 0.5, y = 0.5) 
          } else {}} +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0, alpha = alpha) +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0)+
          coord_cartesian(expand = FALSE) +
          scale_fill_manual(values = palette)+
          xlim(0, max(spatial_data$bcs_information %>%
                        filter(sample == sample_id) %>%
                        dplyr::select(width)))+
          ylim(max(spatial_data$bcs_information %>%
                     filter(sample == sample_id) %>%
                     dplyr::select(height)), 0) +
          # aesthetics plot
          xlab("") +
          ylab("") +
          ggtitle(
            paste(sample_id, 
                  ": ",
                  spatial_data$sample_information[spatial_data$sample_information$sample_ID == sample_id,]$treatment,
                  ", ",
                  spatial_data$sample_information[spatial_data$sample_information$sample_ID == sample_id,]$mouse_genotype,
                  sep = "")
          ) +
          guides(fill = guide_legend(override.aes = list(size=3)))+
          theme_set(theme_bw(base_size = 10))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank()) +
          NoLegend()
      }
    
    plot_list[[sample_id]] <- sample_plot
    
  }
  
  wrap_plot <- wrap_plots(plot_list) +
    plot_layout(ncol = ncol, guides = "collect")
  
  rm(plot_list)
  
  return(wrap_plot)
  
}


spatial_cluster_select(spatial_data = risperidone_st_data_half,
                resolution = 0.8,
                samples = c(samples_saline, samples_risperidone),
                palette = palette_allen,
                size= 1.0,
                tif_image = F,
                # select_clusters = c(1,2),
                ncol = 4)
