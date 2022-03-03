geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

spatial_feature_plot <- function(spatial_data,
                                 type_modification,
                                 peak_id,
                                 samples,
                                 # barcode_data,
                                 # images_tibble,
                                 min_percentile = 0.01,
                                 max_percentile = 0.99,
                                 size = 1.2,
                                 normalization = TRUE){
  
  
  # prepare data to create plot
  preprocessing_data <- spatial_data[[type_modification]]$data[peak_id, ] %>%
    as.data.frame() %>%
    rename(value = ".") %>% 
    rownames_to_column("sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    # left_join(., data_cluster, by = c("barcode", "sample")) %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample")) %>%
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
          geom_spatial(data = images_tibble[spatial_data$images_information$sample == samples[i],],
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
          xlim(0, max(spatial_data$bcs_information %>%
                        filter(sample == samples[i]) %>%
                        select(width)))+
          ylim(max(spatial_data$bcs_information %>%
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
  
  rm(plot_list,
     preprocessing_data)
  
  return(wrap_plot)
  
}




