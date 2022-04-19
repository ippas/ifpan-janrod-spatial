###############################
# function to display barcode #
###############################

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

###############################################
# prepare function to visualize interest peak #
###############################################

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
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == samples[i],],
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
          ggtitle(
            paste(samples_name[i], 
                  ": ",
                  spatial_data$sample_information$treatment[i], 
                  ", ",
                  spatial_data$sample_information$mouse_genotype[i],
                  sep = "")
          ) +
          theme_set(theme_bw(base_size = 10))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank())}
    
    
    plot_list[[samples[i]]] <- sample_plot
    
  }
  
  wrap_plot <- wrap_plots(plot_list)
  
  rm(plot_list,
     preprocessing_data)
  
  return(wrap_plot)
  
}



#########################################
# prepare function to visualize cluster #
#########################################

spatial_cluster <- function(seurat_object,
                            spatial_data,
                            resolution,
                            samples,
                            palette, 
                            size= 1.2,
                            ncol){
  
  # prepare data for plots
  data_cluster <- FindClusters(seurat_object, resolution = resolution) %>% 
    .$seurat_clusters %>%
    as.data.frame() %>% 
    rename(cluster = ".") %>% 
    rownames_to_column(var = "sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
  plot_list <- list()

  # prepare list of plots
  for (i in 1:length(samples)){
    sample_plot <- filter(data_cluster, sample == samples[i]) %>%
      {ggplot(., aes(x = imagecol, y = imagerow, fill = factor(cluster))) +
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == samples[i],],
                       aes(grob = grob),
                       x = 0.5, y = 0.5) +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
          coord_cartesian(expand = FALSE) +
          scale_fill_manual(values = palette)+
          xlim(0, max(spatial_data$bcs_information %>%
                        filter(sample == samples[i]) %>%
                        select(width)))+
          ylim(max(spatial_data$bcs_information %>%
                     filter(sample == samples[i]) %>%
                     select(height)), 0) +
          # aesthetics plot
          xlab("") +
          ylab("") +
          ggtitle(
            paste(samples_name[i], 
                  ": ",
                  spatial_data$sample_information$treatment[i], 
                  ", ",
                  spatial_data$sample_information$mouse_genotype[i],
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
    
    plot_list[[samples[i]]] <- sample_plot
    
  }
  
  wrap_plot <- wrap_plots(plot_list) +
    plot_layout(ncol = ncol, guides = "collect")
  
  rm(plot_list)
  
  return(wrap_plot)
  
}



spatial_cluster(seurat_object = integrated_analysis,
                spatial_data = spatial_transcriptomic_data,
                resolution = 0.2,
                samples = samples_name,
                palette = palette_cluster, 
                size = 1,
                ncol = 4)


############################################################
# prepare function for visualize interest position cluster #
############################################################
spatial_interest_cluster <- function(cluster,
                                     ...){
  mypalette <- rep("#555555", 38)
  names(mypalette) <- as.character(0:37)
  mypalette[cluster] <- "#00FF00"
  
  spatial_cluster(palette = mypalette, 
                  ...)
  
}

spatial_interest_cluster(cluster = 5,
                         seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = samples_name,
                         size= 1,
                         ncol = 4)

########################################################
# prepare function to visualize peaks discribe to gene #
########################################################

spatial_gene_plot <- function(spatial_data,
                              type_modification,
                              ncol,
                              ...,
                              gene){
  
  spatial_data[[type_modification]]$annotate %>%
    filter(gene_name == gene) %>%
    select(peak_id) %>% .[,1] -> vector_peak
  
  for(peak in vector_peak){
    spatial_feature_plot(spatial_data = spatial_data,
                         type_modification = type_modification,
                         peak_id = peak,
                         ...) + 
      plot_layout(ncol = ncol, guides = "collect") + 
      plot_annotation(title = paste(gene, peak, sep = ": ")) -> peak_plot
    print(peak_plot)  
  }
  
  rm(vector_peak,
     peak_plot)
}



# visualize features 
spatial_feature_plot(spatial_data = spatial_transcriptomic_data,
                     type_modification = "range_normalize",
                     peak_id = "merged-samples-peak-19995",
                     samples = samples_name,
                     min_percentile = 0,
                     max_percentile = 1,
                     size = 1,
                     normalization = TRUE)




palette_cluster <- c("0" = "#b2df8a",
                     "1" = "#e41a1c",
                     "2" = "#377eb8",
                     "3" ="#4daf4a",
                     "4" = "#ff7f00",
                     "5" = "gold", 
                     "6" = "#a65628", 
                     "7" = "#999999", 
                     "8" = "black", 
                     "9" = "grey", 
                     "10" = "white", 
                     "11" = "purple",
                     "12" = "red", 
                     "13" = "blue", 
                     "14" = "pink",
                     "15" = "brown",
                     "16" = "green",
                     "17" = "tomato1",
                     "18" = "yellow3",
                     "19" = "violet",
                     "20" = "yellowgreen",
                     "21" = "lightblue1",
                     "22" = "lightblue4",
                     "23" = "lightgoldenrod3",
                     "24" = "lightpink2",
                     "25" = "magenta",
                     "26" = "limegreen",
                     "27" = "maroon",
                     "28" = "mintcream",
                     "29"=  "oldlace",
                     "30" = "tan",
                     "31" = "chartreuse",
                     "32" = "blue4",
                     "33" = "midnightblue",
                     "34" = "slategray4",
                     "35" = "snow3",
                     "36" = "springgreen",
                     "37" = "plum1")



