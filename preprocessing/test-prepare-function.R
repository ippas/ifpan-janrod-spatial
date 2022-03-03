stat <- function(x, treat = colfilt.info$treatment, clusters = as.factor(colfilt.info$integrated_snn_res.2)) {
  idx = 1
  p <- vector()
  for(cluster in levels(clusters)) {
    p[idx] <- 1
    if (sum(clusters == cluster & treat == "ctrl") < 2 |
        sum(clusters == cluster & treat == "treat") < 2) {
      
    } else {
      p[idx] <- wilcox.test(x[clusters == cluster & treat == "ctrl"], x[clusters == cluster & treat == "treat"])$p.value
    }
    
    idx = idx + 1
  }
  names(p) <- levels(clusters)
  p
}


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

plot_feature <- function(data_cluster, peak_id, size){
  plots <- list()
  
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  for (i in 1:length(samples_name)) {
    
    plots[[i]] <- colfilt_norm_data %>% 
      .[peak_id,] %>%
      as.data.frame() %>% 
      rename(value = ".") %>% 
      rownames_to_column(var = "sample.barcode") %>% 
      separate("sample.barcode", c("sample", "barcode"), sep = "_") %>% 
      left_join(., data_cluster, by = c("barcode", "sample")) %>% 
      mutate(percentile0.99 = as.numeric(quantile(value, probs = c(0.99))),
             percentile0.05 = as.numeric(quantile(value, probs = c(0.05)))) %>%
      mutate(value = ifelse(value < percentile0.99, value, percentile0.99),
             value = ifelse(value > percentile0.05, value, percentile0.05)) %>%
      mutate(value = normalize(value)) %>%
      filter(sample == samples_name[i]) %>% 
      ggplot(aes(x=imagecol, y=imagerow,fill= value)) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colours = myPalette(100), limits=c(0,1))+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==samples_name[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==samples_name[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(samples_name[i])+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  }
  
  plot_grid(plotlist = plots)
  
}







plot_clusters <- function(data_cluster, size = 1.2){
  plots <- list()
  
  for (i in 1:length(samples_name)) {
    
    plots[[i]] <- bcs_merge %>% 
      left_join(., {data_cluster %>% 
          select(sample, barcode, cluster) %>%
          as.data.frame()}, by = c("sample", "barcode")) %>% na.omit %>%
      filter(sample ==samples_name[i]) %>%
      filter(tissue == "1") %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(cluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE) +
      scale_fill_manual(values = palette_cluster)+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==samples_name[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==samples_name[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(samples_name[i])+
      labs(fill = "Cluster")+
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
  
  plot_grid(plotlist = plots)
  
}




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

