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
# Custom palette function: gray to red
my_custom_palette <- function(n) {
  colors <- colorRampPalette(c("#EEEEEE", "#ff0000", "#bf0000", "#800000", "#400000"))(n)
  return(colors)
}

spatial_feature_plot <- function(spatial_data,
                                 type_data,
                                 peak_id,
                                 samples,
                                 min_percentile = 0.01,
                                 max_percentile = 0.99,
                                 size = 1.2,
                                 alpha = 1,
                                 tif_image = TRUE,
                                 normalization = TRUE,
                                 show_legend = TRUE,
                                 return_list = FALSE){
  
  # prepare data to create plot
  preprocessing_data <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    dplyr::rename(value = ".") %>% 
    rownames_to_column("sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample")) %>%
    mutate(max_perc = as.numeric(quantile(value, probs = c(max_percentile))),
           min_perc = as.numeric(quantile(value, probs = c(min_percentile))),
           value = ifelse(value < max_perc, value, max_perc),
           value = ifelse(value > min_perc, value, min_perc),
           value = if (normalization == TRUE){
             (value - min(value)) / (max(value) - min(value))
           } else {
             value
           })
  
  plot_list <- list()
  
  # creating plot
  for (sample_id in samples){
    sample_plot <- filter(preprocessing_data, sample == sample_id) %>%
      ggplot(aes(x = imagecol, y = imagerow, fill = value)) +
      {if(tif_image == TRUE){
        geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
                     aes(grob = grob),
                     x = 0.5, y = 0.5) 
      } else {}} +
      geom_point(shape = 21, colour = "black", size = size, stroke = 0, alpha = alpha) +
      coord_cartesian(expand = FALSE) +
      {if(normalization == TRUE){
        scale_fill_gradientn(colours = my_custom_palette(100), limits = c(0,1))
      } else {
        scale_fill_gradientn(colours = my_custom_palette(100))
      }} +
      xlim(0, max(spatial_data$bcs_information %>%
                    filter(sample == sample_id) %>%
                    dplyr::select(width)))+
      ylim(max(spatial_data$bcs_information %>%
                 filter(sample == sample_id) %>%
                 dplyr::select(height)), 0) +
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
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    if (!show_legend) {
      sample_plot <- sample_plot + theme(legend.position = "none")
    }
    
    plot_list[[sample_id]] <- sample_plot
  }
  
  if(return_list == TRUE){
    return(plot_list)
  } else {
    wrap_plot <- wrap_plots(plot_list)
    return(wrap_plot)
  }
  
  rm(preprocessing_data)
}




#########################################
# prepare function to visualize cluster #
#########################################

spatial_cluster <- function(spatial_data,
                            resolution,
                            samples,
                            palette, 
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
    dplyr::left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
  # data_cluster <- FindClusters(seurat_object, resolution = resolution) %>% 
  #   .$seurat_clusters %>%
  #   as.data.frame() %>% 
  #   dplyr::rename(cluster = ".") %>% 
  #   rownames_to_column(var = "sample_barcode") %>%
  #   separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  #   left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
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


# # example to run function
# spatial_cluster(spatial_data = spatial_transcriptomic_data,
#                 # seurat_object = integrated_analysis,
#                 resolution = 1,
#                 samples = samples_name,
#                 palette = palette_cluster, 
#                 size = 1,
#                 ncol = 4)


############################################################
# prepare function for visualize interest position cluster #
############################################################
spatial_interest_cluster <- function(cluster,
                                     ...){
  mypalette <- rep("#555555", 38)
  names(mypalette) <- as.character(0:37)
  # mypalette[cluster + 1] <- "#00FF00"
  mypalette[cluster + 1] <- palette_allen[cluster + 1]
  
  spatial_cluster(palette = mypalette, 
                  ...)
  
}

# # example to run function
# spatial_interest_cluster(cluster = 2,
#                          # seurat_object = integrated_analysis,
#                          spatial_data = spatial_transcriptomic_data,
#                          resolution = 1,
#                          samples = samples_name,
#                          size= 1,
#                          ncol = 4)

########################################################
# prepare function to visualize peaks describe to gene #
########################################################

spatial_gene_plot <- function(spatial_data,
                              type_data,
                              ncol,
                              ...,
                              gene){
  
  spatial_data[[type_data]]$annotate %>%
    filter(gene_name == gene) %>%
    dplyr::select(peak_id) %>% .[,1] -> vector_peak
  
  print(vector_peak)
  
  for(peak in vector_peak){
    spatial_feature_plot(spatial_data = spatial_data,
                         type_data = type_data,
                         peak_id = peak,
                         ...) + 
      plot_layout(ncol = ncol, guides = "collect") + 
      plot_annotation(title = paste(gene, peak, sep = ": ")) -> peak_plot
    print(peak)
    print(peak_plot)
    # return(peak_plot)  
  }
  
  rm(vector_peak,
     peak_plot)
}



# # example run function
# spatial_feature_plot(spatial_data = spatial_transcriptomic_data,
#                      type_data = "range_normalize",
#                      peak_id = "merged-samples-peak-19995",
#                      samples = samples_name,
#                      min_percentile = 0,
#                      max_percentile = 1,
#                      size = 1,
#                      normalization = TRUE)




############################################################
# prepare function to visualize interest peak and clusters #
############################################################
spatial_feature_plot_cluster <- function(spatial_data,
                                         type_data,
                                         peak_id,
                                         samples,
                                         min_percentile = 0.01,
                                         max_percentile = 0.99,
                                         size = 1.2,
                                         resolution = 1,
                                         show_legend = TRUE,
                                         tif_image = TRUE,
                                         clusters,
                                         normalization = TRUE){
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  
  data_cluster <- spatial_data$clusters %>% 
    dplyr::select(sample, barcode, name_column_resolution) %>%
    dplyr::rename(cluster = name_column_resolution)
  
  # prepare data to create plot
  preprocessing_data <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    dplyr::rename(value = ".") %>% 
    rownames_to_column("sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
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
           }) %>%
    left_join(., data_cluster, by = c("barcode", "sample")) %>%
    filter(cluster %in% clusters)
  
  plot_list <- list()
  
  # creating plot
  for (sample_id in samples){
    sample_plot <- filter(preprocessing_data, sample == sample_id) %>%
      # main part of creating plot
      {ggplot(., aes(x = imagecol, y = imagerow, fill = value)) +
          {if(tif_image == TRUE){
            geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
                         aes(grob = grob),
                         x = 0.5, y = 0.5) 
          } else {}} +
          # geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
          #              aes(grob = grob),
          #              x = 0.5, y = 0.5) +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0) +
          coord_cartesian(expand = FALSE) +
          # chose correct scale fill gradientn during condition 
          # {if(normalization == TRUE){
          #   scale_fill_gradientn(colours = myPalette(100), limits = c(0,1))
          # } else {
          #   scale_fill_gradientn(colours = myPalette(100))
          # }} +
          {if(normalization == TRUE){
            scale_fill_gradientn(colours = my_custom_palette(100), limits = c(0,1))
          } else {
            scale_fill_gradientn(colours = my_custom_palette(100))
          }} +
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
          theme_set(theme_bw(base_size = 10))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank())}
    
    if (!show_legend) {
      sample_plot <- sample_plot + theme(legend.position = "none")
    }
    
    plot_list[[sample_id]] <- sample_plot
    
  }
  
  wrap_plot <- wrap_plots(plot_list)
  
  rm(plot_list,
     preprocessing_data)
  
  return(wrap_plot)
  
}


############################################################
# prepare function to visualize interest gene and clusters #
############################################################
spatial_gene_plot_cluster <- function(spatial_data,
                                      type_data,
                                      ncol,
                                      ...,
                                      gene){
  
  spatial_data[[type_data]]$annotate %>%
    filter(gene_name == gene) %>%
    dplyr::select(peak_id) %>% .[,1] -> vector_peak
  
  for(peak in vector_peak){
    spatial_feature_plot_cluster(spatial_data = spatial_data,
                                 type_data = type_data,
                                 peak_id = peak,
                                 ...) + 
      plot_layout(ncol = ncol, guides = "collect") + 
      plot_annotation(title = paste(gene, peak, sep = ": ")) -> peak_plot
    print(peak_plot)  
  }
  
  rm(vector_peak,
     peak_plot)
}


# spatial_gene_plot_cluster(spatial_data = spatial_transcriptomic_data,
#                   type_data = "raw_data",
#                   gene = genetoplot,
#                   samples = samples_name,
#                   min_percentile = 0.0,
#                   max_percentile = 1,
#                   size = 1,
#                   normalization = T,
#                   resolution = 1,
#                   clusters = c(4),
#                   ncol = 4) 


spatial_feature_aggregated_cluster_plot <- function(spatial_data,
                                                    summary_statistics,
                                                    type_data = "raw_data",
                                                    peak_id,
                                                    samples,
                                                    cluster_resolution = "cluster_resolution_0.4",
                                                    summary_metric = "mean",
                                                    normalization = TRUE,
                                                    min_percentile = 0.01,
                                                    max_percentile = 0.99,
                                                    size = 1.2,
                                                    alpha = 1,
                                                    tif_image = TRUE,
                                                    show_legend = TRUE,
                                                    return_list = FALSE) {
  # --- Etap 1: Wstępne przetwarzanie danych agregowanych ---
  cluster_resotution_summary_statistics <- sub("^cluster_", "", cluster_resolution)
  cluster_names <- names(summary_statistics[[type_data]][[cluster_resotution_summary_statistics]])
  
  preprocessing_summary_statistics <- cluster_names %>%
    purrr::set_names() %>%
    purrr::map_dfr(function(clust_name) {
      control_df <- summary_statistics[[type_data]][[cluster_resotution_summary_statistics]][[clust_name]]$control[[summary_metric]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "peak") %>%
        tidyr::pivot_longer(-peak, names_to = "sample", values_to = "summary_value") %>%
        dplyr::mutate(cluster = clust_name)
      
      experiment_df <- summary_statistics[[type_data]][[cluster_resotution_summary_statistics]][[clust_name]]$experiment[[summary_metric]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "peak") %>%
        tidyr::pivot_longer(-peak, names_to = "sample", values_to = "summary_value") %>%
        dplyr::mutate(cluster = clust_name)
      
      dplyr::bind_rows(control_df, experiment_df)
    }, .id = NULL) %>%
    dplyr::filter(peak == peak_id) %>%
    dplyr::filter(sample %in% samples) %>%
    dplyr::mutate(summary_value = if (normalization) {
      (summary_value - min(summary_value, na.rm = TRUE)) / (max(summary_value, na.rm = TRUE) - min(summary_value, na.rm = TRUE))
    } else {
      summary_value
    })
  
  # --- Etap 2: Przygotowanie danych przestrzennych ---
  preprocessing_data <- tibble::enframe(spatial_data[[type_data]]$data[peak_id, ], name = "sample_barcode", value = "value") %>%
    tidyr::separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    dplyr::left_join(spatial_data$bcs_information, by = c("barcode", "sample")) %>%
    dplyr::left_join(spatial_data$clusters %>%
                       dplyr::select(sample, barcode, !!rlang::sym(cluster_resolution)) %>%
                       dplyr::rename(cluster = !!rlang::sym(cluster_resolution)),
                     by = c("barcode", "sample")) %>%
    dplyr::mutate(
      max_perc = as.numeric(stats::quantile(value, probs = max_percentile)),
      min_perc = as.numeric(stats::quantile(value, probs = min_percentile)),
      value = ifelse(value < max_perc, value, max_perc),
      value = ifelse(value > min_perc, value, min_perc),
      value = if (normalization) {
        (value - min(value)) / (max(value) - min(value))
      } else {
        value
      },
      cluster = paste0("cluster_", cluster)
    )
  
  # --- Etap 3: Dołączenie danych agregowanych ---
  preprocessing_data <- dplyr::left_join(
    preprocessing_data,
    preprocessing_summary_statistics %>% dplyr::select(-peak),
    by = c("sample", "cluster")
  ) %>%
    dplyr::mutate(value = summary_value) %>%
    dplyr::select(-summary_value)
  
  # --- Etap 4: Generowanie wykresów ---
  plot_list <- list()
  
  for (sample_id in samples) {
    sample_plot <- preprocessing_data %>%
      dplyr::filter(sample == sample_id) %>%
      ggplot2::ggplot(ggplot2::aes(x = imagecol, y = imagerow, fill = value)) +
      {
        if (tif_image) {
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id, ],
                       ggplot2::aes(grob = grob),
                       x = 0.5, y = 0.5)
        }
      } +
      ggplot2::geom_point(shape = 21, colour = "black", size = size, stroke = 0, alpha = alpha) +
      ggplot2::coord_cartesian(expand = FALSE) +
      {
        if (normalization) {
          ggplot2::scale_fill_gradientn(colours = my_custom_palette(100), limits = c(0, 1))
        } else {
          ggplot2::scale_fill_gradientn(colours = my_custom_palette(100))
        }
      } +
      ggplot2::xlim(0, max(dplyr::filter(spatial_data$bcs_information, sample == sample_id)$width)) +
      ggplot2::ylim(max(dplyr::filter(spatial_data$bcs_information, sample == sample_id)$height), 0) +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::ggtitle(paste0(sample_id, ": ",
                              spatial_data$sample_information[
                                spatial_data$sample_information$sample_ID == sample_id, ]$treatment,
                              ", ",
                              spatial_data$sample_information[
                                spatial_data$sample_information$sample_ID == sample_id, ]$mouse_genotype)) +
      ggplot2::theme_set(ggplot2::theme_bw(base_size = 10)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
    
    if (!show_legend) {
      sample_plot <- sample_plot + ggplot2::theme(legend.position = "none")
    }
    
    plot_list[[sample_id]] <- sample_plot
  }
  
  # --- Etap 5: Zwracanie wyników ---
  if (return_list) {
    return(plot_list)
  } else {
    return(patchwork::wrap_plots(plot_list))
  }
}



spatial_gene_aggregated_cluster_plot <- function(spatial_data,
                                                 summary_statistics,
                                                 type_data,
                                                 cluster_resolution = "cluster_resolution_0.4",
                                                 summary_metric = "mean",
                                                 normalization = TRUE,
                                                 tif_image = TRUE,
                                                 show_legend = TRUE,
                                                 return_list = FALSE,
                                                 samples,
                                                 gene,
                                                 ncol = 4,
                                                 size = 1.2,
                                                 alpha = 1) {
  
  # Wyszukiwanie peak_id przypisanych do genu
  vector_peak <- spatial_data[[type_data]]$annotate %>%
    dplyr::filter(gene_name == gene) %>%
    dplyr::pull(peak_id)
  
  print(vector_peak)
  
  for (peak in vector_peak) {
    peak_plot <- spatial_feature_aggregated_cluster_plot(
      spatial_data = spatial_data,
      summary_statistics = summary_statistics,
      type_data = type_data,
      peak_id = peak,
      samples = samples,
      cluster_resolution = cluster_resolution,
      summary_metric = summary_metric,
      normalization = normalization,
      tif_image = tif_image,
      show_legend = show_legend,
      return_list = return_list,
      size = size,
      alpha = alpha
    ) +
      patchwork::plot_layout(ncol = ncol, guides = "collect") +
      patchwork::plot_annotation(title = paste(gene, peak, sep = ": "))
    
    print(peak)
    print(peak_plot)
  }
  
  # Czyszczenie tymczasowych zmiennych
  rm(vector_peak, peak_plot)
}



spatial_feature_plot_cluster_summary_metric <- function(spatial_data,
                                                        type_data,
                                                        peak_id,
                                                        samples,
                                                        min_percentile = 0.01,
                                                        max_percentile = 0.99,
                                                        size = 1.2,
                                                        resolution = 1,
                                                        show_legend = TRUE,
                                                        tif_image = TRUE,
                                                        clusters,
                                                        normalization = TRUE,
                                                        summary_metric = "mean") {
  
  name_column_resolution <- paste("cluster_resolution", resolution, sep = "_")
  
  data_cluster <- spatial_data$clusters %>% 
    dplyr::select(sample, barcode, !!name_column_resolution) %>%
    dplyr::rename(cluster = !!name_column_resolution)
  
  preprocessing_data <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    dplyr::rename(value = ".") %>% 
    rownames_to_column("sample_barcode") %>%
    tidyr::separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample")) %>%
    dplyr::mutate(
      max_perc = as.numeric(quantile(value, probs = max_percentile)),
      min_perc = as.numeric(quantile(value, probs = min_percentile)),
      value = ifelse(value < max_perc, value, max_perc),
      value = ifelse(value > min_perc, value, min_perc),
      value = if (normalization) {
        (value - min(value)) / (max(value) - min(value))
      } else {
        value
      }
    ) %>%
    left_join(., data_cluster, by = c("barcode", "sample")) %>%
    dplyr::mutate(value = ifelse(cluster %in% clusters, value, 0))
  
  plot_list <- list()
  
  for (sample_id in samples) {
    sample_plot <- dplyr::filter(preprocessing_data, sample == sample_id) %>%
      {
        ggplot(., aes(x = imagecol, y = imagerow, fill = value)) +
          {if (tif_image) {
            geom_spatial(
              data = spatial_data$images_information[spatial_data$images_information$sample == sample_id, ],
              aes(grob = grob),
              x = 0.5, y = 0.5
            )
          }} +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0) +
          coord_cartesian(expand = FALSE) +
          {if (normalization) {
            scale_fill_gradientn(colours = my_custom_palette(100), limits = c(0, 1))
          } else {
            scale_fill_gradientn(colours = my_custom_palette(100))
          }} +
          xlim(0, max(spatial_data$bcs_information %>%
                        dplyr::filter(sample == sample_id) %>%
                        dplyr::pull(width))) +
          ylim(max(spatial_data$bcs_information %>%
                     dplyr::filter(sample == sample_id) %>%
                     dplyr::pull(height)), 0) +
          labs(
            x = NULL, y = NULL,
            title = paste0(
              sample_id, ": ",
              spatial_data$sample_information[
                spatial_data$sample_information$sample_ID == sample_id,
              ]$treatment,
              ", ",
              spatial_data$sample_information[
                spatial_data$sample_information$sample_ID == sample_id,
              ]$mouse_genotype
            )
          ) +
          theme_bw(base_size = 10) +
          theme(
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank()
          )
      }
    
    if (!show_legend) {
      sample_plot <- sample_plot + theme(legend.position = "none")
    }
    
    plot_list[[sample_id]] <- sample_plot
  }
  
  wrap_plot <- patchwork::wrap_plots(plot_list)
  
  rm(plot_list, preprocessing_data)
  
  return(wrap_plot)
}


barplot_feature_expression_clusters <- function(
  summary_data,
  data_type = "raw_data",
  resolution = "resolution_0.4",
  metric = "mean",
  selected_peak,
  test = "wilcoxon_test",
  same_y_scale = TRUE,
  clusters = NULL,     # klastry do pominięcia
  ncol = 5             # liczba kolumn w wrap_plots()
) {
  test_names <- c(
    t_test = "t-test",
    wilcoxon_test = "Wilcoxon test",
    ks_test = "Kolmogorow–Smirnow test"
  )
  
  cluster_ids <- names(summary_data[[data_type]][[resolution]])
  
  if (!is.null(clusters)) {
    excluded <- paste0("cluster_", clusters)
    cluster_ids <- setdiff(cluster_ids, excluded)
  }
  
  if (same_y_scale) {
    max_y_val <- max(sapply(cluster_ids, function(cluster_id) {
      cluster_data <- summary_data[[data_type]][[resolution]][[cluster_id]]
      control_vals <- cluster_data$control[[metric]][selected_peak, ]
      experiment_vals <- cluster_data$experiment[[metric]][selected_peak, ]
      max(c(control_vals, experiment_vals), na.rm = TRUE)
    }))
  }
  
  plot_list <- list()
  pvalue_list <- list()
  
  for (cluster_id in cluster_ids) {
    
    cluster_data <- summary_data[[data_type]][[resolution]][[cluster_id]]
    
    control_mat <- cluster_data$control[[metric]]
    experiment_mat <- cluster_data$experiment[[metric]]
    
    control_long <- control_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("peak") %>%
      pivot_longer(cols = -peak, names_to = "sample", values_to = "value") %>%
      mutate(group = "control")
    
    experiment_long <- experiment_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("peak") %>%
      pivot_longer(cols = -peak, names_to = "sample", values_to = "value") %>%
      mutate(group = "experiment")
    
    long_df <- bind_rows(control_long, experiment_long) %>%
      filter(peak == selected_peak)
    
    if (nrow(long_df) == 0) next
    
    test_df <- tibble::tibble(
      peak = cluster_data$peak,
      gene = cluster_data$gene
    ) %>%
      bind_cols(as.data.frame(cluster_data$statistics[[metric]]))
    
    p_val <- test_df %>%
      filter(peak == selected_peak) %>%
      pull(!!sym(test)) %>%
      { if (length(.) == 0) NA_real_ else . }
    
    pvalue_list[[cluster_id]] <- p_val
    
    p_label <- sprintf(
      "%s, p = %s",
      test_names[[test]],
      ifelse(is.na(p_val), "NA", sprintf("%.5f", p_val))
    )
    
    p <- ggplot(long_df, aes(x = group, y = value, fill = group)) +
      geom_bar(stat = "summary", fun = mean, position = "dodge", width = 0.4, alpha = 0.8) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.1) +
      geom_jitter(
        color = "black", fill = "black",
        width = 0.15, size = 2, shape = 21, stroke = 0.2,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = c(control = "#A1C9F4", experiment = "#FFB482")) +
      theme_minimal() +
      labs(
        title = paste("Cluster", sub("cluster_", "", cluster_id)),
        y = "Expression", x = "Group"
      ) +
      theme(legend.position = "none")
    
    if (same_y_scale) {
      p <- p +
        coord_cartesian(ylim = c(0, max_y_val * 1.1)) +
        annotate("text", x = 1.5, y = max_y_val * 1.05, label = p_label)
    } else {
      p <- p +
        annotate("text", x = 1.5, y = max(long_df$value, na.rm = TRUE) * 1.05, label = p_label)
    }
    
    plot_list[[cluster_id]] <- p
  }
  
  test_table <- tibble::tibble(
    cluster = names(pvalue_list),
    p_value = unlist(pvalue_list)
  )
  
  wrap_plot_list <- wrap_plots(plot_list, ncol = ncol)
  
  wrap_plot_list
}


#####################
# prepared palletes #
#####################
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

palette_allen <- c("0" = "#ed2b2b",
                     "1" = "#c2f702",
                     "2" = "#f0dd30",
                     "3" ="#911432",
                     "4" = "#a0bd11",
                     "5" = "#54630c", 
                     "6" = "#4e57ba", 
                     "7" = "#0a0a0a", 
                     "8" = "#bd721c", 
                     "9" = "#7ff59f", 
                     "10" = "#b64ef2", 
                     "11" = "#f57878",
                     "12" = "#2059d4", 
                     "13" = "#0f20d9", 
                     "14" = "#578a65",
                     "15" = "#e0c8d2",
                     "16" = "#8bc1f0",
                     "17" = "#f5a556",
                     "18" = "#f77d05",
                     "19" = "#24b0bf",
                     "20" = "#5c95e6",
                     "21" = "#f711dd",
                     "22" = "#8721c2",
                     "23" = "#2cdb12",
                     "24" = "#8a8987",
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


