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
                                 type_data,
                                 peak_id,
                                 samples,
                                 min_percentile = 0.01,
                                 max_percentile = 0.99,
                                 size = 1.2,
                                 normalization = TRUE){
  
  
  # prepare data to create plot
  preprocessing_data <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    rename(value = ".") %>% 
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
           })
  
  plot_list <- list()
  
  # creating plot
  for (sample_id in samples){
    sample_plot <- filter(preprocessing_data, sample == sample_id) %>%
      # main part of creating plot
      {ggplot(., aes(x = imagecol, y = imagerow, fill = value)) +
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
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
                        filter(sample == sample_id) %>%
                        select(width)))+
          ylim(max(spatial_data$bcs_information %>%
                     filter(sample == sample_id) %>%
                     select(height)), 0) +
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
    
    
    plot_list[[sample_id]] <- sample_plot
    
  }
  
  wrap_plot <- wrap_plots(plot_list)
  
  rm(plot_list,
     preprocessing_data)
  
  return(wrap_plot)
  
}



#########################################
# prepare function to visualize cluster #
#########################################

spatial_cluster <- function(spatial_data,
                            resolution,
                            # seurat_object,
                            samples,
                            palette, 
                            size= 1.2,
                            ncol){
  
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  # prepare data for plots
  data_cluster <- spatial_data$clusters %>% 
    select(sample, barcode, name_column_resolution) %>%
    rename(cluster = name_column_resolution) %>%
    left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
  # data_cluster <- FindClusters(seurat_object, resolution = resolution) %>% 
  #   .$seurat_clusters %>%
  #   as.data.frame() %>% 
  #   rename(cluster = ".") %>% 
  #   rownames_to_column(var = "sample_barcode") %>%
  #   separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  #   left_join(., spatial_data$bcs_information, by = c("barcode", "sample"))
  
  plot_list <- list()

  # prepare list of plots
  for (sample_id in samples){
    sample_plot <- filter(data_cluster, sample == sample_id) %>%
      {ggplot(., aes(x = imagecol, y = imagerow, fill = factor(cluster))) +
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
                       aes(grob = grob),
                       x = 0.5, y = 0.5) +
          geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
          coord_cartesian(expand = FALSE) +
          scale_fill_manual(values = palette)+
          xlim(0, max(spatial_data$bcs_information %>%
                        filter(sample == sample_id) %>%
                        select(width)))+
          ylim(max(spatial_data$bcs_information %>%
                     filter(sample == sample_id) %>%
                     select(height)), 0) +
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


# example to run function
spatial_cluster(spatial_data = spatial_transcriptomic_data,
                # seurat_object = integrated_analysis,
                resolution = 1,
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
  mypalette[cluster + 1] <- "#00FF00"
  
  spatial_cluster(palette = mypalette, 
                  ...)
  
}

# example to run function
spatial_interest_cluster(cluster = 2,
                         # seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = samples_name,
                         size= 1,
                         ncol = 4)

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
    select(peak_id) %>% .[,1] -> vector_peak
  
  for(peak in vector_peak){
    spatial_feature_plot(spatial_data = spatial_data,
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



# example run function
spatial_feature_plot(spatial_data = spatial_transcriptomic_data,
                     type_data = "range_normalize",
                     peak_id = "merged-samples-peak-19995",
                     samples = samples_name,
                     min_percentile = 0,
                     max_percentile = 1,
                     size = 1,
                     normalization = TRUE)




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
                                         clusters,
                                         normalization = TRUE){
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  
  data_cluster <- spatial_data$clusters %>% 
    select(sample, barcode, name_column_resolution) %>%
    rename(cluster = name_column_resolution)
  
  # prepare data to create plot
  preprocessing_data <- spatial_data[[type_data]]$data[peak_id, ] %>%
    as.data.frame() %>%
    rename(value = ".") %>% 
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
          geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id,],
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
                        filter(sample == sample_id) %>%
                        select(width)))+
          ylim(max(spatial_data$bcs_information %>%
                     filter(sample == sample_id) %>%
                     select(height)), 0) +
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
    select(peak_id) %>% .[,1] -> vector_peak
  
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


spatial_gene_plot_cluster(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = genetoplot,
                  samples = samples_name,
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  normalization = T,
                  resolution = 1,
                  clusters = c(4),
                  ncol = 4) 

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





#####################################################
# function to calculate statistic for interest peak #
#####################################################
statistics_interest_peak <- function(spatial_data,
                                     peak,
                                     stat_peak,
                                     stat_test,
                                     resolution,
                                     type_data,
                                     experiment1,
                                     experiment2 = NULL,
                                     group1,
                                     group2){
  
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  
  spatial_data$clusters %>% 
    select(name_column_resolution) %>% 
    .[,1] %>% levels() -> clusters_vector
  
  if(is.null(experiment2)){
    experiment2 <- experiment1
  }
  
  # samples for group1
  sample_vector_group1 <- spatial_data$clusters %>%
    filter(experiment == experiment1,
           treatment == group1) %>%
    .[, 1] %>% unique()
  
  
  # samples for group2
  sample_vector_group2 <- spatial_data$clusters %>%
    filter(experiment == experiment2,
           treatment == group2) %>% 
    unique() %>%
    .[, 1]
  
  # create vector to storage pvalue
  pvalue_vector <- vector()
  
  for (cluster in clusters_vector){
    
    # prepare barcode for group1
    spatial_transcriptomic_data[[type_data]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_transcriptomic_data[[type_data]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group2) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    spatial_data[[type_data]]$data[peak,barcode_group1] -> peak_group1
    spatial_data[[type_data]]$data[peak, barcode_group2] -> peak_group2
    
    
    if(stat_test == "t.test"){
      pvalue_vector[cluster] <- t.test(peak_group1,
                                       peak_group2,
                                       var.equal = TRUE)$p.value
    } else if(stat_test == "wilcox.test"){
      pvalue_vector[cluster] <- wilcox.test(peak_group1, peak_group2)$p.value
    }
    
  }
  
  pvalue_vector
  
}

# example to run function
start_time <- Sys.time()
statistics_interest_peak(spatial_data = spatial_transcriptomic_data,
                         experiment1 = "risperidone",
                         stat_test = "wilcox.test",
                         group1 = "risperidone",
                         group2 = "saline",
                         peak = "merged-samples-peak-94135",
                         type_data = "range_normalize",
                         resolution = 0.1)
end_time <- Sys.time()

end_time - start_time
##

#####################################################
# function to calculate statistic for interest gene #
#####################################################
statistics_interest_gene <- function(spatial_data,
                                     gene,
                                     type_data,
                                     ...){
  spatial_data[[type_data]]$annotate %>%
    filter(gene_name == gene) %>%
    select(peak_id) %>% .[,1] -> peaks_vector
  
  pvalue_list <- list()
  
  
  for(peak in peaks_vector){
    pvalue_list[[peak]] <- statistics_interest_peak(spatial_data = spatial_data,
                                                    type_data = type_data,
                                                    peak = peak,
                                                    ...)
  }
  
  pvalue_list
  
}

#########################
# function for statistics
#
#########################
############################################
# prepare function to calculate statistics #
############################################
statistics <- function(spatial_data,
                       stat_test,
                       type_data,
                       resolution,
                       per,
                       samples_vector1,
                       samples_vector2,
                       save_spatial_data = FALSE){
  
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  
  spatial_data$clusters[name_column_resolution] %>%
    # select(name_column_resolution) %>% 
    .[,1] %>% levels() -> clusters_vector
  
  # create vector to storage pvalue
  pvalue_list <- list()
  
  for (cluster in c(clusters_vector)){
    
    # prepare barcode for group1
    spatial_data[[type_data]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_data[[type_data]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector2) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    if(per == "spot"){
      spatial_data[[type_data]]$data[ ,barcode_group1] -> peak_group1
      spatial_data[[type_data]]$data[ ,barcode_group2] -> peak_group2
    } 
    else if(per == "sample"){
      sapply(samples_vector1, 
             function(sample_id) {spatial_data[[type_data]]$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, mean)}) -> peak_group1
      
      sapply(samples_vector2, 
             function(sample_id) {spatial_data[[type_data]]$data[ , barcode_group2[grepl(sample_id, barcode_group2)]] %>% apply(., 1, mean)}) -> peak_group2
    }
    
    pvalue_vector <- vector()
    
    if(stat_test == "t.test"){
      # pvalue_vector <- apply(spatial_data[[type_data]]$data, 1, function(index) t.test(index[barcode_group1], index[barcode_group2])$p.value)
      for(index in 1:nrow(spatial_data[[type_data]]$data)){
        if(length(unique(peak_group1[index, ][!is.na(peak_group1[index,])])) == 1 & length(unique(peak_group2[index, ][!is.na(peak_group2[index,])])) == 1){
          pvalue_vector[index] <- 1
        } else {
          pvalue_vector[index] <- t.test(peak_group1[index, ],
                                         peak_group2[index, ],
                                         var.equal = TRUE)$p.value
        }
        # pvalue_vector[index] <- t.test(peak_group1[index, ],
        #                                peak_group2[index, ],
        #                                var.equal = TRUE)$p.value
      }
    } else if(stat_test == "fold_change"){
      for(index in 1:nrow(spatial_data[[type_data]]$data)){
        pvalue_vector[index] <- t.test(peak_group1[index, ],
                                       peak_group2[index, ],
                                       var.equal = TRUE)$estimate %>% 
          diff() %>% abs %>% {if(. == 0) {print(.)} else {log2(.)}} %>% as.vector()
      }
    } else if(stat_test == "wilcox.test"){
      for(index in 1:nrow(peak_group1)){
        pvalue_vector[index] <- t.test(peak_group1[index,],
                                       peak_group2[index,])$p.value
      }
    }
    
    pvalue_list[[paste("cluster", cluster, sep = "_")]] <- pvalue_vector
    
  }
  
  pvalue_matrix <- pvalue_list %>% do.call(rbind, .) %>% t
  rownames(pvalue_matrix) <- rownames(peak_group1)
  
  if(save_spatial_data == TRUE){
    spatial_data[[type_data]]$statistics[[name_column_resolution]][[per]][[stat_test]] <- pvalue_matrix
    spatial_data
  }
  else if(save_spatial_data == FALSE){
    pvalue_matrix
  }
  
}

# example to run function
samples_vector1 <- meta_data %>% 
  filter(treatment == "saline" & mouse_genotype == "tif-mutant") %>%
  .[, 1]


samples_vector2 <- meta_data %>% 
  filter(treatment == "ldopa" & mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>%
  .[, 1]

start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
           resolution = 1,
           type_data = "range_normalize",
           stat_test = "t.test",
           per = "sample",
           samples_vector1 = samples_vector1,
           samples_vector2 = samples_vector2,
           save_spatial_data = TRUE) -> spatial_transcriptomic_data
end_time <- Sys.time()


################################################
# prepare function to extract n top statistics #
################################################
n_top_statistics <- function(spatial_data,
                             type_data,
                             resolution,
                             cluster,
                             stat_test,
                             per,
                             n){
  
  
  # fold_change <- spatial_data[[type_data]]$statistics[[paste("cluster_resolution", resolution, sep = "_")]][[per]]$fold_change %>%
  #   as.data.frame() %>%
  #   select(cluster) %>%
  #   rownames_to_column(var = "peak_id") %>%
  #   rename(fold_change = cluster)
  
  name_cluster <- paste("cluster", cluster, sep = "_")
  
  out <- spatial_data[[type_data]]$statistics[[paste("cluster_resolution", resolution, sep = "_")]][[per]][[stat_test]] %>%
    # wilcox_saline_ldopa_range_normalize  %>% 
    .[order(.[, name_cluster]), ] %>% 
    as.data.frame() %>% 
    select(name_cluster) %>% 
    rownames_to_column(var = "peak_id") %>% 
    left_join(., 
              spatial_data[[type_data]]$annotate[, c("peak_id", "gene_name")], 
              by = "peak_id") %>%
    filter(get(name_cluster) != 1)
    out$FDRO <- p.adjust(out[,name_cluster], method = "BH")
    #head(n) %>%
    # left_join(.,
    #           fold_change,
    #           by = "peak_id") %>%
    # select(c(peak_id, gene_name, cluster, FDR, fold_change))
    out %>% head(n) %>%
    select(c(peak_id, gene_name, name_cluster, "FDRO"))  
  
}

# example to run function
n_top_statistics(spatial_data = spatial_transcriptomic_data,
                 type_data = "range_normalize",
                 stat_test = "t.test",
                 resolution = 1,
                 cluster = 3,
                 per = "sample",
                 n = 40) 