#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################
install.packages("doParallel")

require(doParallel)
require(rlang)


# prepare vector with sample_id to calculate statistcs
samples_vector1 <- meta_data %>% 
  filter(treatment == "saline" & mouse_genotype == "tif-mutant") %>%
  .[, 1]


samples_vector2 <- meta_data %>% 
  filter(treatment == "ldopa" & mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>%
  .[, 1]

# execute t.test for range_normalize
start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
           resolution = 1,
           type_data = "range_normalize",
           stat_test = "t.test",
           per = "sample",
           samples_vector1 = samples_vector1,
           samples_vector2 = samples_vector2,
           save_spatial_data = F) -> spatial_transcriptomic_data2
end_time <- Sys.time()

end_time - start_time


start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
            resolution = 1,
            type_data = "range_normalize",
            stat_test = "log2ratio",
            per = "sample",
            samples_vector1 = samples_vector1,
            samples_vector2 = samples_vector2,
            save_spatial_data = TRUE) -> spatial_transcriptomic_data
end_time <- Sys.time()

end_time - start_time


start_time <- Sys.time()
permuate_save_spatial_data(spatial_data = spatial_transcriptomic_data,
                           resolution = 1,
                           type_data = "range_normalize",
                           permutate = 100,
                           pvalue_threshold = 0.05,
                           samples_vector1 = samples_vector1,
                           samples_vector2 = samples_vector2) -> spatial_transcriptomic_data
end_time <- Sys.time()

end_time - start_time



# extract n top gene from spatial_data
n_top_statistics(spatial_data = spatial_transcriptomic_data,
                 type_data = "range_normalize",
                 stat_test = "t.test",
                 resolution = 1,
                 cluster = 0,
                 per = "sample",
                 n = 40) 


genetoplot = "Egr2"

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "range_normalize",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0.00,
                  max_percentile = 0.99,
                  size = 1,
                  normalization = T,
                  ncol = 3) 

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0,
                  max_percentile = 1,
                  size = 1,
                  normalization = T,
                  ncol = 3) 

spatial_cluster(spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                palette = palette_cluster, 
                size= 1, 
                ncol = 3)

spatial_interest_cluster(cluster = 7,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                         size= 1, 
                         ncol = 3)


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
    # dplyr::select(name_column_resolution) %>% 
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
      dplyr::select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_data[[type_data]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector2) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      dplyr::select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    if(per == "spot"){
      spatial_data[[type_data]]$data[ ,barcode_group1] -> peak_group1
      spatial_data[[type_data]]$data[ ,barcode_group2] -> peak_group2
    } 
    else if(per == "sample"){
      sapply(samples_vector1, 
             function(sample_id) {spatial_data[[type_data]]$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, mean, trim = 0.05)}) -> peak_group1
      
      sapply(samples_vector2, 
             function(sample_id) {spatial_data[[type_data]]$data[ , barcode_group2[grepl(sample_id, barcode_group2)]] %>% apply(., 1, mean, trim = 0.05)}) -> peak_group2
    }
    
    result_vector <- vector()
    
    if(stat_test == "t.test"){
      # result_vector <- apply(spatial_data[[type_data]]$data, 1, function(index) t.test(index[barcode_group1], index[barcode_group2])$p.value)
      for(index in 1:nrow(spatial_data[[type_data]]$data)){
        if(length(unique(peak_group1[index, ][!is.na(peak_group1[index,])])) == 1 & length(unique(peak_group2[index, ][!is.na(peak_group2[index,])])) == 1){
          result_vector[index] <- 1
        } else {
          if(mean(peak_group1[index,], na.rm = TRUE) > 0.8 | mean(peak_group2[index,], na.rm = TRUE) > 0.8){
            result_vector[index] <- t.test(peak_group1[index, ],
                                           peak_group2[index, ],
                                           var.equal = TRUE)$p.value
          } else {
            result_vector[index] <- 1
          }
          
        }
        # result_vector[index] <- t.test(peak_group1[index, ],
        #                                peak_group2[index, ],
        #                                var.equal = TRUE)$p.value
      }
    } else if(stat_test == "log2ratio"){
      for(index in 1:nrow(spatial_data[[type_data]]$data)){
        if(length(unique(peak_group1[index, ][!is.na(peak_group1[index,])])) == 1 & length(unique(peak_group2[index, ][!is.na(peak_group2[index,])])) == 1){
          result_vector[index] <- 0
        } else {
          result_vector[index] <- t.test(peak_group1[index, ],
                                         peak_group2[index, ],
                                         var.equal = TRUE)$estimate %>% 
            log2() %>% diff()
          # {if(diff(.) == 0) {.} else {log2(.)}} %>% diff %>% abs %>% as.vector() 
          # diff() %>% abs %>% {if(. == 0) {print(.)} else {log2(.)}} %>% as.vector()
        }
        # result_vector[index] <- t.test(peak_group1[index, ],
        #                                peak_group2[index, ],
        #                                var.equal = TRUE)$estimate %>% 
        #   diff() %>% abs %>% {if(. == 0) {print(.)} else {log2(.)}} %>% as.vector()
      }
    } else if(stat_test == "wilcox.test"){
      for(index in 1:nrow(peak_group1)){
        result_vector[index] <- t.test(peak_group1[index,],
                                       peak_group2[index,])$p.value
      }
    }
    
    pvalue_list[[paste("cluster", cluster, sep = "_")]] <- result_vector
    
  }
  
  result_matrix <- pvalue_list %>% do.call(rbind, .) %>% t
  rownames(result_matrix) <- rownames(peak_group1)
  
  if(save_spatial_data == TRUE){
    spatial_data[[type_data]]$statistics[[name_column_resolution]][[per]][[stat_test]] <- result_matrix
    spatial_data
  }
  else if(save_spatial_data == FALSE){
    result_matrix
  }
  
}

