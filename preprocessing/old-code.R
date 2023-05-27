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

  for (cluster in c(clusters_vector[c(1:25)])){

    # prepare barcode for group1
    spatial_data[[type_data]]$metadata[, c(1:2)] %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1

    # prepare barcode for group2
    spatial_data[[type_data]]$metadata[, c(1:2)] %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
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


# execute t.test for range_normalize
start_time <- Sys.time()
statistics(spatial_data = risperidone_st_data_half,
           resolution = 1,
           type_data = "range_normalize",
           stat_test = "t.test",
           per = "sample",
           samples_vector1 = samples_saline,
           samples_vector2 = samples_risperidone,
           save_spatial_data = F) -> tmp
end_time <- Sys.time()

end_time - start_time


tmp %>% colnames() -> colnames_vector

sapply(colnames_vector, function(x){tmp %>% as.data.frame() %>% filter(get(x) < 0.05) %>% rownames()}) -> signif_peaks_list

for(peaks_vector in signif_peaks_list){
  risperidone_st_data_half$raw_data$annotate %>% 
    select(c(gene_name, peak_id)) %>% 
    filter(peak_id %in% peaks_vector) %>% print
} 

sapply(signif_peaks_list, function(x){
  spatial_transcriptomic_data$raw_data$annotate %>% 
    select(c(gene_name, peak_id)) %>% 
    filter(peak_id %in% x) %>% select(gene_name) %>% .[1]
}) -> signif_genes_list

spatial_interest_cluster(cluster = 4,
                         # seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = c(samples_vector1, samples_vector2),
                         size= 1,
                         ncol = 4)

# visualize features 
spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = "Fos",
                  samples =  c(samples_vector1, samples_vector2),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)

tmp %>% as.data.frame() %>% filter(cluster_16 < 0.05) %>% rownames()


spatial_transcriptomic_data$raw_data$annotate %>% select(c(gene_name, peak_id)) %>% 
  filter(peak_id %in% c("risperidone-peak-54964", "risperidone-peak-104847", "risperidone-peak-110939", "risperidone-peak-116990", "risperidone-peak-150996", "risperidone-peak-185604"))

DimPlot(data_cluster, reduction = "umap", 
        ncol = 4)

# visualize clusters
spatial_cluster(spatial_data = spatial_transcriptomic_data,
                resolution = 0.5,
                samples = c(samples_vector1, samples_vector2),
                palette = palette_allen, 
                size= 1.0, 
                ncol = 4)


install.packages("ape")
library(ape)
integrated_analysis <- integrated_analysis %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(integrated_analysis, ndims = 50)


test_dims <- function(integrated_analysis, dims_range) {
  integrated_analysis <- integrated_analysis %>% 
    RunUMAP(reduction = "pca", dims = dims_range) %>%
    FindNeighbors(reduction = "pca", dims = dims_range) %>%
    FindClusters(resolution = 0.5)
  
  integrated_analysis <- BuildClusterTree(integrated_analysis)
  return(integrated_analysis)
}

# Test different ranges of dimensions
cluster_tree_1_to_10 <- test_dims(integrated_analysis, 1:10)
cluster_tree_1_to_20 <- test_dims(integrated_analysis, 1:20)
cluster_tree_1_to_30 <- test_dims(integrated_analysis, 1:30)
# Add more tests if needed

PlotClusterTree(cluster_tree_1_to_10)
PlotClusterTree(cluster_tree_1_to_20)
PlotClusterTree(cluster_tree_1_to_30)

# Visualize the cluster trees
plot(cluster_tree_1_to_10)
plot(cluster_tree_1_to_20)
plot(cluster_tree_1_to_30)


compute_data_summary <-
  function(spatial_data,
           resolution = 0.8,
           trim = 0.05,
           num_cores = 24,
           control_samples,
           experiment_samples,
           data_type = "quantile_normalize") {
    
    resolution_column <- paste0("cluster_resolution_", resolution)
    unique_clusters <-
      unique(spatial_data$clusters[[resolution_column]]) %>% sort
    
    compute_data_summary_cluster <- function(cluster) {
      all_barcodes <- spatial_data[[data_type]]$metadata %>%
        right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
        filter((!!sym(resolution_column)) == cluster) %>%
        mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
        mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
        select(group, sample_barcode)
      
      # Split barcodes into control and experiment groups based on group label
      barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
      
      # Extract control and experiment barcodes from the list
      control_barcodes <- barcodes_list[["control"]]
      experiment_barcodes <- barcodes_list[["experiment"]]
      
      control_expression_spot <- sapply(control_samples, function(sample_id) {
        # Selecting the control data based on the sample ID
        control_data <- spatial_data[[data_type]]$data[, control_barcodes[grepl(paste0(sample_id, "_"), control_barcodes)]]
      })
      
      control_sum <- sapply(control_expression_spot, function(x){apply(x, 1, sum)})
      control_mean <- sapply(control_expression_spot, function(x){apply(x, 1, mean, trim = trim, na.rm = TRUE)})
      # control_geometric_mean <- sapply(control_expression_spot, function(x){apply(x, 1, psych::geometric.mean)})
      control_median <- sapply(control_expression_spot, function(x){apply(x, 1, median)})
      control_IQR <- sapply(control_expression_spot, function(x){apply(x, 1, IQR)})
      control_diff_range <- sapply(control_expression_spot, function(x){apply(x, 1, function(x){range(x) %>% diff()})})
      control_var <- sapply(control_expression_spot, function(x){apply(x, 1, var)})
      control_skewness <- sapply(control_expression_spot, function(x){apply(x, 1, skewness)})
      control_kurtosis <- sapply(control_expression_spot, function(x){apply(x, 1, kurtosis)})
      
      control <- list(expression_spot = control_expression_spot,
                      sum = control_sum,
                      mean = control_mean,
                      # geometric_mean = control_geometric_mean,
                      median = control_median,
                      IQR = control_IQR,
                      diff_range = control_diff_range,
                      var = control_var,
                      skewness = control_skewness,
                      kurtosis = control_kurtosis)
      
      rm(
        control_expression_spot, control_mean, control_median, control_IQR,
        control_dff_range, control_var, control_skewness, control_kurtosis,
        control_geometric_mean, control_sum
      )
      
      experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
        # Selecting the experiment data based on the sample ID
        experiment_data <- spatial_data[[data_type]]$data[, experiment_barcodes[grepl(paste0(sample_id, "_"), experiment_barcodes)]]
      })
      
      experiment_sum <- sapply(experiment_expression_spot, function(x){apply(x, 1, sum)})
      experiment_mean <- sapply(experiment_expression_spot, function(x){apply(x, 1, mean, trim = trim, na.rm = TRUE)})
      # experiment_geometric_mean <- sapply(experiment_expression_spot, function(x){apply(x, 1, psych::geometric.mean)})
      experiment_median <- sapply(experiment_expression_spot, function(x){apply(x, 1, median)})
      experiment_IQR <- sapply(experiment_expression_spot, function(x){apply(x, 1, IQR)})
      experiment_diff_range <- sapply(experiment_expression_spot, function(x){apply(x, 1, function(x){range(x) %>% diff()})})
      experiment_var <- sapply(experiment_expression_spot, function(x){apply(x, 1, var)})
      experiment_skewness <- sapply(experiment_expression_spot, function(x){apply(x, 1, skewness)})
      experiment_kurtosis <- sapply(experiment_expression_spot, function(x){apply(x, 1, kurtosis)})
      
      
      experiment <- list(expression_spot = experiment_expression_spot,
                         sum = experiment_sum,
                         mean = experiment_mean,
                         # geometric_mean = experiment_geometric_mean,
                         median = experiment_median,
                         IQR = experiment_IQR,
                         diff_range = experiment_diff_range,
                         var = experiment_var,
                         skewness = experiment_skewness,
                         kurtosis = experiment_kurtosis)
      
      rm(
        experiment_expression_spot, experiment_mean, experiment_median, experiment_IQR,
        experiment_diff_range, experiment_var, experiment_skewness, experiment_kurtosis,
        experiment_geometric_mean, experiment_sum
      )
      
      return(list(peak = spatial_data[[data_type]]$annotate$peak_id,
                  gene = spatial_data[[data_type]]$annotate$gene_name,
                  control = control,
                  experiment = experiment
      ))
    }
    
    start_time <- Sys.time()
    # Use mclapply instead of lapply
    results <- mclapply(unique_clusters, compute_data_summary_cluster, mc.cores = num_cores)
    # Set names for the results
    results <- setNames(results, paste0("cluster_", unique_clusters))
    diff_time <- end_time <- Sys.time()
    print(diff_time)
    
    return(results)
  }
