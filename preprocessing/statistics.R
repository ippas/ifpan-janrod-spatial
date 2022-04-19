#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################
install.packages("doParallel")

require(doParallel)
require(rlang)


spatial_transcriptomic_data$raw_data$metadata %>%
  .[, c(1, 2, 7, 8, 9)] -> cluster_df


# prepare data frame with cluster for differ resolution value
for(resolution in seq(0.1, 2, 0.1)){
  
  tmp_column_name <- paste("cluster_resolution", resolution, sep = "_")
  
  
  cluster_df <- FindClusters(integrated_analysis, resolution = resolution) %>%
    .$seurat_clusters %>%
    as.data.frame() %>%
    rename({{tmp_column_name}} := ".") %>%
    rownames_to_column(var = "sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    right_join(cluster_df, ., by = c("barcode", "sample"))
  rm(tmp_column_name)
  
}

# add cluser to 
spatial_transcriptomic_data$clusters <- cluster_df




# # old function to calculate statystic
# stat <- function(x, treat = colfilt.info$treatment, clusters = as.factor(colfilt.info$integrated_snn_res.2)) {
#   idx = 1
#   p <- vector()
#   for(cluster in levels(clusters)) {
#     p[idx] <- 1
#     if (sum(clusters == cluster & treat == "ctrl") < 2 |
#         sum(clusters == cluster & treat == "treat") < 2) {
#       
#     } else {
#       p[idx] <- wilcox.test(x[clusters == cluster & treat == "ctrl"], x[clusters == cluster & treat == "treat"])$p.value
#     }
#     
#     idx = idx + 1
#   }
#   names(p) <- levels(clusters)
#   p
# }



#####################################################
# function to calculate statistic for interest peak #
#####################################################
statistics_interest_peak <- function(spatial_data,
                                   peak,
                                   stat_peak,
                                   stat_test,
                                   resolution,
                                   type_modyfication,
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
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group2) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    spatial_data[[type_modyfication]]$data[peak,barcode_group1] -> peak_group1
    spatial_data[[type_modyfication]]$data[peak, barcode_group2] -> peak_group2
    

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


start_time <- Sys.time()
statistics_interest_peak(spatial_data = spatial_transcriptomic_data,
                       experiment1 = "risperidone",
                       stat_test = "wilcox.test",
                       group1 = "risperidone",
                       group2 = "saline",
                       peak = "merged-samples-peak-94135",
                       type_modyfication = "range_normalize",
                       resolution = 0.1)
end_time <- Sys.time()

end_time - start_time


###############################
# statistic for interest gene #
###############################
statistics_interest_gene <- function(spatial_data,
                                     gene,
                                     type_modyfication,
                                     ...){
  spatial_data[[type_modyfication]]$annotate %>%
    filter(gene_name == gene) %>%
    select(peak_id) %>% .[,1] -> peaks_vector
  
  pvalue_list <- list()
  

  for(peak in peaks_vector){
    pvalue_list[[peak]] <- statistics_interest_peak(spatial_data = spatial_data,
                                                    type_modyfication = type_modyfication,
                                                    peak = peak,
                                                    ...)
  }
  
  pvalue_list
  
}



spatial_interest_cluster(cluster = 27,
                         seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = samples_name,
                         size= 1)


statistics <- function(spatial_data,
                       stat_test,
                       resolution,
                       type_modyfication,
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
  pvalue_list <- list()
  
  for (cluster in clusters_vector){
    
    # prepare barcode for group1
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% sample_vector_group2) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    spatial_data[[type_modyfication]]$data[ ,barcode_group1] -> peak_group1
    spatial_data[[type_modyfication]]$data[, barcode_group2] -> peak_group2
    
    pvalue_vector <- vector()
    
    
    if(stat_test == "t.test"){
      for(index in 1:nrow(peak_group1)){
        pvalue_vector[index] <- t.test(peak_group1[index, ],
                                       peak_group2[index, ],
                                       var.equal = TRUE)$p.value
      }
    } else if(stat_test == "wilcox.test"){
      for(index in 1:nrow(peak_group1)){
        pvalue_vector[index] <- t.test(peak_group1[index,],
                                            peak_group2[index,])$p.value
      }
    }

    pvalue_list[[paste("cluster", cluster, sep = "_")]] <- pvalue_vector

  }
  
  pvalue_matrix <- pvalue_list %>% as.data.frame() %>% as.matrix()
  rownames(pvalue_matrix) <- rownames(peak_group1)
  pvalue_matrix
  
}



start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
           experiment1 = "risperidone",
           stat_test = "wilcox.test",
           group1 = "risperidone",
           group2 = "saline",
           type_modyfication = "range_normalize",
           resolution = 0.1) -> pvalue_matrix
end_time <- Sys.time()

end_time - start_time

filter_gene_statistic <- function(spatial_data,
                                  type_modyfication,
                                  pvalue_data,
                                  gene){
  spatial_data[[type_modyfication]]$annotate %>% 
    filter(gene_name == gene) %>% 
    select(peak_id) %>% .[, 1] -> peaks_vector
  
  pvalue_data[peaks_vector,]
}

filter_gene_statistic(spatial_data = spatial_transcriptomic_data,
                      type_modyfication = "range_normalize",
                      pvalue_data = pvalue_matrix,
                      gene = "Homer1")




statistics_data <- list()

statistic_automate_resolution <- function(spatial_data,
                                          statistics_data,
                                          experiment1,
                                          experiment2,
                                          group1,
                                          group2){
  
  comparision <- paste(paste(experiment1, group1, sep = "."),
                       paste(experiment2, group2, sep = "."),
                       sep = "_")
  
  resolution_vector <- seq(0.1, 2, 0.1)
  # type_modification_vector <- c("raw_data", "range_normalize", "seurat")
  # type_modification_vector <- c("raw_data")
  # stat_test_vector <- c("t.test", "wilcox.test")
  
  
  registerDoParallel(16)
  foreach(resolution=resolution_vector)%dopar% {
    type_modification_vector <- c("raw_data", "range_normalize", "seurat")
    foreach(type_modification = type_modification_vector, .final = function(x) setNames(x, type_modification_vector))%do% {
      stat_test_vector <- c("t.test", "wilcox.test")
      foreach(stat_test = stat_test_vector, .final = function(x) setNames(x, stat_test_vector))%do% {
         # paste(resolution, type_modification, stat_test)
         statistics(spatial_data = spatial_data,
           experiment1 = experiment1,
           experiment2 = experiment2,
           stat_test = stat_test,
           group1 = group1,
           group2 = group2,
           type_modyfication = type_modification,
           resolution = resolution)
      }
    }
  }
}

start_time <- Sys.time()
statistic_automate_resolution(spatial_data = spatial_transcriptomic_data,
                              statistics_data = statistics_data,
                              experiment1 = "risperidone",
                              experiment2 = "risperidone",
                              group1 = "risperidone",
                              group2 = "saline") -> statistics_data

end_time <- Sys.time()

end_time - start_time




start_time <- Sys.time()
for(resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)){
  print(resolution)
  statistics(spatial_data = spatial_transcriptomic_data,
             experiment1 = "risperidone",
             stat_test = "t.test",
             group1 = "risperidone",
             group2 = "saline",
             type_modyfication = "range_normalize",
             resolution = resolution) -> pvalue_matrix
}
end_time <- Sys.time()

end_time - start_time

start_time <- Sys.time()
registerDoParallel(4)
foreach(resolution = c(0.1), .final = function(x) setNames(x, c(0.1))) %:% 
  # type_modification_vector <- c("raw_data", "range_normalize")
  foreach(type_modification =  c("raw_data", "range_normalize"), .final = function(x) setNames(x, c("raw_data", "range_normalize"))) %dopar% {
  # for(type_modification in type_modification_vector){
    # paste(resolution, type_modification)
    statistics(spatial_data = spatial_transcriptomic_data,
               experiment1 = "risperidone",
               stat_test = "t.test",
               group1 = "risperidone",
               group2 = "saline",
               type_modyfication = type_modification,
               resolution = resolution)
  }
  

end_time <- Sys.time()

end_time - start_time


foreach(i=c("test", "test1"), .final = function(x) setNames(x, c("a", "b"))) %:%
  foreach(j=c("some", "some2"), .final = function(x) setNames(x, c("name1", "name2"))) %dopar% {
    type_modification_vector <- c("raw_data", "range_normalize", "seurat")
    paste(i, j)
  }

start_time <- Sys.time()
registerDoParallel(16)
foreach(row = iter(tmp_stat_info, by='row'), 
        .final = function(x) setNames(x, 
                                      paste(tmp_stat_info$resolution, 
                                            tmp_stat_info$type_modification, 
                                            tmp_stat_info$test, 
                                            sep = "-"))) %dopar% {
                                              statistics(spatial_data = spatial_transcriptomic_data,
                                                         experiment1 = "risperidone",
                                                         stat_test = row[,3],
                                                         group1 = "risperidone",
                                                         group2 = "saline",
                                                         type_modyfication = row[,2],
                                                         resolution = row[,1]) -> statistics_data
                                              }
end_time <- Sys.time()

end_time - start_time



start_time <- Sys.time()
registerDoParallel(16)
foreach(row = iter(tmp_stat_info, by='row')) %dopar% {
                                              statistics(spatial_data = spatial_transcriptomic_data,
                                                         experiment1 = "risperidone",
                                                         stat_test = row[,3],
                                                         group1 = "risperidone",
                                                         group2 = "saline",
                                                         type_modyfication = row[,2],
                                                         resolution = row[,1])
                                            }
end_time <- Sys.time()

end_time - start_time





data.frame(resolution = rep(seq(0.1, 2, 0.1), each = 6), 
           type_modification = rep(rep(c("raw_data", "range_normalize", "seurat"), each = 2), 20), 
           test = rep(c("t.test", "wilcox.test"), 60)) %>%
  filter(test == "t.test") -> tmp_stat_info





statistics(spatial_data = spatial_transcriptomic_data,
           experiment1 = "risperidone",
           stat_test = "t.test",
           group1 = "risperidone",
           group2 = "saline",
           type_modyfication = "range_normalize",
           resolution = 1) -> pvalue_matrix


# prepare function to calculate, statistic for insterest resolution
statistics(spatial_data = spatial_transcriptomic_data,
           experiment1 = "tif-ldopa",
           experiment2 = "tif-ldopa",
           stat_test = "t.test",
           group1 = "saline",
           group2 = "ldopa",
           type_modyfication = "range_normalize",
           resolution = 1) -> ttest_saline_ldopa_range_normalize

start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
           experiment1 = "tif-ldopa",
           experiment2 = "tif-ldopa",
           stat_test = "wilcox.test",
           group1 = "saline",
           group2 = "ldopa",
           type_modyfication = "range_normalize",
           resolution = 1) -> wilcox_saline_ldopa_range_normalize
end_time <- Sys.time()
end_time - start_time




######################################################
# function to calculate statistics between two group #
######################################################
statistics2 <- function(spatial_data,
                        stat_test,
                        type_modyfication,
                        resolution,
                        samples_vector1,
                        samples_vector2){
  
  name_column_resolution <- paste("cluster_resolution",
                                  resolution, 
                                  sep = "_")
  
  spatial_data$clusters %>% 
    select(name_column_resolution) %>% 
    .[,1] %>% levels() -> clusters_vector
  
  # create vector to storage pvalue
  pvalue_list <- list()
  
  for (cluster in clusters_vector){
    
    # prepare barcode for group1
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector1) %>%
      filter((!!sym(name_column_resolution)) %in% cluster) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group1
    
    # prepare barcode for group2
    spatial_transcriptomic_data[[type_modyfication]]$metadata[, c(1:2)] %>% 
      right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
      filter(sample %in% samples_vector2) %>%
      mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
      select(sample_barcode) %>% .[,1] -> barcode_group2
    
    # prepare data to statistic test
    spatial_data[[type_modyfication]]$data[ ,barcode_group1] -> peak_group1
    spatial_data[[type_modyfication]]$data[, barcode_group2] -> peak_group2
    
    pvalue_vector <- vector()
    
    
    if(stat_test == "t.test"){
      for(index in 1:nrow(peak_group1)){
        pvalue_vector[index] <- t.test(peak_group1[index, ],
                                       peak_group2[index, ],
                                       var.equal = TRUE)$p.value
      }
    } else if(stat_test == "wilcox.test"){
      for(index in 1:nrow(peak_group1)){
        pvalue_vector[index] <- t.test(peak_group1[index,],
                                       peak_group2[index,])$p.value
      }
    }
    
    pvalue_list[[paste("cluster", cluster, sep = "_")]] <- pvalue_vector
    
  }
  
  pvalue_matrix <- pvalue_list %>% as.matrix()
  rownames(pvalue_matrix) <- rownames(peak_group1)
  pvalue_matrix

}


samples_vector1 <- meta_data %>% 
  filter(treatment == "saline" & mouse_genotype == "tif-mutant") %>%
  .[, 1]


samples_vector2 <- meta_data %>% 
  filter(treatment == "ldopa" & mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>%
  .[, 1]

# wilcox test for resolution equal 1, and range_normalize
start_time <- Sys.time()
statistics2(spatial_data = spatial_transcriptomic_data,
            resolution = 1,
            type_modyfication = "range_normalize",
            stat_test = "t.test",
            samples_vector1 = samples_vector1,
            samples_vector2 = samples_vector2) -> wilcox_saline_ldopa_range_normalize
end_time <- Sys.time()

end_time - start_time


for(stat_test2 in c("t.test", "wilcox.test")){
  for(type_modyfication2 in c("range_normalize", "seurat")){
    statistics2(spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                type_modyfication = type_modyfication2,
                stat_test = stat_test2,
                samples_vector1 = samples_vector1,
                samples_vector2 = samples_vector2) #-> spatial_transcriptomic_data[[type_modyfication]]$statistics$resolution_1[[stat_test]]
  }
}

spatial_transcriptomic_data$range_normalize$statistics$resolution_1$`wilcox.test` <- wilcox_saline_ldopa_range_normalize

spatial_transcriptomic_data$range_normalize$annotate[, c("peak_id", "gene_name")] %>% head



n_top_statistics <- function(spatial_data,
                             type_modyfication,
                             resolution,
                             cluster,
                             stat_test,
                             n){

  spatial_data[[type_modyfication]]$statistics[[paste("resolution", resolution, sep = "_")]][[stat_test]] %>%
  # wilcox_saline_ldopa_range_normalize  %>% 
    .[order(.[, cluster]), ] %>% 
    head(n) %>% 
    as.data.frame() %>% 
    select(cluster) %>% 
    rownames_to_column(var = "peak_id") %>% 
    left_join(., 
              spatial_data[[type_modyfication]]$annotate[, c("peak_id", "gene_name")], 
              by = "peak_id") %>%
    select(c(peak_id, gene_name, cluster))
  
}

n_top_statistics(spatial_data = spatial_transcriptomic_data,
                 type_modyfication = "range_normalize",
                 stat_test = "wilcox.test",
                 resolution = 1,
                 cluster = "cluster_4",
                 n = 100)



spatial_interest_cluster(cluster = 4,
                         seurat_object = integrated_analysis,
                         spatial_data = spatial_transcriptomic_data,
                         resolution = 1,
                         samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[,1]},
                         size= 1, 
                         ncol = 3)






# fold change
t.test(c(1,2,3),c(4,5,6), var.equal = T)$estimate %>% log2 %>% func %>% abs
