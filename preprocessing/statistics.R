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


start_time <- Sys.time()
statistics(spatial_data = spatial_transcriptomic_data,
            resolution = 1,
            type_data = "seurat",
            stat_test = "t.test",
            per = "sample",
            samples_vector1 = samples_vector1,
            samples_vector2 = samples_vector2,
            save_spatial_data = TRUE) -> spatial_transcriptomic_data
end_time <- Sys.time()

end_time - start_time



start_time <- Sys.time()
for(index in 1:nrow({stat_arg_df %>% filter(type_data == "range_normalize")})){
  row <- {stat_arg_df %>% filter(type_data == "range_normalize")} %>% .[index,]
  
  print(row)
  statistics(spatial_data = spatial_transcriptomic_data,
              resolution = 1,
              type_data = row[,1],
              stat_test = row[,3],
              per = row[,2],
              samples_vector1 = samples_vector1,
              samples_vector2 = samples_vector2,
              save_spatial_data = TRUE) -> spatial_transcriptomic_data
  
}
end_time <- Sys.time()

end_time - start_time



# extract n top gene from spatial_data
n_top_statistics(spatial_data = spatial_transcriptomic_data,
                 type_data = "seurat",
                 stat_test = "t.test",
                 resolution = 1,
                 cluster = 0,
                 per = "sample",
                 n = 40) 

for(number in 0:24) {
  n_top_statistics(spatial_data = spatial_transcriptomic_data,
                   type_data = "range_normalize",
                   stat_test = "t.test",
                   resolution = 1,
                   cluster = number,
                   per = "sample",
                   n = 20) %>% print
}


genetoplot = "Junb"

spatial_gene_plot(spatial_data = spatial_transcriptomic_data,
                  type_data = "raw_data",
                  gene = genetoplot,
                  samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                  min_percentile = 0.0,
                  max_percentile = 1,
                  size = 1,
                  normalization = T,
                  ncol = 3) 


spatial_cluster(spatial_data = spatial_transcriptomic_data,
                resolution = 1,
                # seurat_object = integrated_analysis,
                samples = samples_name,
                palette = palette_cluster, 
                size= 1, 
                ncol = 4)

spatial_interest_cluster(cluster = 0,
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



sapply(samples_vector1, function(sample_id) {spatial_transcriptomic_data$range_normalize$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, mean)} )

# create data.frame with condition to calculate statistics
rep(c("raw_data", "range_normalize", "seurat"), each = 6)
rep(rep(c("spot", "sample"), each = 3), 3)
rep(c("t.test", "fold_change", "wilcox.test"), 6) 

data.frame(
  type_data = rep(c("raw_data", "range_normalize", "seurat"), each = 6),
  per = rep(rep(c("spot", "sample"), each = 3), 3),
  test = rep(c("t.test", "fold_change", "wilcox.test"), 6) 
) -> stat_arg_df

# Example for uses foreach to paralize statistics
start_time <- Sys.time()
registerDoParallel(16)
foreach(row = iter(stat_arg_df, by='row'), .final = function(x) setNames(x, {stat_arg_df %>% mutate(name = paste(type_data, per, test, sep = "_")) %>% .[,4]})) %dopar% {
  
  statistics(spatial_data = spatial_transcriptomic_data,
              resolution = 0.1,
              type_data = row[,1],
              stat_test = row[,3],
              per = row[,2],
              samples_vector1 = samples_vector1,
              samples_vector2 = samples_vector2) 
} -> statistics_data
end_time <- Sys.time()

end_time - start_time




######################################################
statistics(spatial_data = spatial_transcriptomic_data,
           resolution = 1,
           type_data = "range_normalize",
           stat_test = "t.test",
           per = "sample",
           samples_vector1 = samples_vector1,
           samples_vector2 = samples_vector2,
           save_spatial_data = TRUE) -> spatial_transcriptomic_data_mean

statistics(spatial_data = spatial_transcriptomic_data_mean,
           resolution = 1,
           type_data = "seurat",
           stat_test = "t.test",
           per = "sample",
           samples_vector1 = samples_vector1,
           samples_vector2 = samples_vector2,
           save_spatial_data = TRUE) -> spatial_transcriptomic_data_mean

n_top_statistics(spatial_data = spatial_transcriptomic_data_mean,
                 type_data = "range_normalize",
                 stat_test = "t.test",
                 resolution = 1,
                 cluster = 0,
                 per = "sample",
                 n = 40) 

# code for extract info for insterst gene
junb <- "merged-samples-peak-201965"
type_data <- "range_normalize"
resolution <- 1
name_column_resolution <- paste("cluster_resolution", resolution, sep = "_")
cluster <- 0


spatial_transcriptomic_data[[type_data]]$metadata[, c(1:2)] %>% 
  right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
  filter(sample %in% samples_vector1) %>%
  filter((!!sym(name_column_resolution)) %in% cluster) %>%
  mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
  select(sample_barcode) %>% .[,1] -> barcode_group1

# prepare barcode for group2
spatial_transcriptomic_data[[type_data]]$metadata[, c(1:2)] %>% 
  right_join(spatial_transcriptomic_data$clusters, ., by = c("barcode", "sample")) %>%
  filter(sample %in% samples_vector2) %>%
  filter((!!sym(name_column_resolution)) %in% cluster) %>%
  mutate(sample_barcode = paste(.$sample, .$barcode, sep = "_")) %>%
  select(sample_barcode) %>% .[,1] -> barcode_group2

sapply(samples_vector1, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, median)}) -> peak_group1
peak_group1 %>% .[junb,]
 
sapply(samples_vector2, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group2[grepl(sample_id, barcode_group2)]] %>% apply(., 1, median)}) -> peak_group2
peak_group2 %>% .[junb,]



sapply(samples_vector1, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, mean)}) -> peak_group1
peak_group1 %>% .[junb,]

sapply(samples_vector2, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group2[grepl(sample_id, barcode_group2)]] %>% apply(., 1, mean)}) -> peak_group2
peak_group2 %>% .[junb,]


sapply(samples_vector1, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group1[grepl(sample_id, barcode_group1)]] %>% apply(., 1, length)}) -> peak_group1
peak_group1 %>% .[junb,]

sapply(samples_vector2, 
       function(sample_id) {spatial_transcriptomic_data[[type_data]]$data[ , barcode_group2[grepl(sample_id, barcode_group2)]] %>% apply(., 1, length)}) -> peak_group2
peak_group2 %>% .[junb,]



