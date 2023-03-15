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
           save_spatial_data = TRUE) -> spatial_transcriptomic_data
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






