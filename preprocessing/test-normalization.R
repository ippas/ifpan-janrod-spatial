######################
# Normalization data #
######################
results <- data.frame(median = pbapply(sample_data_all, 1, median))
results$mean <- pbapply(sample_data_all, 1, mean)

sample_anno_all$strand_agree <- as.numeric(paste(sample_anno_all$strand_coverage, "1", sep = "")) == sample_anno_all$strand_gene

filtered_wh <- which(sample_anno_all$strand_agree & results$mean > 0.03)
filtered_anno <- sample_anno_all[filtered_wh,]
filtered_data <- sample_data_all[filtered_wh,]
filtered_info <- sample_info_all

filtered_info$over_2 <- apply(filtered_data, 2, function(x){sum(x > 2)})

threshold <- 200

colfilt_wh <- which(filtered_info$over_2 > 0)
colfilt_anno <- filtered_anno
colfilt_data <- filtered_data[,colfilt_wh]
colfilt_info <- filtered_info[colfilt_wh,]



normalize <- function(x, range = threshold) {
  order <- order(x, decreasing = T)
  out <- rep(0, length(x))
  out[order[1:range]] <- range:1
  out
}

range_normalize <- function(x, range = 2000) {
  wh <- which(x > 1)
  order <- order(x[wh], decreasing = T)
  len <- length(wh)
  out <- x
  out[wh[order[1:(len - 1)]]] <- (len):2 
  out
}

colfilt_norm_data <- apply(colfilt_data, 1, normalize, threshold) %>% t

colnames(colfilt_norm_data) <- colnames(colfilt_data)



resolution = 0.2
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)


data_cluster <- FindClusters(integrated_analysis, resolution = resolution) %>% 
  .$seurat_clusters %>%
  as.data.frame() %>% 
  rename(cluster = ".") %>% 
  rownames_to_column(var = "sample_barcode") %>%
  separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  left_join(., bcs_merge, by = c("barcode", "sample"))

plot_clusters(data_cluster)



colfilt_anno %>%
  filter(gene_name == "Jun") %>%
  .[,4] -> peak_gene_vector

for(peak in peak_gene_vector){
  plot_feature(data_cluster = data_cluster,
               peak_id = peak,
               size = 1) %>% print()
  rm(peak)
}



# test sparse matrix
sample_data_all %>%  object.size()

sample_data_all %>% Matrix(., sparse = TRUE) %>% object.size()

spatial_list <- list()

spatial_list$sample <- samples_name

spatial_list$raw_data <- list()
spatial_list$raw_data$metadata <- sample_info_all
spatial_list$raw_data$annotate <- sample_info_all
spatial_list$raw_data$data <- sample_data_all %>% Matrix(., sparse = TRUE)
spatial_list$raw_data$simply_stat <- results <- data.frame(median = pbapply(sample_data_all, 1, median)) %>% 
  mutate(mean = pbapply(sample_data_all, 1, mean))

filter_peaks <- function(data, metric, threshold){
  which(data$raw_data$sample_anno_all$strand_agree & data$raw_data$simply_stat$mean > threshold)
}

filter_peaks(spatial_list, "mean", 0.03)


filtered_wh <- which(sample_anno_all$strand_agree & results$mean > 0.03)
filtered_anno <- sample_anno_all[filtered_wh,]
filtered_data <- sample_data_all[filtered_wh,]
filtered_info <- sample_info_all

filtered_info$over_2 <- apply(filtered_data, 2, function(x){sum(x > 2)})

threshold <- 200

colfilt_wh <- which(filtered_info$over_2 > 0)
colfilt_anno <- filtered_anno
colfilt_data <- filtered_data[,colfilt_wh]
colfilt_info <- filtered_info[colfilt_wh,]


sample_data_all_sparse %>%
  pbapply(., 1, mean)
  .[1:10, 1:10] %>% apply(., 1, median)
  .["merged-samples-peak-165",]

  
  
sample_anno_all %>% filter(gene_name == "Egr2")
sample_data_all["merged-samples-peak-19993",] %>% table



## prepare object to store data to spatial transcriptomic
# data contain:
# sample name - vector
# tible witm images
# dataframe contain inforamation about barcode information
# list with annotate peaks, metadata, and data from expresion
# dataframe with simply statistic: media and mean needed to normalization
# list for normalization data: annotate peak, metadata, and data from expression after normalization

spatial_transcriptomic_data <- list()

# add vector of samples 
spatial_transcriptomic_data$samples <- samples_name




spatial_list$sample <- samples_name

spatial_list$raw_data <- list()
spatial_list$raw_data$metadata <- sample_info_all
spatial_list$raw_data$annotate <- sample_info_all
spatial_list$raw_data$data <- sample_data_all %>% Matrix(., sparse = TRUE)
spatial_list$raw_data$simply_stat <- results <- data.frame(median = pbapply(sample_data_all, 1, median)) %>% 
  mutate(mean = pbapply(sample_data_all, 1, mean))




