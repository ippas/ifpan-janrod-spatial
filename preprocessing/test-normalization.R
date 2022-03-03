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

# add information about barcode
spatial_transcriptomic_data$bcs_information <- bcs_merge

# add tibble with images
spatial_transcriptomic_data$images_information <- images_tibble

# add list with raw data contain:
# dataframe with metadata for peaks
spatial_transcriptomic_data$raw_data$metadata <- integrated_analysis@meta.data %>% 
  left_join(., meta_data, by = c("sample" = "sample_ID"))

# dataframe with annotate peaks
spatial_transcriptomic_data$raw_data$annotate <- info_peaks[match(str_replace_all(rownames(integrated_analysis@assays$RNA@counts ), "_", "-"), 
                                                                  str_replace_all(info_peaks$peak_id, "_", "-")),] %>%
  # check strange agree between coverage and gene
  mutate(peak_id = str_replace_all(peak_id, "_", "-")) %>% 
  mutate(strand_agree = as.numeric(paste(strand_coverage, "1", sep = "")) == strand_gene)

# sparse matrix contain expression data
spatial_transcriptomic_data$raw_data$data <- integrated_analysis@assays$RNA@counts 

# create dataframe with simple statistic mean and median
spatial_list$raw_data$simple_statistics <- data.frame(median = pbapply(spatial_transcriptomic_data$raw_data$data, 1, median)) %>% 
  mutate(mean = pbapply(spatial_transcriptomic_data$raw_data$data, 1, mean))


# add list with filter data contain:
# create vector with index to filter 
spatial_transcriptomic_data$filtered_data$filter_index <- which(spatial_transcriptomic_data$raw_data$annotate$strand_agree & spatial_list$raw_data$simple_statistics$mean > 0.05)

spatial_transcriptomic_data$filtered_data$metadata <- spatial_transcriptomic_data$raw_data$metadata 

spatial_transcriptomic_data$filtered_data$annotate <- spatial_transcriptomic_data$raw_data$annotate[spatial_transcriptomic_data$filtered_data$filter_index, ]

spatial_transcriptomic_data$filtered_data$data <- spatial_transcriptomic_data$raw_data$data[spatial_transcriptomic_data$filtered_data$filter_index, ]


# add list with filter data contain:
spatial_transcriptomic_data$colfilt_data$col_index <- which(apply(spatial_transcriptomic_data$filtered_data$data, 2, function(x){sum(x > 2)}) > 0)

spatial_transcriptomic_data$colfilt_data$metadata <- spatial_transcriptomic_data$filtered_data$metadata[spatial_transcriptomic_data$colfilt_data$col_index,]

spatial_transcriptomic_data$colfilt_data$annotate <- spatial_transcriptomic_data$filtered_data$annotate

spatial_transcriptomic_data$colfilt_data$data <- spatial_transcriptomic_data$filtered_data$data[, spatial_transcriptomic_data$colfilt_data$col_index]

# add normalize data
range_normalize <- function(x, range = 500) {
  wh <- which(x > 1)
  order <- order(x[wh], decreasing = T)
  len <- length(wh)
  out <- x
  out[wh[order[1:(len - 1)]]] <- (len):2 
  out
}

threshold = 500

spatial_transcriptomic_data$range_normalize$data <- apply(
  spatial_transcriptomic_data$colfilt_data$data, 1, 
  normalize) %>% t

colnames(spatial_transcriptomic_data$range_normalize$data) <- colnames(spatial_transcriptomic_data$colfilt_data$data)



################
# creating PCA #
################
pca_normalize_data <- prcomp(spatial_transcriptomic_data$range_normalize$data, center = TRUE,scale. = TRUE)

pca_var <- pca_normalize_data$sdev^2
pve <- pca_var / sum(pca_var)

plot(pve , xlab = "Principal Component",
     ylab = "Proportion of Variance Expained",
     ylim = c(0, 1), type = "b")


pca_normalize_data$rotation %>% .[1:10, 1:10]


pca_normalize_data$rotation %>% .[1:10, 1:1000]

install.packages("umap")
require(umap)

umap_data <- umap(pca_normalize_data$rotation %>% .[, 1:20])



umap_data$layout %>% 
  as.data.frame() %>% 
  rename(umap1 = "V1",
         umap2 = "V2") %>% 
  # filter(umap1 < 5,
  #        umap2 > -5) %>%
  kmeans(., 15) %>% .$cluster %>% 
  as.data.frame() %>% rename(cluster = ".")  -> umap_cluster

umap_data$layout %>% 
as.data.frame() %>% 
rename(umap1 = "V1",
       umap2 = "V2") %>% 
# filter(umap1 < 5,
#        umap2 > -5) %>%
  cbind(., umap_cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>%
ggplot(aes(x = umap1, y = umap2, col = cluster)) + 
geom_point()


data_cluster %>% select(-cluster) %>% head %>%
  left_join(., umap_cluster)

umap_cluster %>%
  rownames_to_column(var = "sample_barcode") %>%
  separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  left_join(., data_cluster[, -3], by = c("barcode", "sample")) -> data_cluster_kmeans

plot_clusters(data_cluster_kmeans)




#####################################
# test seurat method to integration #
#####################################


