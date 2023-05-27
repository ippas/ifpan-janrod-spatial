#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

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

# add information about samples
spatial_transcriptomic_data$sample_information <- meta_data

# add information about barcode
spatial_transcriptomic_data$bcs_information <- bcs_merge

# add tibble with images
spatial_transcriptomic_data$images_information <- images_tibble

# add list with raw data contain:
# dataframe with metadata for peaks
# spatial_transcriptomic_data$raw_data$metadata <- integrated_analysis@meta.data 
spatial_transcriptomic_data$raw_data$metadata <- integrated_analysis@meta.data %>% .[, 2:3] %>%
   rownames_to_column(var = "sample_barcode") %>%
   separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
   left_join(., meta_data, by = c("sample" = "sample_ID"))

# dataframe with annotate peaks
spatial_transcriptomic_data$raw_data$annotate <- info_peaks[match(str_replace_all(rownames(integrated_analysis@assays$RNA@counts ), "_", "-"), 
                                                                  str_replace_all(info_peaks$peak_id, "_", "-")),] %>%
   # check strange agree between coverage and gene
   mutate(peak_id = str_replace_all(peak_id, "_", "-")) # %>% to remove, because strand peak agree with gene in each peak 
   # mutate(strand_agree = as.numeric(paste(strand_coverage, "1", sep = "")) == strand_gene)

# sparse matrix contain expression data
spatial_transcriptomic_data$raw_data$data <- integrated_analysis@assays$RNA@counts 

# create dataframe with simple statistic mean and median
spatial_transcriptomic_data$raw_data$simple_statistics <- data.frame(median = pbapply(spatial_transcriptomic_data$raw_data$data, 1, median)) %>% 
   mutate(mean = pbapply(spatial_transcriptomic_data$raw_data$data, 1, mean))


# add list with filter data contain:
# create vector with index to filter
# spatial_transcriptomic_data$filtered_data$filter_index <- which(spatial_transcriptomic_data$raw_data$annotate$strand_agree & spatial_transcriptomic_data$raw_data$simple_statistics$mean > 0.05)

spatial_transcriptomic_data$filtered_data$filter_index <- which(spatial_transcriptomic_data$raw_data$simple_statistics$mean > 0.05)

spatial_transcriptomic_data$filtered_data$metadata <- spatial_transcriptomic_data$raw_data$metadata

spatial_transcriptomic_data$filtered_data$annotate <- spatial_transcriptomic_data$raw_data$annotate[spatial_transcriptomic_data$filtered_data$filter_index, ]

spatial_transcriptomic_data$filtered_data$data <- spatial_transcriptomic_data$raw_data$data[spatial_transcriptomic_data$filtered_data$filter_index, ]


# add list with filter data contain:
spatial_transcriptomic_data$colfilt_data$col_index <- which(apply(spatial_transcriptomic_data$filtered_data$data, 2, function(x){sum(x > 2)}) > 0)

spatial_transcriptomic_data$colfilt_data$metadata <- spatial_transcriptomic_data$filtered_data$metadata[spatial_transcriptomic_data$colfilt_data$col_index,]

spatial_transcriptomic_data$colfilt_data$annotate <- spatial_transcriptomic_data$filtered_data$annotate

spatial_transcriptomic_data$colfilt_data$data <- spatial_transcriptomic_data$filtered_data$data[, spatial_transcriptomic_data$colfilt_data$col_index]

# # add normalize data
# range_normalize <- function(x, range = 500) {
#    wh <- which(x > 1)
#    order <- order(x[wh], decreasing = T)
#    len <- length(wh)
#    out <- x
#    out[wh[order[1:(len - 1)]]] <- (len):2 
#    out
# }

range_normalize <- function(x, range = 1500, flatten = 1) {
   wh <- which(x > flatten)
   order <- order(x[wh], decreasing = T)
   out <- x
   len <- length(wh)
   rangecount <- min(range, len)
   maxrange <- range
   minrange <- range - rangecount + 2
   c(rangecount, maxrange, minrange)
   out[out > flatten] <- flatten
   out[wh[order[(rangecount - 1)]]] <- maxrange:minrange 
   out
}

threshold = 500

spatial_transcriptomic_data$range_normalize$data <- apply(
   as.matrix(spatial_transcriptomic_data$colfilt_data$data), 1, 
   range_normalize) %>% t

colnames(spatial_transcriptomic_data$range_normalize$data) <- colnames(spatial_transcriptomic_data$colfilt_data$data)

spatial_transcriptomic_data$range_normalize$metadata <- spatial_transcriptomic_data$colfilt_data$metadata 

spatial_transcriptomic_data$range_normalize$annotate <- spatial_transcriptomic_data$colfilt_data$annotate



spatial_transcriptomic_data$seurat$data <- integrated_analysis@assays$RNA@data 

spatial_transcriptomic_data$seurat$annotate <- spatial_transcriptomic_data$raw_data$annotate

spatial_transcriptomic_data$seurat$metadata <- spatial_transcriptomic_data$raw_data$metadata


spatial_transcriptomic_data$raw_data$metadata %>%
   .[, c(1, 2, 7, 8, 9)] -> cluster_df

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




