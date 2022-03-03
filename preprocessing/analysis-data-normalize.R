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
spatial_transcriptomic_data$raw_data$simple_statistics <- data.frame(median = pbapply(spatial_transcriptomic_data$raw_data$data, 1, median)) %>% 
   mutate(mean = pbapply(spatial_transcriptomic_data$raw_data$data, 1, mean))


# add list with filter data contain:
# create vector with index to filter 
spatial_transcriptomic_data$filtered_data$filter_index <- which(spatial_transcriptomic_data$raw_data$annotate$strand_agree & spatial_transcriptomic_data$raw_data$simple_statistics$mean > 0.05)

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
   range_normalize) %>% t

colnames(spatial_transcriptomic_data$range_normalize$data) <- colnames(spatial_transcriptomic_data$colfilt_data$data)



# For all features
sample_data_all <- integrated_analysis@assays$RNA@counts %>% as.matrix()

# sample.data.all.normalize <- integrated.analysis.cluster@assays$RNA@data %>% as.matrix()

sample_anno_all <- info_peaks[match(str_replace_all(rownames(sample_data_all), "_", "-"), 
                                    str_replace_all(info_peaks$peak_id, "_", "-")),] %>%
  mutate(peak_id = str_replace_all(peak_id, "_", "-"))

sample_info_all <- integrated_analysis@meta.data %>% 
  left_join(., meta_data, by = c("sample" = "sample_ID"))




######################
# Normalization data #
######################
results <- data.frame(median = pbapply(sample_data_all, 1, median))
results$mean <- pbapply(sample_data_all, 1, mean)

sample_anno_all$strand_agree <- as.numeric(paste(sample_anno_all$strand_coverage, "1", sep = "")) == sample_anno_all$strand_gene

filtered_wh <- which(sample_anno_all$strand_agree & results$mean > 0.05)
filtered_anno <- sample_anno_all[filtered_wh,]
filtered_data <- sample_data_all[filtered_wh,]
filtered_info <- sample_info_all

filtered_info$over_2 <- apply(filtered_data, 2, function(x){sum(x > 2)})

threshold <- 500

colfilt_wh <- which(filtered_info$over_2 > 0)
colfilt_anno <- filtered_anno
colfilt_data <- filtered_data[,colfilt_wh]
colfilt_info <- filtered_info[colfilt_wh,]

range_normalize <- function(x, range = 500) {
  wh <- which(x > 1)
  order <- order(x[wh], decreasing = T)
  len <- length(wh)
  out <- x
  out[wh[order[1:(len - 1)]]] <- (len):2 
  out
}

colfilt_norm_data <- apply(colfilt_data, 1, range_normalize, threshold) %>% t

colnames(colfilt_norm_data) <- colnames(colfilt_data)

rm(threshold,
   colfilt_anno,
   colfilt_info,
   colfilt_data,
   filtered_anno,
   filtered.data,
   filtered.info,
   filtered.wh,
   colfilt.wh
   )



