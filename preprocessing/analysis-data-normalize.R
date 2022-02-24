#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

meta_data <- read.table("data/samples-spatial-metadata.tsv",
                        header = TRUE,
                        sep = "\t")

# Require source files
info_peaks <- read.table("data/gene-annotation/peaks-annotate-sort.bed", 
                         header = FALSE,
                         sep = "\t",
                         col.names = c('chr_peak',
                                       'start_peak',
                                       'end_peak',
                                       'peak_id',
                                       'score_int(-10*log10pvalue)',
                                       'strand_coverage',
                                       'fold_change_peak_summit',
                                       '-log10pvalue_peak_summit',
                                       '-log10qvalue_peak_summit',
                                       'relative_summit_position_peak_start',
                                       'type_peak',
                                       'chr_gene',
                                       'start_gene',
                                       'end_gene',
                                       'gene_id',
                                       'gene_name',
                                       'strand_gene')) 



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



