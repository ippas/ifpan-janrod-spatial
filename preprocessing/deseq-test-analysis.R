# przeprowadzić DESeq dla clastra 5 -> na odczytach dla spotów, zobaczyć wyniki co wychodzi, znormalizowane dane
# przeprowadzić DESeq dla clastra 5 -> na średnich 


# This function normalizes gene expression data based on quantiles within each cluster for spatial transcriptomics data.
# add_quantile_norm_data <- function(spatial_data, resolution = 0.8, num_cores = 24, data_type = "raw_data"){
#   
#   # This function normalizes gene expression data based on quantiles within each cluster for spatial transcriptomics data.
#   # The normalization is performed on a given resolution level.
#   # 'spatial_data': A list containing spatial data with raw_data (gene expression data, metadata, and annotation data) and cluster data.
#   # 'resolution': The resolution level at which the data will be clustered. Default is 0.8.
#   # 'num_cores': The number of cores to use for parallel processing. Default is 24.
#   # 'data_type': The type of data to be used for normalization. Default is "raw_data".
#   
#   # Function to normalize a given cluster.
#   # 'cluster': The cluster id that should be normalized.
#   quantile_norm_cluster <- function(cluster){
#     
#     # Get all barcodes from the cluster to be normalized.
#     # The !!sym() function allows us to programmatically use the resolution column name.
#     all_barcodes <- spatial_data[[data_type]]$metadata %>%
#       right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
#       filter((!!sym(resolution_column)) == cluster) %>%
#       mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
#       pull(sample_barcode)
#     
#     # Extract the data matrix for the specific cluster.
#     cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] #%>% as.matrix()
#     
#     
#     
#     
#     # Normalize the data for the cluster.
#     norm_data <- normalize.quantiles(cluster_data)
#     
#     # Keep the column and row names identical to the original data.
#     colnames(norm_data) <- colnames(cluster_data)
#     rownames(norm_data) <- rownames(cluster_data)
#     
#     # Return the normalized data.
#     return(norm_data)
#   } 
#   
#   # Define the column name based on the resolution parameter.
#   resolution_column <- paste0("cluster_resolution_", resolution)
#   
#   # Get a sorted list of unique cluster ids at the given resolution.
#   unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
#   
#   # Start the timer to measure the computation time.
#   start_time <- Sys.time()
#   
#   # Apply the quantile_norm_cluster function to each unique cluster in parallel using 'num_cores'.
#   results <- mclapply(unique_clusters, quantile_norm_cluster,  mc.cores = num_cores)
#   
#   # Name each result (i.e., each normalized cluster) by its cluster id.
#   results <- setNames(results, paste0("cluster_", unique_clusters))
#   
#   # Stop the timer and calculate the elapsed time.
#   end_time <- Sys.time()
#   diff_time <- end_time - start_time
#   
#   # Print the computation time.
#   print(diff_time)
#   
#   # Combine the normalized clusters back into a single matrix.
#   quantile_normalize_data <- do.call(cbind, results) %>% Matrix()
#   
#   # Add the normalized data and the original metadata and annotation data back into the spatial data.
#   spatial_data$quantile_normalize$metadata <- spatial_data$raw_data$metadata
#   spatial_data$quantile_normalize$annotate <- spatial_data$raw_data$annotate
#   spatial_data$quantile_normalize$data <- quantile_normalize_data
#   
#   # Return the updated spatial data.
#   return(spatial_data)
# }

#############################################33
spatial_data$sample_information %>% 
  filter(sample_ID %in% c(samples_saline, samples_risperidone)) %>%
  select(c(sample_ID, treatment)) %>%
  dplyr::rename(sample = "sample_ID") -> df_condition



spatial_data$sample_information %>% 
  filter(sample_ID %in% c(samples_saline, samples_risperidone)) %>%
  select(c(sample_ID, treatment)) %>%
  dplyr::rename(sample = "sample_ID") -> df_condition



deseq_norm_cluster(cluster = 1)
##########################################################################
deseq_norm_cluster <- function(cluster){
  
  # Get all barcodes from the cluster to be normalized.
  # The !!sym() function allows us to programmatically use the resolution column name.
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    pull(sample_barcode)
  
  # Extract the data matrix for the specific cluster.
  cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] %>% as.matrix()
  
  data.frame(barcodes =  all_barcodes) %>%
    mutate(sample = str_split(barcodes, "_")) %>% separate(barcodes,
                                                           into = c("sample", "barcodes"),
                                                           sep = "_") %>% mutate(barcodes = paste0(sample, "_", barcodes)) %>%
    left_join(., df_condition, by = "sample") -> col_data
  
  dds <- DESeqDataSetFromMatrix(
    countData =  cluster_data,
    colData = col_data,
    design = ~condition)
  
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  # # normalized_counts <- normalized_counts[all_barcodes]
  # all_barcodes_ordered <- intersect(colnames(normalized_counts), all_barcodes)
  # normalized_counts <- normalized_counts[, all_barcodes_ordered]
  
  # Return the normalized data.
  return(normalized_counts)
} 

start_time <- Sys.time()
deseq_norm_cluster(cluster = 1) -> tmp
end_time <- Sys.time()
start_time - end_time

add_deseq2_norm_data <- function(spatial_data, resolution = 0.8, num_cores = 24, data_type = "raw_data", control_samples, experiment_samples){
  
  spatial_data$sample_information %>% 
    filter(sample_ID %in% c(control_samples, experiment_samples)) %>%
    select(c(sample_ID, treatment)) %>%
    dplyr::rename(sample = "sample_ID") %>%
    dplyr::rename(condition = "treatment") -> df_condition
  
  deseq_norm_cluster <- function(cluster){
    
    # Get all barcodes from the cluster to be normalized.
    # The !!sym() function allows us to programmatically use the resolution column name.
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      pull(sample_barcode)
    
    # Extract the data matrix for the specific cluster.
    cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] %>% as.matrix()
    
    cluster_data <- cluster_data +1
    
    data.frame(barcodes =  all_barcodes) %>%
      mutate(sample = str_split(barcodes, "_")) %>% separate(barcodes,
                                                             into = c("sample", "barcodes"),
                                                             sep = "_") %>% mutate(barcodes = paste0(sample, "_", barcodes)) %>%
      left_join(., df_condition, by = "sample") -> col_data
    
    dds <- DESeqDataSetFromMatrix(
      countData =  cluster_data,
      colData = col_data,
      design = ~condition)
    
    # Run the DESeq pipeline
    dds <- DESeq(dds)
    
    # compute size factors
    dds <- estimateSizeFactors(dds)
    
    normalized_counts <- counts(dds, normalized=TRUE)
    
    # normalized_counts <- normalized_counts[all_barcodes]
    # all_barcodes_ordered <- intersect(colnames(normalized_counts), all_barcodes)
    # normalized_counts <- normalized_counts[, all_barcodes_ordered]
    
    # Return the normalized data.
    return(normalized_counts)
  } 

  
  # Define the column name based on the resolution parameter.
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  # Get a sorted list of unique cluster ids at the given resolution.
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  # Start the timer to measure the computation time.
  start_time <- Sys.time()
  
  # Apply the quantile_norm_cluster function to each unique cluster in parallel using 'num_cores'.
  results <- mclapply(unique_clusters, deseq_norm_cluster,  mc.cores = num_cores)
  
  # Name each result (i.e., each normalized cluster) by its cluster id.
  results <- setNames(results, paste0("cluster_", unique_clusters))
  
  # Stop the timer and calculate the elapsed time.
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  
  # Print the computation time.
  print(diff_time)

  # Combine the normalized clusters back into a single matrix.
  deseq2_normalize_data <- do.call(cbind, results) %>% Matrix()

  # Add the normalized data and the original metadata and annotation data back into the spatial data.
  spatial_data$deseq2_normalize$metadata <- spatial_data$raw_data$metadata
  spatial_data$deseq2_normalize$annotate <- spatial_data$raw_data$annotate
  spatial_data$deseq2_normalize$data <- deseq2_normalize_data

  # Return the updated spatial data.
  return(spatial_data)
  # return(results)
  
}




add_deseq2_norm_data(spatial_data = risperidone_st_data_half,
                     resolution = 0.8,
                     num_cores = 24,
                     data_type = "raw_data",
                     control_samples = samples_saline,
                     experiment_samples = samples_risperidone) -> tmp


spatial_gene_plot(spatial_data = tmp,
                  type_data = "quantile_normalize",
                  gene = "Cmtm5",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T)

tmp$deseq2_normalize$data["risperidone-peak-19962",]
##########################################################################


control_expression <- spatial_data[[data_type]]$data[, results_raw2$cluster_5$control_barcodes]
experiment_expression <- spatial_data[[data_type]]$data[, results_raw2$cluster_5$experiment_barcodes]

countdata <- cbind(control_expression, experiment_expression) 

df_condition <- data.frame(sample = c(samples_saline, samples_risperidone),
                condition = c(rep("saline", 6), rep("risperidone", 6)))

data.frame(
  barcodes = c(
    results_raw2$cluster_5$control_barcodes,
    results_raw2$cluster_5$experiment_barcodes
  )
) %>%
  mutate(sample = str_split(barcodes, "_")) %>% separate(barcodes,
                                                         into = c("sample", "barcodes"),
                                                         sep = "_") %>% mutate(barcodes = paste0(sample, "_", barcodes)) %>%
  left_join(., df_condition, by = "sample") -> coldata


dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(cbind(control_expression, experiment_expression)),
  colData = coldata,
  design = ~condition)
  
# Run the DESeq pipeline
dds <- DESeq(dds)

# compute size factors
dds <- estimateSizeFactors(dds)


normalized_counts <- counts(dds, normalized=TRUE)

# Create boxplots of raw counts
# boxplot(countdata, main = "Raw Counts", las = 2)

# Create boxplots of normalized counts
plotMA(dds, ylim=c(-2,2))

# First, apply variance stabilizing transformation
vst_data <- varianceStabilizingTransformation(dds)

# Now, plot the PCA
plotPCA(vst_data, intgroup = "condition")


res <- results(dds)

res %>% 
  as.data.frame() %>%
  filter(padj < 0.0001) %>%
  rownames_to_column(var = "peak_id") %>% 
  left_join(., risperidone_st_data_half$raw_data$annotate %>%
              select(c(peak_id, gene_name)), by = "peak_id") -> deseq_signif

risperidone_st_data_half$raw_data$annotate %>%
  select(c(peak_id, gene_name))


perform_enrichment_analysis(genes = deseq_genes, database = "Mouse_Gene_Atlas", top_rows)

perform_enrichment_analysis(genes = deseq_genes, database = "GO_Molecular_Function_2021", top_rows)

enrichr_database %>% filter(libraryName == "KEGG_2021_Human")


normalized_counts %>% .[1:10, 1:20]

mean_sample_cluster <- function(data, samples, min_spot = 0, trim = 0){
  
  barcodes <- data %>% colnames()
  
  sapply(samples, function(sample_id) {
    # Selecting the control data based on the sample ID
    sample_data <- data[, barcodes[grepl(paste0(sample_id, "_"), barcodes)]]
    if (ncol(sample_data) >= min_spots) { 
      # Calculate the mean when the number of columns is greater than or equal to the minimum spots 
      apply(sample_data, 1, mean, trim = trim) 
    } else { 
      # Return a vector of NAs with names from the control data rows when the number of columns is less than minimum spots
      setNames(rep(NA, nrow(sample_data)), rownames(sample_data)) 
    }
  })
}

heatmap_signif_genes <- function(data, signif_data){
  data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "peak_id") %>%
    right_join(., signif_data[, c("peak_id", "gene_name")]) %>%
    mutate(gene_peak = paste0(gene_name, "_", peak_id)) %>%
    select(-c(peak_id, gene_name)) %>% 
    column_to_rownames(var = "gene_peak") %>%
    as.matrix() %>%
    gplots::heatmap.2(.,
                      distfun = function(x) as.dist(1-cor(t(x))),
                      # col=rev(morecols(50)),trace="none",
                      Colv = FALSE,
                      scale="row",
                      sepwidth = c(0.3,0.3),
                      # labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
                      srtCol = 0,
                      # cexRow = sizeRow,
                      offsetCol = 0.1,
                      ColSideColors = c(rep("red", 6), rep("blue", 6)),
                      symkey=FALSE,
                      cexCol = 1.1,
                      col = colorRampPalette(c("blue", "white", "red"))(n = 100),
                      # lhei = lhei,
                      margins = c(4, 7))
}


mean_sample_cluster(data = normalized_counts, samples =  c(control_samples, experiment_samples)) -> norm_deseq_mean

heatmap_signif_genes <- function(data, signif_data){
  data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "peak_id") %>%
    right_join(., signif_data[, c("peak_id", "gene_name")]) %>%
    mutate(gene_peak = paste0(gene_name, "_", peak_id)) %>%
    select(-c(peak_id, gene_name)) %>% 
    column_to_rownames(var = "gene_peak") %>%
    as.matrix() %>%
    gplots::heatmap.2(.,
                      distfun = function(x) as.dist(1-cor(t(x))),
                      # col=rev(morecols(50)),trace="none",
                      Colv = FALSE,
                      scale="row",
                      sepwidth = c(0.3,0.3),
                      # labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
                      srtCol = 0,
                      # cexRow = sizeRow,
                      offsetCol = 0.1,
                      ColSideColors = c(rep("red", 6), rep("blue", 6)),
                      symkey=FALSE,
                      cexCol = 1.1,
                      col = colorRampPalette(c("blue", "white", "red"))(n = 100),
                      # lhei = lhei,
                      margins = c(4, 7))
}

heatmap_signif_genes(data = norm_deseq_mean, signif_data = deseq_signif)

heatmap_signif_genes(
  data = cbind(
    results_raw2$cluster_0$control_expression,
    results_raw2$cluster_0$experiment_expression
  ),
  signif_data = deseq_signif
)


tmp <- cbind(
  results_raw2$cluster_0$control_expression,
  results_raw2$cluster_0$experiment_expression
)

quantiles_mean <- normalize.quantiles(tmp)

colnames(quantiles_mean) <- colnames(tmp)
rownames(quantiles_mean) <- rownames(tmp)

heatmap_signif_genes(data = norm_deseq_mean, signif_data = deseq_signif)

heatmap_signif_genes(data = quantiles_mean, signif_data = deseq_signif)




# visualization norm_data 




{results_raw2$cluster_5$control_expression *100} %>% 
  as.data.frame() %>%
  mutate_all(function(x) as.integer(x))


## comparison on mean from raw data
control_expression <- {results_raw2$cluster_5$control_expression *100} %>% 
  as.data.frame() %>%
  mutate_all(function(x) as.integer(x))


experiment_expression <- {results_raw2$cluster_5$experiment_expression * 100} %>% 
  as.data.frame() %>%
  mutate_all(function(x) as.integer(x))

countdata <- cbind(control_expression, experiment_expression) %>% as.matrix()

df_condition <- data.frame(sample = c(samples_saline, samples_risperidone),
                           condition = c(rep("saline", 6), rep("risperidone", 6)))

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = df_condition,
  design = ~condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

# Create boxplots of raw counts
# boxplot(countdata, main = "Raw Counts", las = 2)

# Create boxplots of normalized counts
plotMA(dds, ylim=c(-2,2))

# First, apply variance stabilizing transformation
vst_data <- varianceStabilizingTransformation(dds)

# Now, plot the PCA
plotPCA(vst_data, intgroup = "condition")


res <- results(dds)

res %>% 
  as.data.frame() %>%
  filter(pvalue < 0.01) %>%
  filter(log2FoldChange > 1) %>%
  rownames_to_column(var = "peak_id") %>% 
  left_join(., risperidone_st_data_half$raw_data$annotate %>%
              select(c(peak_id, gene_name)), by = "peak_id")

