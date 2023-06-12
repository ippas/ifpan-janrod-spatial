# data_type <-  "quantile_normalize"
# resolution_column <- "cluster_resolution_0.1"
# spatial_data <- risperidone_st_data_half
# control_samples <- samples_saline
# experiment_samples <- samples_risperidone
# cluster <- 2


# Function to plot UMI reads per cluster for a given sample
plot_umi_reads_cluster <- function(spatial_data, data_type, resolution, cluster) {
  
  # create a string to represent the column for the desired resolution
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  # Get the unique clusters sorted
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  # Join metadata and clusters data frames on "barcode" and "sample"
  # Filter for the desired cluster
  # Create a new column "sample_barcode" that is the combination of "sample" and "barcode"
  # Select the "sample_barcode" column and transform it into a character vector
  all_barcodes <- spatial_data[[data_type]]$metadata %>%
    right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>% 
    filter((!!sym(resolution_column)) == cluster) %>%
    mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
    select(sample_barcode) %>% .[,1]
  
  
  # Subsets the UMI count matrix to include only barcodes present in the cluster
  # Calculates the total UMI counts per cell
  cell_umi_counts <- spatial_data[[data_type]]$data[, all_barcodes] %>%
    as.matrix() %>%
    apply(2, sum)
  
  # Convert cell UMI counts to a dataframe
  # Rename the column to "ncounts"
  # Split the "sample_barcode" column into "sample" and "barcode"
  # Join with the sample_information dataframe to get treatment information
  df <- cell_umi_counts %>%
    as.data.frame() %>%
    dplyr::rename(ncounts = ".")  %>%
    rownames_to_column(var = "sample_barcode") %>%
    separate(sample_barcode, into = c("sample", "barcode"), sep = "_") %>%
    left_join(., spatial_data$sample_information, by = c("sample" = "sample_ID"))
  
  # Calculate and print mean and median for each sample
  sample_stats <- df %>%
    group_by(sample, treatment) %>%
    summarize(mean = mean(ncounts), median = median(ncounts), .groups = 'drop') %>%
    as.data.frame() %>% 
    arrange(treatment)
  
  
  # Create a boxplot of UMI counts per sample, colored by treatment
  ggplot(df, aes(x = sample, y = ncounts, fill = treatment)) +
    geom_boxplot() +
    labs(title = paste0("UMI counts per sample: ", data_type), x = "Sample", y = "UMI count") +
    theme_minimal()
  
}


# plot_umi_reads_cluster(spatial_data = risperidone_st_data_half, 
#                        data_type = "raw_data",
#                        resolution = 0.1,
#                        cluster = 1) -> p1
# 
# plot_umi_reads_cluster(spatial_data = risperidone_st_data_half, 
#                        data_type = "quantile_normalize_resolution_0.1",
#                        resolution = 0.1,
#                        cluster = 1) -> p2
# p1/p2
