# Function to install and load packages
install_and_load <- function(package_name, bioc = FALSE, github = FALSE) {
  if (!is.element(package_name, installed.packages()[,1])) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(package_name)
    } else if (github) {
      if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
      remotes::install_github(package_name)
    } else {
      install.packages(package_name)
    }
  }
  require(package_name, character.only = TRUE)
}

#install_and_load <- function(package_name) {
#  if (!is.element(package_name, installed.packages()[,1])) {
#    message(paste("Package", package_name, "not found"))
#  } else {
#    success <- require(package_name, character.only = TRUE)
#    if (!success) {
#      message(paste("Failed to load", package_name))
#    }
#  }
#}


# List of CRAN packages
cran_packages <- c("ggplot2", "Matrix", "rjson", "cowplot", "RColorBrewer", "Seurat", 
                   "grid", "readbitmap", "dplyr", "data.table", "doSNOW", "hdf5r",
                   "remotes", "BiocManager", "tidyr", "tibble", "gridExtra", "parallel", 
                   "preprocessCore", "enrichR", "gplots", "e1071", "psych", "DT", "shinydashboard")

# List of Bioconductor packages
bioc_packages <- c("rhdf5", "DESeq2", "edgeR")

# List of GitHub packages
github_packages <- c("satijalab/seurat-data")

# Install and load CRAN packages
for (package in cran_packages) {
  install_and_load(package)
}

# Install and load Bioconductor packages
for (package in bioc_packages) {
  install_and_load(package, bioc = TRUE)
}

# # Install and load GitHub packages
# for (package in github_packages) {
#   install_and_load(package, github = TRUE)
# }

# Load additional packages if necessary
# additional_packages <- c("SeuratData", "patchwork", "pbapply", "stringr")
additional_packages <- c("patchwork", "pbapply", "stringr")
for (package in additional_packages) {
  if (!require(package, character.only = TRUE)) {
    stop(paste("Package", package, "not found."))
  }
}

##########################################################
# functions to spatial transcriptomics analysis - seurat #
##########################################################

# function for seurat analysis
integrate_data_seurat <- function(path_to_data, nfeatures = 2000, dims = 1:30) {
  # This function integrates multiple single-cell RNA-seq datasets using the Seurat package in R.
  # It reads and processes each sample, merges them, normalizes the data, finds variable features,
  # identifies integration anchors, and finally integrates the data. After data integration,
  # it performs downstream analysis, including scaling, PCA, UMAP, and neighbor finding.
  # The function returns a Seurat object containing the integrated analysis.
  #
  # Input variables:
  # path_to_data: The path to the directory containing the single-cell RNA-seq datasets.
  # nfeatures: The number of variable features to consider for integration (default: 2000).
  # dims: The principal components to use for downstream analysis (default: 1:30).
  
  library(Seurat)
  
  # Get list of samples
  samples_name <- list.files(path = path_to_data)
  
  # Read and process each sample
  samples_list <- lapply(samples_name, function(sample_name) {
    directory <- paste(path_to_data, sample_name, "/outs/filtered_feature_bc_matrix/", sep = "")
    sample_data <- Read10X(data.dir = directory)
    sample <- CreateSeuratObject(counts = sample_data, project = sample_name)
    AddMetaData(sample, metadata = sample_name, col.name = "sample")
  })
  
  # Merge samples
  merged_samples <- merge(samples_list[[1]], samples_list[-1], add.cell.ids = samples_name, project = "merged.sample")
  
  # Split the dataset into a list
  merged_samples_list <- SplitObject(merged_samples, split.by = "sample")
  
  # Normalize and identify variable features for each dataset independently
  merged_samples_list <- lapply(merged_samples_list, function(x) {
    x %>%
      NormalizeData(normalization.method = "LogNormalize") %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
  })
  
  # Select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = merged_samples_list, nfeatures = nfeatures)
  
  # Find integration anchors
  merged_anchors <- FindIntegrationAnchors(object.list = merged_samples_list, anchor.features = features)
  
  # Create an 'integrated' data assay
  integrated_data <- IntegrateData(anchorset = merged_anchors)
  
  # Set the default assay to the corrected (integrated) data
  DefaultAssay(integrated_data) <- "integrated"
  
  # Perform downstream analysis on the integrated data
  integrated_analysis <- integrated_data %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = length(dims), verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = dims) %>%
    FindNeighbors(reduction = "pca", dims = dims)
  
  return(integrated_analysis)
}

# Read metadata from a file and filter for specified treatments
read_metadata <- function(file_path, treatments) {
  # Read metadata from a file and filter for specified treatments
  #
  # Args:
  #   file_path: A character string specifying the path to the metadata file
  #   treatments: A character vector of treatment names to filter the metadata
  #
  # Returns:
  #   A data frame containing the filtered metadata
  
  metadata <- read.table(file_path,
                         header = TRUE,
                         sep = "\t") %>%
    filter(treatment %in% treatments)
  
  return(metadata)
}

# Read gene annotation information from a file
read_gene_annotation <- function(file_path) {
  # Read gene annotation information from a file
  #
  # Args:
  #   file_path: A character string specifying the path to the gene annotation file
  #
  # Returns:
  #   A data frame containing the gene annotation information
  
  info_peaks <- read.table(file_path,
                           header = TRUE,
                           sep = "\t")
  
  return(info_peaks)
}

# function to read images for spatial transcriptomics
create_images_tibble <- function(path_to_data, samples_name) {
  # This function reads low-resolution tissue images, converts them into
  # graphical objects (grobs), and stores them in a tibble, along with
  # information about the sample name, image height, and image width.
  #
  # Input parameters:
  # path_to_data: The path to the directory containing the tissue images.
  # samples_name: A vector of sample names used to build the image file paths.
  
  
  # Load the necessary libraries
  library(grid)
  library(tibble)
  library(png)
  library(gridExtra)
  
  # Read low-resolution tissue images and store them in a list
  images_list <- lapply(samples_name, function(x) {
    file_path <- paste(path_to_data, x, "/outs/spatial/tissue_lowres_image.png", sep = "")
    readPNG(file_path)
  })
  
  # Convert the images to grobs (graphical objects) and store them in a tibble
  images_tibble <- lapply(images_list, function(image) {
    rasterGrob(image, width = unit(1, "npc"), height = unit(1, "npc"))
  }) %>%
    tibble(sample = factor(samples_name), grob = .) %>%
    mutate(height = rapply(images_list, function(x) data.frame(height = nrow(x))),
           width = rapply(images_list, function(x) data.frame(weight = ncol(x))))
  
  return(images_tibble)
}

# function to read barcode spot information
create_barcode_data <- function(path_to_data, sample_names, images_tibble) {
  # This function creates a merged table containing barcode spot information
  # for all samples, including image dimensions and adjusted row/column positions.
  #
  # Input parameters:
  # path_to_data: The path to the directory containing the tissue images and barcode data.
  # samples_name: A vector of sample names used to build the file paths.
  # images_tibble: A tibble containing the image information, including sample name, height, and width.
  
  # Load the necessary libraries
  library(dplyr)
  library(rjson)
  
  # Create a merged table containing barcode spot information for all samples
  barcode_data <- lapply(sample_names, function(sample_name) {
    # Load tissue position data
    tissue_positions <- read.csv(
      paste(
        path_to_data,
        sample_name,
        "/outs/spatial/tissue_positions_list.csv",
        sep = ""
      ),
      col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),
      header = FALSE
    )
    
    # Calculate image row and column positions using scale factors
    scale_factors <-
      paste(path_to_data,
            sample_name,
            "/outs/spatial/scalefactors_json.json",
            sep = "") %>%
      rjson::fromJSON(file = .)
    tissue_lowres_scalef <- scale_factors$tissue_lowres_scalef
    
    # Add image dimensions and adjusted row/column positions to the data
    tissue_positions %>%
      mutate(
        imagerow = imagerow * tissue_lowres_scalef,
        imagecol = imagecol * tissue_lowres_scalef,
        tissue = as.factor(tissue),
        height = images_tibble[images_tibble$sample == sample_name, ]$height,
        width = images_tibble[images_tibble$sample == sample_name, ]$width
      )
  }) %>%
    setNames(sample_names) %>%
    bind_rows(., .id = "sample")
  
  return(barcode_data)
}

# Create a spatial transcriptomic data object
create_spatial_data <- function(sample_names, metadata, barcode_info, images_info, integrated_data, peaks_info) {
  # Create a spatial transcriptomic data object
  #
  # Args:
  #   sample_names: A character vector of sample names
  #   metadata: A data frame containing sample metadata
  #   barcode_info: A data frame containing barcode information
  #   images_info: A tibble with image information
  #   integrated_data: A Seurat object with integrated single-cell data
  #   peaks_info: A data frame containing peak annotation information
  #
  # Returns:
  #   A list containing the spatial transcriptomic data
  
  spatial_data <- list()
  
  # Add vector of samples
  spatial_data$samples <- sample_names
  
  # Add information about samples
  spatial_data$sample_information <- metadata
  
  # Add information about barcode
  spatial_data$bcs_information <- barcode_info
  
  # Add tibble with images
  spatial_data$images_information <- images_info
  
  # Add raw data to spatial data
  spatial_data <- add_raw_data(spatial_data, integrated_data, peaks_info)
  
  return(spatial_data)
}

# Add raw data to the spatial_data object
add_raw_data <- function(spatial_data, integrated_data, peaks_info) {
  # Add raw data to the spatial_data object
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   integrated_data: A Seurat object with integrated single-cell data
  #   peaks_info: A data frame containing peak annotation information
  #
  # Returns:
  #   A list containing the spatial_data object with raw data added
  
  # Add metadata to raw_data
  spatial_data$raw_data$metadata <- integrated_data@meta.data %>% .[, 2:3] %>%
    rownames_to_column(var = "sample_barcode") %>%
    separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
    left_join(., spatial_data$sample_information, by = c("sample" = "sample_ID"))
  
  # Add peak annotations to raw_data
  spatial_data$raw_data$annotate <- peaks_info[match(str_replace_all(rownames(integrated_data@assays$RNA@counts ), "_", "-"),
                                                     str_replace_all(peaks_info$peak_id, "_", "-")),] %>%
    mutate(peak_id = str_replace_all(peak_id, "_", "-"))
  
  # Add raw counts data to raw_data
  spatial_data$raw_data$data <- integrated_data@assays$RNA@counts
  
  # Add simple statistics to raw_data
  spatial_data$raw_data$simple_statistics <- data.frame(median = pbapply(spatial_data$raw_data$data, 1, median)) %>%
    mutate(mean = pbapply(spatial_data$raw_data$data, 1, mean))
  
  return(spatial_data)
}

# Add Seurat data to the spatial_data object
add_seurat_data <- function(spatial_data, integrated_data) {
  # Add Seurat data to the spatial_data object
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   integrated_data: A Seurat object with integrated single-cell data
  #
  # Returns:
  #   A list containing the spatial_data object with Seurat data added
  
  # Add Seurat data to the spatial_data object
  spatial_data$seurat$data <- integrated_data@assays$RNA@data
  
  # Add peak annotations to the Seurat data
  spatial_data$seurat$annotate <- spatial_data$raw_data$annotate
  
  # Add metadata to the Seurat data
  spatial_data$seurat$metadata <- spatial_data$raw_data$metadata
  
  return(spatial_data)
}

# Add filtered data to the spatial_data object
add_filtered_data <- function(spatial_data, mean_expression_threshold = 0.5) {
  # Add filtered data to the spatial_data object
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   threshold: A numeric value representing the minimum mean expression to filter the data (default is 0.05)
  #
  # Returns:
  #   A list containing the spatial_data object with filtered data added
  
  # Determine which rows have a mean expression greater than the threshold
  filter_index <- which(spatial_data$raw_data$simple_statistics$mean > mean_expression_threshold)
  
  # Keep only the filtered metadata
  spatial_data$filtered_data$metadata <- spatial_data$raw_data$metadata
  
  # Keep only the filtered annotation data
  spatial_data$filtered_data$annotate <- spatial_data$raw_data$annotate[filter_index, ]
  
  # Keep only the filtered expression data
  spatial_data$filtered_data$data <- spatial_data$raw_data$data[filter_index, ]
  
  return(spatial_data)
}

# Function to add colfilt (column filtered) data to the spatial_transcriptomic_data object
add_colfilt_data <- function(spatial_data, min_spot_threshold = 0, expression_threshold = 2) {
  # Function to add colfilt (column filtered) data to the spatial_transcriptomic_data object
  # This function filters columns based on the number of genes with expression greater than a given threshold
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   min_spot_threshold: Minimum number of genes with expression greater than the expression threshold required to keep a column
  #   expression_threshold: Expression threshold for filtering columns
  #
  # Returns:
  #   spatial_data: The updated spatial_data list with the colfilt_data added
  
  # Determine which columns have more than (min_genes_expressed - 1) genes with expression greater than the expression_threshold
  col_index <- which(apply(spatial_data$filtered_data$data, 2, function(x) { sum(x > expression_threshold) }) > (min_spot_threshold))
  
  # Keep only the colfilt metadata
  spatial_data$colfilt_data$metadata <- spatial_data$filtered_data$metadata[col_index, ]
  
  # Keep the colfilt annotation data (same as filtered_data$annotate)
  spatial_data$colfilt_data$annotate <- spatial_data$filtered_data$annotate
  
  # Keep only the colfilt expression data
  spatial_data$colfilt_data$data <- spatial_data$filtered_data$data[, col_index]
  
  return(spatial_data)
}

# Function to add range-normalized data to the spatial_transcriptomic_data object
add_range_normalize_data <- function(spatial_data, range = 1500, flatten = 1, threshold = 500) {
  # Function to add range-normalized data to the spatial_transcriptomic_data object
  # The range normalization is applied to the colfilt_data and the result is stored in the range_normalize_data field
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   range: Maximum range for range normalization (default: 1500)
  #   flatten: Value to flatten the expression below the threshold (default: 1)
  #   threshold: Threshold for range normalization (default: 500)
  #
  # Returns:
  #   spatial_data: The updated spatial_data list with the range_normalize_data added
  
  range_normalize <- function(x, range, flatten) {
    wh <- which(x > flatten)
    order <- order(x[wh], decreasing = T)
    out <- x
    len <- length(wh)
    rangecount <- min(range, len)
    maxrange <- range
    minrange <- range - rangecount + 2
    out[out > flatten] <- flatten
    out[wh[order[(rangecount - 1)]]] <- maxrange:minrange 
    out
  }
  
  # Apply range normalization to colfilt data
  spatial_data$range_normalize$data <- apply(
    as.matrix(spatial_data$colfilt_data$data), 1, 
    range_normalize, range = range, flatten = flatten) %>% t
  
  # Keep the same column names as colfilt_data
  colnames(spatial_data$range_normalize$data) <- colnames(spatial_data$colfilt_data$data)
  
  # Keep the same metadata as colfilt_data
  spatial_data$range_normalize$metadata <- spatial_data$colfilt_data$metadata 
  
  # Keep the same annotation data as colfilt_data
  spatial_data$range_normalize$annotate <- spatial_data$colfilt_data$annotate
  
  return(spatial_data)
}


# Function to add clusters data to the spatial_transcriptomic_data object
add_clusters_data <- function(spatial_data, integrated_data, resolution_start = 0.1, resolution_end = 2, resolution_step = 0.1) {
  # Function to add clusters data to the spatial_transcriptomic_data object
  # The clusters are calculated for the given range of resolutions and stored in the clusters field
  #
  # Args:
  #   spatial_data: A list containing the spatial transcriptomic data
  #   integrated_data: An integrated Seurat object
  #   resolution_start: The starting value of the resolution range (default: 0.1)
  #   resolution_end: The ending value of the resolution range (default: 2)
  #   resolution_step: The step value for the resolution range (default: 0.1)
  #
  # Returns:
  #   spatial_data: The updated spatial_data list with the clusters added
  spatial_data$raw_data$metadata %>%
    .[, c(1, 2, 7, 8, 9)] -> cluster_df
  
  resolution_range <- seq(resolution_start, resolution_end, resolution_step)
  
  for(resolution in resolution_range){
    
    tmp_column_name <- paste("cluster_resolution", resolution, sep = "_")
    
    cluster_df <- FindClusters(integrated_data, resolution = resolution) %>%
      .$seurat_clusters %>%
      as.data.frame() %>%
      dplyr::rename({{tmp_column_name}} := ".") %>%
      rownames_to_column(var = "sample_barcode") %>%
      separate("sample_barcode", c("sample", "barcode"), sep = "_") %>%
      right_join(cluster_df, ., by = c("barcode", "sample"))
    rm(tmp_column_name)
  }
  
  # Add clusters to the spatial_data list
  spatial_data$clusters <- cluster_df
  
  return(spatial_data)
}
# maybe is better way to storage cluster info

# Function to evaluate clustering stability for a range of resolutions using silhouette scores.
evaluate_clustering_stability <- function(spatial_data, seurat_object, resolution_start = 0.1, resolution_end = 2, resolution_step = 0.1) {
  # Function to evaluate clustering stability for a range of resolutions using silhouette scores.
  # Silhouette scores measure the quality of clustering by calculating the similarity of objects 
  # within the same cluster compared to objects in other clusters. Scores range from -1 to 1, 
  # with higher values indicating better clustering quality.
  # Returns the spatial_data object with the added stability_results data.frame.
  #
  # Args:
  #  spatial_data: A list containing the spatial transcriptomic data and associated metadata.
  #  resolution_start: A numeric value specifying the start of the resolution range (default = 0.1).
  #  resolution_end: A numeric value specifying the end of the resolution range (default = 2).
  #  resolution_step: A numeric value specifying the step size between resolution values (default = 0.1).
  
  # Generate a sequence of resolution values
  resolution_values <- seq(resolution_start, resolution_end, resolution_step)
  clustering_results <- list()
  
  # Loop through each resolution value and compute cluster assignments
  for (resolution in resolution_values) {
    tmp_object <- FindClusters(seurat_object, resolution = resolution)
    clustering_results[[as.character(resolution)]] <- tmp_object@meta.data$seurat_clusters
  }
  
  # Load the cluster package for silhouette score calculation
  library(cluster)
  
  silhouette_scores <- c()
  num_clusters <- c()
  
  # Loop through each resolution value and calculate silhouette scores and number of clusters
  for (resolution in names(clustering_results)) {
    print(resolution)
    predicted_labels <- clustering_results[[resolution]]
    
    # Convert the factor to numeric
    predicted_labels_numeric <- as.numeric(as.character(predicted_labels))
    
    # Compute the dissimilarity matrix
    dissimilarity_matrix <- dist(seurat_object@reductions$umap@cell.embeddings)
    
    # Use the numeric labels in the silhouette function
    silhouette_result <- silhouette(predicted_labels_numeric, dissimilarity_matrix)
    
    # Calculate the average silhouette score for each resolution
    avg_silhouette_score <- mean(silhouette_result[, 3])
    silhouette_scores <- c(silhouette_scores, avg_silhouette_score)
    
    # Calculate the number of clusters for each resolution
    num_clusters <- c(num_clusters, length(unique(predicted_labels_numeric)))
  }
  
  print("create data frame")
  print(silhouette_scores)
  # Create a data.frame with resolution_values, silhouette_scores, and num_clusters
  stability_results <- data.frame(resolution = resolution_values, 
                                  silhouette_score = silhouette_scores, 
                                  num_clusters = num_clusters)
  
  # Add stability_results to the spatial_data object
  spatial_data$stability_results <- stability_results
  
  return(spatial_data)
}

# This function normalizes gene expression data based on quantiles within each cluster for spatial transcriptomics data.
add_quantile_norm_data <- function(spatial_data, resolution = 0.8, num_cores = 24, data_type = "raw_data"){
  
  # This function normalizes gene expression data based on quantiles within each cluster for spatial transcriptomics data.
  # The normalization is performed on a given resolution level.
  # 'spatial_data': A list containing spatial data with raw_data (gene expression data, metadata, and annotation data) and cluster data.
  # 'resolution': The resolution level at which the data will be clustered. Default is 0.8.
  # 'num_cores': The number of cores to use for parallel processing. Default is 24.
  # 'data_type': The type of data to be used for normalization. Default is "raw_data".
  
  # Function to normalize a given cluster.
  # 'cluster': The cluster id that should be normalized.
  quantile_norm_cluster <- function(cluster){
    
    # Get all barcodes from the cluster to be normalized.
    # The !!sym() function allows us to programmatically use the resolution column name.
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      pull(sample_barcode)
    
    # Extract the data matrix for the specific cluster.
    cluster_data <- spatial_data[[data_type]]$data[, all_barcodes] %>% as.matrix()
    
    # Normalize the data for the cluster.
    norm_data <- normalize.quantiles(cluster_data)
    
    # Keep the column and row names identical to the original data.
    colnames(norm_data) <- colnames(cluster_data)
    rownames(norm_data) <- rownames(cluster_data)
    
    # Return the normalized data.
    return(norm_data)
  } 
  
  # Define the column name based on the resolution parameter.
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  name_data <- paste0("quantile_normalize_", "resolution_", resolution)
  
  # Get a sorted list of unique cluster ids at the given resolution.
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort
  
  # Start the timer to measure the computation time.
  start_time <- Sys.time()
  
  # Apply the quantile_norm_cluster function to each unique cluster in parallel using 'num_cores'.
  results <- mclapply(unique_clusters, quantile_norm_cluster,  mc.cores = num_cores)
  
  # Name each result (i.e., each normalized cluster) by its cluster id.
  results <- setNames(results, paste0("cluster_", unique_clusters))
  
  # Stop the timer and calculate the elapsed time.
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  
  # Print the computation time.
  print(diff_time)
  
  # Combine the normalized clusters back into a single matrix.
  quantile_normalize_data <- do.call(cbind, results) %>% Matrix()
  
  # Add the normalized data and the original metadata and annotation data back into the spatial data.
  spatial_data[[name_data]]$metadata <- spatial_data$raw_data$metadata
  spatial_data[[name_data]]$annotate <- spatial_data$raw_data$annotate
  spatial_data[[name_data]]$data <- quantile_normalize_data
  
  # Return the updated spatial data.
  return(spatial_data)
}



####################################### 
# experimental functions to developed #
#######################################
# Function to calculate the number of barcodes (spots) for each cluster for each sample using a specified resolution.
calculate_spots_per_cluster_matrix <- function(seurat_object, spatial_data, resolution = 0.5) {
  # Function to calculate the number of barcodes (spots) for each cluster for each sample using a specified resolution.
  # Returns a data.frame with samples as rows, clusters as columns, and the number of spots assigned to each
  # cluster in the corresponding cells.
  #
  # Args:
  #  seurat_object: A Seurat object containing the integrated spatial transcriptomic data.
  #  resolution: A numeric value specifying the resolution to use when identifying clusters (default = 0.5).
  
  # Find neighbors and clusters using the specified resolution
  seurat_object <- FindClusters(seurat_object, resolution = resolution)
  
  # Extract the sample, cluster, and barcode information from the Seurat object
  sample_cluster_info <- data.frame(sample = seurat_object@meta.data$orig.ident,
                                    cluster = seurat_object@meta.data$seurat_clusters,
                                    barcode = rownames(seurat_object@meta.data))
  
  # Group the data by sample and cluster and count the number of barcodes in each group
  spots_per_cluster <- sample_cluster_info %>%
    dplyr::group_by(sample, cluster) %>%
    dplyr::summarize(num_barcodes = n()) %>%
    tidyr::spread(cluster, num_barcodes, fill = 0)
  
  # Merge the treatment information from spatial_data$sample_information with the spots_per_cluster data.frame
  spots_per_cluster_with_treatment <-
    right_join({spatial_data$sample_information[, c("sample_ID", "treatment")] %>% rename(sample = sample_ID)},
               spots_per_cluster,
               by = c("sample"))
  
  return(spots_per_cluster_with_treatment)
}






