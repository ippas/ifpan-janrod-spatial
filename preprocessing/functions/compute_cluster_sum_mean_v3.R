
compute_cluster_sum_mean_v3 <- function(
  spatial_data,
  resolution,
  samples,
  data_type = "raw_data",
  min_number_spots = 15,
  num_cores = 24,
  verbose = TRUE
) {
  
  # ============================================================
  # 1. Basic input checks
  # ============================================================
  
  if (!data_type %in% names(spatial_data)) {
    stop(
      "data_type '", data_type,
      "' was not found in spatial_data."
    )
  }
  
  if (is.null(spatial_data$clusters)) {
    stop("spatial_data$clusters is missing.")
  }
  
  if (is.null(spatial_data[[data_type]]$data)) {
    stop("spatial_data[[data_type]]$data is missing.")
  }
  
  if (is.null(spatial_data[[data_type]]$metadata)) {
    stop("spatial_data[[data_type]]$metadata is missing.")
  }
  
  if (length(samples) == 0) {
    stop("Provide at least one sample in `samples`.")
  }
  
  samples <- unique(as.character(samples))
  
  resolution_column <- paste0(
    "cluster_resolution_",
    resolution
  )
  
  if (!resolution_column %in% colnames(spatial_data$clusters)) {
    stop(
      "Resolution column not found in spatial_data$clusters: ",
      resolution_column
    )
  }
  
  # ============================================================
  # 2. Extract core objects
  # ============================================================
  
  expression_data <- spatial_data[[data_type]]$data
  expression_metadata <- spatial_data[[data_type]]$metadata
  clusters <- spatial_data$clusters
  
  if (!is.null(spatial_data[[data_type]]$annotate)) {
    annotate <- spatial_data[[data_type]]$annotate
  } else if (!is.null(spatial_data$annotate)) {
    annotate <- spatial_data$annotate
  } else {
    stop(
      "Annotation table not found. Expected either ",
      "spatial_data[[data_type]]$annotate or spatial_data$annotate."
    )
  }
  
  # ============================================================
  # 3. Validate required columns
  # ============================================================
  
  required_cluster_columns <- c(
    "barcode",
    "sample",
    resolution_column
  )
  
  required_metadata_columns <- c(
    "barcode",
    "sample"
  )
  
  missing_cluster_columns <- setdiff(
    required_cluster_columns,
    colnames(clusters)
  )
  
  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(expression_metadata)
  )
  
  if (length(missing_cluster_columns) > 0) {
    stop(
      "Missing columns in spatial_data$clusters: ",
      paste(missing_cluster_columns, collapse = ", ")
    )
  }
  
  if (length(missing_metadata_columns) > 0) {
    stop(
      "Missing columns in spatial_data[[data_type]]$metadata: ",
      paste(missing_metadata_columns, collapse = ", ")
    )
  }
  
  if (is.null(colnames(expression_data))) {
    stop(
      "Expression matrix must contain column names ",
      "in sample_barcode format."
    )
  }
  
  if (nrow(annotate) != nrow(expression_data)) {
    stop(
      "Number of rows in annotation table does not match ",
      "number of rows in expression matrix."
    )
  }
  
  # ============================================================
  # 4. Prepare feature annotation
  # ============================================================
  
  peak_ids <- if ("peak_id" %in% colnames(annotate)) {
    as.character(annotate$peak_id)
  } else {
    rownames(expression_data)
  }
  
  gene_names <- if ("gene_name" %in% colnames(annotate)) {
    as.character(annotate$gene_name)
  } else {
    rownames(expression_data)
  }
  
  if (is.null(peak_ids)) {
    peak_ids <- paste0(
      "feature_",
      seq_len(nrow(expression_data))
    )
  }
  
  if (is.null(gene_names)) {
    gene_names <- peak_ids
  }
  
  # ============================================================
  # 5. Construct unique sample_barcode identifiers
  # ============================================================
  
  expression_metadata <- expression_metadata %>%
    mutate(
      sample = as.character(sample),
      barcode = as.character(barcode),
      sample_barcode = paste(sample, barcode, sep = "_")
    )
  
  clusters <- clusters %>%
    mutate(
      sample = as.character(sample),
      barcode = as.character(barcode),
      sample_barcode = paste(sample, barcode, sep = "_")
    )
  
  if (anyDuplicated(expression_metadata$sample_barcode)) {
    stop(
      "Duplicated sample_barcode values found in ",
      "spatial_data[[data_type]]$metadata."
    )
  }
  
  if (anyDuplicated(clusters$sample_barcode)) {
    stop(
      "Duplicated sample_barcode values found in spatial_data$clusters."
    )
  }
  
  missing_expression_columns <- setdiff(
    expression_metadata$sample_barcode,
    colnames(expression_data)
  )
  
  if (length(missing_expression_columns) > 0) {
    stop(
      "Some metadata sample_barcodes are missing from expression matrix columns. ",
      "Example: ",
      missing_expression_columns[1]
    )
  }
  
  available_samples <- expression_metadata %>%
    pull(sample) %>%
    unique()
  
  missing_samples <- setdiff(
    samples,
    available_samples
  )
  
  if (length(missing_samples) > 0) {
    stop(
      "These requested samples are not present in expression metadata: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  # ============================================================
  # 6. Keep selected samples and valid cluster spots
  # ============================================================
  
  spot_info <- clusters %>%
    filter(
      sample %in% samples,
      !is.na(.data[[resolution_column]])
    ) %>%
    select(
      sample,
      sample_barcode,
      cluster = all_of(resolution_column)
    ) %>%
    filter(
      sample_barcode %in% expression_metadata$sample_barcode
    )
  
  if (nrow(spot_info) == 0) {
    stop(
      "No spots remained after filtering selected samples and clusters."
    )
  }
  
  unique_clusters <- spot_info %>%
    pull(cluster) %>%
    unique() %>%
    sort()
  
  if (verbose) {
    cat(">>> START compute_cluster_sum_mean_v3()\n")
    cat(">>> Resolution column: ", resolution_column, "\n", sep = "")
    cat(">>> Number of selected samples: ", length(samples), "\n", sep = "")
    cat(">>> Samples: ", paste(samples, collapse = ", "), "\n", sep = "")
    cat(">>> Number of clusters: ", length(unique_clusters), "\n", sep = "")
    cat(">>> Clusters: ", paste(unique_clusters, collapse = ", "), "\n", sep = "")
    cat(">>> Minimum spots per sample × cluster: ",
        min_number_spots, "\n\n", sep = "")
  }
  
  # ============================================================
  # 7. Helper for sparse / dense matrices
  # ============================================================
  
  row_sum_safe <- function(x) {
    
    if (inherits(x, "Matrix")) {
      return(as.numeric(Matrix::rowSums(x)))
    }
    
    return(as.numeric(rowSums(x)))
  }
  
  # ============================================================
  # 8. Compute sum and mean for one cluster
  # ============================================================
  
  compute_one_cluster <- function(cluster_id) {
    
    cluster_spots <- spot_info %>%
      filter(cluster == cluster_id)
    
    n_spots_per_sample <- sapply(
      samples,
      function(sample_id) {
        sum(cluster_spots$sample == sample_id)
      }
    )
    
    sample_qc <- data.frame(
      cluster = rep(cluster_id, length(samples)),
      sample = samples,
      n_spots = as.integer(n_spots_per_sample),
      pass_min_spots = n_spots_per_sample >= min_number_spots,
      stringsAsFactors = FALSE
    )
    
    sum_matrix <- matrix(
      NA_real_,
      nrow = nrow(expression_data),
      ncol = length(samples),
      dimnames = list(peak_ids, samples)
    )
    
    mean_matrix <- matrix(
      NA_real_,
      nrow = nrow(expression_data),
      ncol = length(samples),
      dimnames = list(peak_ids, samples)
    )
    
    for (sample_id in samples) {
      
      current_barcodes <- cluster_spots %>%
        filter(sample == sample_id) %>%
        pull(sample_barcode)
      
      current_n_spots <- length(current_barcodes)
      
      if (current_n_spots < min_number_spots) {
        next
      }
      
      current_columns <- match(
        current_barcodes,
        colnames(expression_data)
      )
      
      current_expression <- expression_data[
        ,
        current_columns,
        drop = FALSE
      ]
      
      current_sum <- row_sum_safe(current_expression)
      
      sum_matrix[, sample_id] <- current_sum
      
      mean_matrix[, sample_id] <- current_sum / current_n_spots
    }
    
    list(
      peak = peak_ids,
      gene = gene_names,
      sample_qc = sample_qc,
      sum = sum_matrix,
      mean = mean_matrix
    )
  }
  
  # ============================================================
  # 9. Run across clusters
  # ============================================================
  
  start_time <- Sys.time()
  
  results <- parallel::mclapply(
    unique_clusters,
    function(cluster_id) {
      
      if (verbose) {
        cat(
          ">>> Processing cluster: ",
          cluster_id,
          "\n",
          sep = ""
        )
      }
      
      tryCatch(
        compute_one_cluster(cluster_id),
        error = function(e) {
          structure(
            list(
              cluster = cluster_id,
              error_message = e$message
            ),
            class = "cluster_error"
          )
        }
      )
    },
    mc.cores = min(num_cores, length(unique_clusters))
  )
  
  names(results) <- paste0(
    "cluster_",
    unique_clusters
  )
  
  error_clusters <- vapply(
    results,
    function(x) inherits(x, "cluster_error"),
    logical(1)
  )
  
  if (any(error_clusters)) {
    warning(
      "Errors occurred in clusters: ",
      paste(
        names(results)[error_clusters],
        collapse = ", "
      )
    )
  }
  
  end_time <- Sys.time()
  
  if (verbose) {
    runtime_seconds <- round(
      as.numeric(
        difftime(
          end_time,
          start_time,
          units = "secs"
        )
      ),
      1
    )
    
    cat(
      "\n>>> Total runtime: ",
      runtime_seconds,
      " seconds\n",
      sep = ""
    )
  }
  
  # ============================================================
  # 10. Add metadata to whole result object
  # ============================================================
  
  attr(results, "parameters") <- list(
    resolution = resolution,
    data_type = data_type,
    samples = samples,
    min_number_spots = min_number_spots,
    metrics = c("sum", "mean")
  )
  
  class(results) <- c(
    "cluster_sum_mean_summary",
    class(results)
  )
  
  return(results)
}