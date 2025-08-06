compute_data_summary_v2 <- function(spatial_data,
                                    resolution = 0.8,
                                    trim = 0.05,
                                    num_cores = 24,
                                    control_samples,
                                    experiment_samples,
                                    data_type = "raw_data",
                                    min_number_spots = 20,
                                    metrics = c("expression_spot", "sum", "mean", "median", "IQR", "diff_range", "var", "skewness", "kurtosis"),
                                    verbose = FALSE) {
  
  #’ Compute summary statistics for spatial expression data by cluster
  #’
  #’ This function calculates a variety of summary metrics (e.g., sum, mean, median, IQR, variance, skewness, kurtosis) 
  #’ for each cluster in a spatial transcriptomics dataset, comparing control versus experimental groups.
  #’ It operates at a specified clustering resolution, processes only genes/spots with a minimum number of observations,
  #’ and supports parallel execution.
  #’
  #’ @param spatial_data A list containing at least three entries:
  #’   \itemize{
  #’     \item{\code{clusters}}{ — a data.frame of cluster assignments, with column names like \code{cluster_resolution_<value>}.}
  #’     \item{\code{raw_data}}{ (or other \code{data_type} you select) — a list with 
  #’       \code{$data}: a matrix of expression values (rows = peaks/genes, cols = sample_barcode), 
  #’       and \code{$metadata}: a data.frame with columns \code{barcode} and \code{sample}.}
  #’     \item{\code{annotate}}{ — a data.frame with \code{peak_id} and \code{gene_name}.}
  #’   }
  #’ @param resolution Numeric scalar. The clustering resolution to use (e.g. \code{0.8}), 
  #’   which corresponds to the column \code{cluster_resolution_<resolution>} in \code{spatial_data$clusters}.
  #’ @param trim Numeric between 0 and 0.5. The fraction (per tail) to trim when computing trimmed means.
  #’ @param num_cores Integer. Number of CPU cores to use for parallel processing via \code{mclapply}.
  #’ @param control_samples Character vector. Sample IDs designated as “control”.
  #’ @param experiment_samples Character vector. Sample IDs designated as “experiment”.
  #’ @param data_type String. The name of the assay in \code{spatial_data} to summarize (default \code{"raw_data"}).
  #’ @param min_number_spots Integer. Minimum number of spots (observations) required to compute any metric; 
  #’   otherwise returns \code{NA}.
  #’ @param metrics Character vector. Which summary metrics to compute. Options include:
  #’   \code{"expression_spot"}, \code{"sum"}, \code{"mean"}, \code{"median"}, \code{"IQR"}, 
  #’   \code{"diff_range"}, \code{"var"}, \code{"skewness"}, \code{"kurtosis"}.
  #’ @param verbose Logical. If \code{TRUE}, prints progress messages.
  #’
  #’ @return A named list of length equal to the number of clusters. Each element is itself a list with:
  #’   \itemize{
  #’     \item{\code{peak}}{ — vector of peak IDs from \code{spatial_data[[data_type]]$annotate$peak_id}.}
  #’     \item{\code{gene}}{ — vector of gene names from \code{spatial_data[[data_type]]$annotate$gene_name}.}
  #’     \item{\code{control}}{ — a list of metric matrices for control samples.}
  #’     \item{\code{experiment}}{ — a list of metric matrices for experimental samples.}
  #’   }
  #’   If an error occurs in a given cluster, that element will have class \code{"cluster_error"} 
  #’   with an \code{error_message}.
  #’
  #’ @examples
  #’ \dontrun{
  #’   results <- compute_data_summary_v2(
  #’     spatial_data = my_spatial_obj,
  #’     resolution = 0.8,
  #’     control_samples = c("S1","S2"),
  #’     experiment_samples = c("S3","S4"),
  #’     verbose = TRUE
  #’   )
  #’ }
  #’
  #’ @export
  
  # Construct the column name for the specified clustering resolution
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- unique(spatial_data$clusters[[resolution_column]]) %>% sort()
  
  if (verbose) {
    cat(">>> START compute_data_summary_v2()\n")
    cat(">>> Using resolution column: ", resolution_column, "\n")
    cat(">>> Identified clusters: ", paste(unique_clusters, collapse = ", "), "\n")
    cat(">>> Selected metrics: ", paste(metrics, collapse = ", "), "\n\n")
  }
  
  # Internal helper to compute summary for a single cluster
  compute_data_summary_cluster <- function(cluster) {
    if (verbose) cat(">>> Processing cluster: ", cluster, "\n")
    
    # Join cluster assignments with expression metadata for selected data type
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, ., by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(sample_barcode = paste(sample, barcode, sep = "_")) %>%
      mutate(group = ifelse(sample %in% control_samples, "control", "experiment")) %>%
      select(group, sample_barcode)
    
    barcodes_list <- split(all_barcodes$sample_barcode, all_barcodes$group)
    control_barcodes <- barcodes_list[["control"]]
    experiment_barcodes <- barcodes_list[["experiment"]]
    
    # Extract expression matrices for control and experiment groups
    control_expression_spot <- sapply(control_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, grep(paste0(sample_id, "_"), control_barcodes), drop = FALSE]
    }, simplify = FALSE)
    experiment_expression_spot <- sapply(experiment_samples, function(sample_id) {
      spatial_data[[data_type]]$data[, grep(paste0(sample_id, "_"), experiment_barcodes), drop = FALSE]
    }, simplify = FALSE)
    
    control_list <- list()
    experiment_list <- list()
    
    # Compute each requested metric
    for (metric in metrics) {
      if (metric == "expression_spot") {
        control_list$expression_spot <- control_expression_spot
        experiment_list$expression_spot <- experiment_expression_spot
        next
      }
      
      # Compute metric across spots for each sample
      control_list[[metric]] <- sapply(control_expression_spot, function(x) {
        apply(x, 1, function(y) {
          if (length(y) < min_number_spots) {
            NA
          } else {
            switch(metric,
                   sum        = sum(y, na.rm = TRUE),
                   mean       = mean(y, trim = trim, na.rm = TRUE),
                   median     = median(y, na.rm = TRUE),
                   IQR        = IQR(y, na.rm = TRUE),
                   diff_range = diff(range(y, na.rm = TRUE)),
                   var        = var(y, na.rm = TRUE),
                   skewness   = skewness(y),
                   kurtosis   = kurtosis(y),
                   NA)
          }
        })
      })
      
      experiment_list[[metric]] <- sapply(experiment_expression_spot, function(x) {
        apply(x, 1, function(y) {
          if (length(y) < min_number_spots) {
            NA
          } else {
            switch(metric,
                   sum        = sum(y, na.rm = TRUE),
                   mean       = mean(y, trim = trim, na.rm = TRUE),
                   median     = median(y, na.rm = TRUE),
                   IQR        = IQR(y, na.rm = TRUE),
                   diff_range = diff(range(y, na.rm = TRUE)),
                   var        = var(y, na.rm = TRUE),
                   skewness   = skewness(y),
                   kurtosis   = kurtosis(y),
                   NA)
          }
        })
      })
    }
    
    # Return summary results for this cluster
    list(
      peak    = spatial_data[[data_type]]$annotate$peak_id,
      gene    = spatial_data[[data_type]]$annotate$gene_name,
      control = control_list,
      experiment = experiment_list
    )
  }
  
  start_time <- Sys.time()
  
  # Parallel computation across clusters
  results <- mclapply(unique_clusters, function(cl) {
    if (verbose) cat("\n>>> STARTING CLUSTER: ", cl, "\n")
    tryCatch({
      compute_data_summary_cluster(cl)
    }, error = function(e) {
      if (verbose) message("!! ERROR in cluster ", cl, ": ", e$message)
      structure(list(error_message = e$message), class = "cluster_error")
    })
  }, mc.cores = num_cores)
  names(results) <- paste0("cluster_", unique_clusters)
  
  # Report clusters with errors
  error_clusters <- sapply(results, function(x) inherits(x, "cluster_error"))
  if (any(error_clusters)) {
    message("\n!!! WARNING: errors in clusters: ", paste(names(results)[error_clusters], collapse = ", "))
  }
  
  end_time <- Sys.time()
  if (verbose) cat("\n>>> Total computation time: ", round(end_time - start_time, 2), "seconds\n")
  
  return(results)
}

perform_statistical_tests_v2 <- function(
  spatial_data,
  summary_data,
  metric = "mean",
  resolution = 0.8,
  num_cores = 24,
  mean_threshold = 0,
  control_samples,
  experiment_samples,
  quantile_normalization = FALSE,
  verbose = FALSE
) {
  #' Perform statistical testing for each cluster in spatial transcriptomic data
  #'
  #' This function performs three statistical tests (Student's t-test, Wilcoxon rank-sum test,
  #' and Kolmogorov–Smirnov test) to compare control and experimental samples in each cluster
  #' for a given summary metric. It returns p-values and descriptive statistics for each peak/gene.
  #'
  #' @param spatial_data A list with clustering assignments under \code{$clusters}. The element
  #'   \code{$clusters} must contain a column \code{cluster_resolution_<value>} for each resolution.
  #' @param summary_data A named list of cluster-wise summary statistics,
  #'   typically returned by \code{compute_data_summary_v2()}.
  #' @param metric Character string specifying which metric to test (e.g., \code{"mean"}, \code{"median"}).
  #' @param resolution Numeric. Clustering resolution used to determine cluster column name.
  #' @param num_cores Integer. Number of CPU cores to use for parallelization via \code{mclapply()}.
  #' @param mean_threshold Numeric. If both control and experiment means are below this value,
  #'   statistical tests are skipped and result is marked as \code{"low_mean"}.
  #' @param control_samples Character vector of sample IDs representing the control group.
  #' @param experiment_samples Character vector of sample IDs representing the experimental group.
  #' @param quantile_normalization Logical. If \code{TRUE}, quantile normalization is applied before testing.
  #' @param verbose Logical. If \code{TRUE}, prints progress information for each cluster.
  #'
  #' @return A named list of cluster results. For each cluster:
  #'   \itemize{
  #'     \item If successful, returns a list including test p-values, log2 fold-change, and descriptive statistics.
  #'     \item If error occurs, returns an object of class \code{"cluster_error"} with an \code{error_message}.
  #'   }
  #'
  #' @export
  
  # Define the name of the clustering column
  resolution_column <- paste0("cluster_resolution_", resolution)
  unique_clusters <- sort(unique(spatial_data$clusters[[resolution_column]]))
  cluster_names <- paste0("cluster_", unique_clusters)
  
  if (verbose) {
    cat(">>> START perform_statistical_tests_v2()\n")
    cat(">>> Metric:", metric, "\n")
    cat(">>> Clusters:", paste(cluster_names, collapse = ", "), "\n\n")
  }
  
  # Internal function to process a single cluster
  perform_statistical_tests_cluster <- function(cluster_name) {
    if (verbose) cat(">>> Testing cluster:", cluster_name, "\n")
    
    tryCatch({
      # Retrieve the summary data for the current cluster
      cluster_data <- summary_data[[cluster_name]]
      
      # Extract expression matrix for control and experimental samples
      control_expression <- cluster_data$control[[metric]][, control_samples, drop = FALSE]
      experiment_expression <- cluster_data$experiment[[metric]][, experiment_samples, drop = FALSE]
      
      # Skip clusters with no data in both groups
      if (all(is.na(control_expression)) && all(is.na(experiment_expression))) {
        stop("All values are NA in both groups")
      }
      
      # Apply quantile normalization if requested
      if (quantile_normalization) {
        expr_raw <- cbind(control_expression, experiment_expression)
        expr_norm <- preprocessCore::normalize.quantiles(expr_raw)
        colnames(expr_norm) <- colnames(expr_raw)
        rownames(expr_norm) <- rownames(expr_raw)
        control_expression <- expr_norm[, control_samples, drop = FALSE]
        experiment_expression <- expr_norm[, experiment_samples, drop = FALSE]
      }
      
      # Initialize result vectors
      n <- nrow(control_expression)
      t_test <- numeric(n)
      wilcoxon_test <- numeric(n)
      ks_test <- numeric(n)
      condition <- character(n)
      
      # Loop over each peak/gene
      for (i in seq_len(n)) {
        ctrl_vals <- control_expression[i, ]
        exp_vals <- experiment_expression[i, ]
        
        if (all(is.na(ctrl_vals)) || all(is.na(exp_vals))) {
          # All values missing — skip
          condition[i] <- "all_NA"
          t_test[i] <- 1
          wilcoxon_test[i] <- 1
          ks_test[i] <- 1
        } else if (
          mean(ctrl_vals, na.rm = TRUE) < mean_threshold ||
          mean(exp_vals, na.rm = TRUE) < mean_threshold
        ) {
          # Means too low — skip
          condition[i] <- "low_mean"
          t_test[i] <- 1
          wilcoxon_test[i] <- 1
          ks_test[i] <- 1
        } else {
          # Perform tests
          condition[i] <- "yes_mean"
          
          t_res <- tryCatch(t.test(ctrl_vals, exp_vals, var.equal = TRUE), error = function(e) list(p.value = 1))
          wilc_res <- tryCatch(wilcox.test(ctrl_vals, exp_vals), error = function(e) list(p.value = 1))
          ks_res <- tryCatch(ks.test(ctrl_vals, exp_vals), error = function(e) list(p.value = 1))
          
          t_test[i] <- t_res$p.value
          wilcoxon_test[i] <- wilc_res$p.value
          ks_test[i] <- ks_res$p.value
        }
      }
      
      # Compute summary statistics
      control_mean <- rowMeans(control_expression, na.rm = TRUE)
      experiment_mean <- rowMeans(experiment_expression, na.rm = TRUE)
      control_sd <- apply(control_expression, 1, sd, na.rm = TRUE)
      experiment_sd <- apply(experiment_expression, 1, sd, na.rm = TRUE)
      control_se <- control_sd / sqrt(ncol(control_expression))
      experiment_se <- experiment_sd / sqrt(ncol(experiment_expression))
      log2ratio <- log2(experiment_mean / control_mean)
      
      # Assemble result
      stats_list <- list(
        t_test = t_test,
        wilcoxon_test = wilcoxon_test,
        ks_test = ks_test,
        log2ratio = unname(log2ratio),
        condition = condition,
        control_mean = control_mean,
        experiment_mean = experiment_mean,
        control_sd = control_sd,
        experiment_sd = experiment_sd,
        control_se = control_se,
        experiment_se = experiment_se
      )
      
      # Add optional skewness and kurtosis
      if (!is.null(cluster_data$control$skewness) && !is.null(cluster_data$experiment$skewness)) {
        stats_list$control_skewness <- rowMeans(cluster_data$control$skewness, na.rm = TRUE)
        stats_list$experiment_skewness <- rowMeans(cluster_data$experiment$skewness, na.rm = TRUE)
        stats_list$control_kurtosis <- rowMeans(cluster_data$control$kurtosis, na.rm = TRUE)
        stats_list$experiment_kurtosis <- rowMeans(cluster_data$experiment$kurtosis, na.rm = TRUE)
      }
      
      # Store results inside cluster_data object
      cluster_data$statistics[[metric]] <- stats_list
      return(cluster_data)
      
    }, error = function(e) {
      # If anything fails, return an error object instead of result
      message("!! ERROR in ", cluster_name, ": ", e$message)
      structure(list(error_message = e$message, cluster = cluster_name), class = "cluster_error")
    })
  }
  
  # Run in parallel across all clusters
  start_time <- Sys.time()
  results_list <- parallel::mclapply(cluster_names, perform_statistical_tests_cluster, mc.cores = num_cores)
  names(results_list) <- cluster_names
  end_time <- Sys.time()
  
  if (verbose) {
    cat(">>> Execution time:", round(end_time - start_time, 2), "seconds\n")
  }
  
  return(results_list)
}



summarize_and_test_v2 <- function(spatial_data,
                                  data_params_df,
                                  control_samples,
                                  experiment_samples,
                                  summary_metrics = c("mean", "median"),
                                  statistic_metrics = c("mean", "median"),
                                  mean_threshold = 0,
                                  trim = 0.05,
                                  num_cores = 20,
                                  verbose = TRUE) {
  
  #' Summarize data at multiple resolutions and perform statistical tests
  #'
  #' This function iterates over a data-parameters data.frame, computes summary statistics
  #' for each combination of data type and clustering resolution, and then runs statistical
  #' tests on the requested metrics. Progress is reported with a text progress bar.
  #'
  #' @param spatial_data A list-like object containing spatial transcriptomics data; 
  #'   passed on to \code{compute_data_summary_v2()} and \code{perform_statistical_tests_v2()}.
  #' @param data_params_df A data.frame with one row per analysis step. Columns should be in this order:
  #'   \itemize{
  #'     \item{\strong{data_type} (chr)}{ — assay name to summarize (e.g. “raw_data”).}
  #'     \item{\strong{resolution} (num)}{ — clustering resolution (e.g. 0.8).}
  #'     \item{\strong{data_type_name} (chr)}{ — label under which to store results.}
  #'     \item{\strong{quantile_normalization} (logical)}{ — whether to apply quantile normalization before testing.}
  #'   }
  #' @param control_samples   Character vector of control sample IDs.
  #' @param experiment_samples Character vector of experimental sample IDs.
  #' @param summary_metrics   Character vector of metrics to compute when summarizing 
  #'   (passed as \code{metrics} to \code{compute_data_summary_v2()}; default \code{c("mean","median")}).
  #' @param statistic_metrics Character vector of metrics to test (passed as \code{metric} to 
  #'   \code{perform_statistical_tests_v2()}; default \code{c("mean","median")}).
  #' @param mean_threshold    Numeric. Minimum mean value below which tests are skipped.
  #' @param trim              Numeric (0–0.5). Fraction to trim in trimmed-mean calculations.
  #' @param num_cores         Integer. Number of cores for parallel computation.
  #' @param verbose           Logical. If \code{TRUE}, prints progress and status messages.
  #'
  #' @details
  #' The function does the following for each row of \code{data_params_df}:
  #' \enumerate{
  #'   \item Calls \code{compute_data_summary_v2()} to get control vs. experiment summary stats.
  #'   \item Stores the summary under \code{result_summary_statistics[[data_type_name]][["resolution_<resolution>"]]}.
  #'   \item Iterates over \code{statistic_metrics}, calling \code{perform_statistical_tests_v2()} to add p-values, fold-changes, etc.
  #'   \item Removes any cluster that failed during testing (i.e., returned an object of class \code{"cluster_error"}).
  #' }
  #' A text progress bar tracks overall progress.
  #'
  #' @return A nested list:
  #'   \code{result_summary_statistics[[data_type_name]][["resolution_<res>"]]} 
  #'   contains the output of \code{perform_statistical_tests_v2()} for that combination.
  #'
  #' @seealso
  #' \code{\link{compute_data_summary_v2}} for summary computation  
  #' \code{\link{perform_statistical_tests_v2}} for per-metric hypothesis testing
  #'
  #' @examples
  #' \dontrun{
  #' params <- data.frame(
  #'   data_type = c("raw_data","normalized_data"),
  #'   resolution = c(0.8, 1.0),
  #'   data_type_name = c("raw","norm"),
  #'   quantile_normalization = c(FALSE, TRUE),
  #'   stringsAsFactors = FALSE
  #' )
  #' result <- summarize_and_test_v2(
  #'   spatial_data        = my_spatial_obj,
  #'   data_params_df      = params,
  #'   control_samples     = c("C1","C2"),
  #'   experiment_samples  = c("E1","E2"),
  #'   summary_metrics     = c("mean","median"),
  #'   statistic_metrics   = c("mean","median"),
  #'   mean_threshold      = 0.1,
  #'   trim                = 0.05,
  #'   num_cores           = 4,
  #'   verbose             = TRUE
  #' )
  #' }
  #'
  #' @export
  
  result_summary_statistics <- list()
  n_steps <- nrow(data_params_df)
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  
  for (i in seq_len(n_steps)) {
    # Extract parameters
    data_type  <- data_params_df[i, 1, drop = TRUE]
    resolution <- data_params_df[i, 2, drop = TRUE]
    data_type_name <- data_params_df[i, 3, drop = TRUE]
    quantile_normalization <- data_params_df[i, 4, drop = TRUE]
    
    if (verbose) {
      cat("\n=============================\n")
      cat(">>> ROW", i, "/", n_steps, ": ", data_type, "@", resolution, "\n")
      cat("=============================\n")
    }
    
    # Step 1: Compute summary statistics
    summary_data <- compute_data_summary_v2(
      spatial_data        = spatial_data,
      resolution          = resolution,
      trim                = trim,
      num_cores           = num_cores,
      control_samples     = control_samples,
      experiment_samples  = experiment_samples,
      data_type           = data_type,
      metrics             = summary_metrics,
      verbose             = verbose
    )
    
    # Store initial result
    result_summary_statistics[[data_type_name]][[paste0("resolution_", resolution)]] <- summary_data
    
    # Step 2: Test each requested metric
    for (metric in statistic_metrics) {
      if (verbose) cat("\n>>> Testing metric:", metric, "\n")
      
      summary_data <- perform_statistical_tests_v2(
        spatial_data           = spatial_data,
        summary_data           = summary_data,
        metric                 = metric,
        resolution             = resolution,
        num_cores              = num_cores,
        mean_threshold         = mean_threshold,
        control_samples        = control_samples,
        experiment_samples     = experiment_samples,
        quantile_normalization = quantile_normalization,
        verbose                = verbose
      )
      
      # Identify and remove cluster-level errors
      error_mask <- sapply(summary_data, function(x) inherits(x, "cluster_error"))
      error_names <- names(summary_data)[error_mask]
      summary_data <- summary_data[!error_mask]
      
      if (verbose && length(error_names) > 0) {
        cat(">>> Removed", length(error_names), "cluster(s) due to errors in metric:", metric, "\n")
        cat(">>> Removed clusters:", paste(error_names, collapse = ", "), "\n")
      }
    }
    
    # Save final result
    result_summary_statistics[[data_type_name]][[paste0("resolution_", resolution)]] <- summary_data
    
    # Progress bar update
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(result_summary_statistics)
}






