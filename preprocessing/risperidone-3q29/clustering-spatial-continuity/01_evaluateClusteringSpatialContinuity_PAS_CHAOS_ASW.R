#!/usr/bin/env Rscript

# ==============================================================================
# 01_evaluateClusteringSpatialContinuity_PAS_CHAOS_ASW_v4.R
#
# Dataset:
#   risperidone-3q29
#
# Metrics calculated separately for each sample and clustering resolution:
#   1. transcriptomic_ASW - PCA space, dims 1:30; higher = better
#   2. spatial_ASW        - physical spot coordinates; higher = better
#   3. CHAOS              - physical spot coordinates; lower = better
#   4. PAS                - physical spot coordinates; lower = better
#
# PAS definition:
#   k = 10 spatial nearest neighbours;
#   abnormal spot = >= 6/10 neighbours assigned to another cluster.
#
# CHAOS definition:
#   nearest same-cluster spatial-neighbour distance, normalized by median
#   nearest-neighbour spot pitch within each tissue section.
#
# Output table:
#   resolution | method | sample_1 ... sample_26 | median | mean | SD
#
# Output directory:
#   results/risperidone-3q29/clustering-spatial-continuity/
#
# Output files:
#   ris3q29_clusteringSpatialContinuity_PAS_CHAOS_ASW.tsv
#   ris3q29_clusteringSpatialContinuity_PAS_CHAOS_ASW.xlsx
#   ris3q29_clusteringSpatialContinuity_PAS_CHAOS_ASW_allResolutions.pdf
#   ris3q29_clusteringSpatialContinuity_PAS_CHAOS_ASW_resolutionStep010.pdf
#
# PDFs:
#   Each PDF contains four pages:
#   page 1 - transcriptomic ASW
#   page 2 - spatial ASW
#   page 3 - CHAOS
#   page 4 - PAS
#
#   PDF 1 - all tested clustering resolutions.
#   PDF 2 - only resolutions in steps of 0.10.
#
# Plotting:
#   black points only; each point is the mean across samples;
#   no connecting lines and no SD/error bars are shown on plots.
#   Median, mean and SD are retained in TSV/XLSX outputs.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ifpan-janrod-spatial"

integrated_rdata_file <- file.path(
  project_root,
  "results",
  "risperidone-3q29",
  "ris3q29_integrate.RData"
)

spatial_rdata_file <- file.path(
  project_root,
  "results",
  "risperidone-3q29",
  "ris3q29_st_data.RData"
)

output_dir <- file.path(
  project_root,
  "results",
  "risperidone-3q29",
  "clustering-spatial-continuity"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_prefix <- "ris3q29_clusteringSpatialContinuity_PAS_CHAOS_ASW"
output_tsv <- file.path(output_dir, paste0(output_prefix, ".tsv"))
output_xlsx <- file.path(output_dir, paste0(output_prefix, ".xlsx"))
output_pdf_all_resolutions <- file.path(
  output_dir,
  paste0(output_prefix, "_allResolutions.pdf")
)

output_pdf_resolution_step_010 <- file.path(
  output_dir,
  paste0(output_prefix, "_resolutionStep010.pdf")
)

transcriptomic_reduction <- "pca"
transcriptomic_dims <- 1:30
transcriptomic_asw_max_spots_per_sample <- 4000L
transcriptomic_asw_minimum_spots_per_cluster <- 20L

spatial_asw_max_spots_per_sample <- 3000L
spatial_asw_minimum_spots_per_cluster <- 10L

pas_k <- 10L
pas_minimum_different_neighbours <- 6L

normalize_chaos_by_spot_pitch <- TRUE
random_seed <- 7L

method_order <- c(
  "transcriptomic_ASW",
  "spatial_ASW",
  "CHAOS",
  "PAS"
)

# ==============================================================================
# 2. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "cluster",
  "RANN",
  "ggplot2",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# ==============================================================================
# 3. General utilities
# ==============================================================================

workflow_message <- function(...) {
  cat(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )
  flush.console()
  invisible(NULL)
}

format_elapsed_time <- function(seconds) {
  seconds <- max(0, as.numeric(seconds))
  sprintf(
    "%02d:%02d:%02d",
    as.integer(floor(seconds / 3600)),
    as.integer(floor((seconds %% 3600) / 60)),
    as.integer(floor(seconds %% 60))
  )
}

load_expected_object <- function(rdata_file, expected_object_name) {
  if (!file.exists(rdata_file)) {
    stop("Input RData file does not exist:\n", rdata_file, call. = FALSE)
  }
  
  load_environment <- new.env(parent = globalenv())
  loaded_object_names <- load(rdata_file, envir = load_environment)
  
  if (!expected_object_name %in% loaded_object_names) {
    stop(
      "Expected object `", expected_object_name, "` was not found in:\n",
      rdata_file,
      "\nObjects found: ",
      paste(loaded_object_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  get(expected_object_name, envir = load_environment, inherits = FALSE)
}

write_tsv <- function(data, filename) {
  utils::write.table(
    x = data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(filename)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

# ==============================================================================
# 4. Common spot selection for exact ASW
# ==============================================================================

select_common_validation_indices <- function(
  cluster_table,
  cluster_columns,
  max_spots,
  minimum_spots_per_cluster,
  seed
) {
  n_spots <- nrow(cluster_table)
  
  if (n_spots < 2L) {
    stop("At least two spots are required.", call. = FALSE)
  }
  
  if (n_spots <= max_spots) {
    return(seq_len(n_spots))
  }
  
  set.seed(as.integer(seed))
  
  all_indices <- seq_len(n_spots)
  mandatory_indices <- integer(0)
  
  for (cluster_column in cluster_columns) {
    cluster_labels <- as.character(cluster_table[[cluster_column]])
    cluster_levels <- unique(cluster_labels)
    
    for (cluster_level in cluster_levels) {
      available_indices <- which(cluster_labels == cluster_level)
      number_to_select <- min(
        length(available_indices),
        as.integer(minimum_spots_per_cluster)
      )
      
      if (number_to_select > 0L) {
        selected_indices <- sample(
          available_indices,
          size = number_to_select,
          replace = FALSE
        )
        mandatory_indices <- union(mandatory_indices, selected_indices)
      }
    }
  }
  
  target_number <- min(n_spots, as.integer(max_spots))
  
  if (length(mandatory_indices) < target_number) {
    remaining_indices <- setdiff(all_indices, mandatory_indices)
    number_to_add <- min(
      target_number - length(mandatory_indices),
      length(remaining_indices)
    )
    
    if (number_to_add > 0L) {
      additional_indices <- sample(
        remaining_indices,
        size = number_to_add,
        replace = FALSE
      )
      selected_indices <- c(mandatory_indices, additional_indices)
    } else {
      selected_indices <- mandatory_indices
    }
  } else {
    selected_indices <- sample(
      mandatory_indices,
      size = target_number,
      replace = FALSE
    )
  }
  
  sort(unique(selected_indices))
}

# ==============================================================================
# 5. ASW
# ==============================================================================

calculate_asw_from_distance <- function(distance_object, cluster_labels) {
  cluster_labels <- as.character(cluster_labels)
  n_observations <- attr(distance_object, "Size")
  
  if (is.null(n_observations)) {
    stop("`distance_object` is not a valid dist object.", call. = FALSE)
  }
  
  if (n_observations != length(cluster_labels)) {
    stop(
      "Distance matrix and cluster-label lengths differ.",
      call. = FALSE
    )
  }
  
  if (n_observations < 3L) return(NA_real_)
  
  cluster_factor <- factor(cluster_labels)
  
  if (nlevels(cluster_factor) < 2L) return(NA_real_)
  if (nlevels(cluster_factor) >= n_observations) return(NA_real_)
  
  silhouette_object <- cluster::silhouette(
    x = as.integer(cluster_factor),
    dist = distance_object
  )
  
  silhouette_values <- as.numeric(silhouette_object[, "sil_width"])
  safe_mean(silhouette_values)
}

# ==============================================================================
# 6. Spatial spot pitch
# ==============================================================================

estimate_spatial_spot_pitch <- function(coordinate_matrix) {
  coordinate_matrix <- as.matrix(coordinate_matrix)
  
  if (nrow(coordinate_matrix) < 2L) {
    stop(
      "At least two spots are required to estimate spatial spot pitch.",
      call. = FALSE
    )
  }
  
  nearest_neighbour_result <- RANN::nn2(
    data = coordinate_matrix,
    query = coordinate_matrix,
    k = 2L
  )
  
  nearest_neighbour_distances <- nearest_neighbour_result$nn.dists[, 2L]
  nearest_neighbour_distances <- nearest_neighbour_distances[
    is.finite(nearest_neighbour_distances) & nearest_neighbour_distances > 0
  ]
  
  if (length(nearest_neighbour_distances) == 0L) {
    stop(
      "Could not estimate nearest-neighbour spatial spot distance.",
      call. = FALSE
    )
  }
  
  stats::median(nearest_neighbour_distances)
}

# ==============================================================================
# 7. CHAOS
# ==============================================================================

calculate_chaos_score <- function(
  coordinate_matrix,
  cluster_labels,
  spot_pitch,
  normalize_by_spot_pitch = TRUE
) {
  coordinate_matrix <- as.matrix(coordinate_matrix)
  cluster_labels <- as.character(cluster_labels)
  
  if (nrow(coordinate_matrix) != length(cluster_labels)) {
    stop(
      "Coordinate and cluster-label lengths differ in CHAOS calculation.",
      call. = FALSE
    )
  }
  
  if (!is.finite(spot_pitch) || spot_pitch <= 0) {
    stop("Invalid spatial spot pitch.", call. = FALSE)
  }
  
  if (isTRUE(normalize_by_spot_pitch)) {
    coordinate_matrix <- coordinate_matrix / spot_pitch
  }
  
  cluster_levels <- unique(cluster_labels)
  nearest_same_cluster_distances <- rep(NA_real_, nrow(coordinate_matrix))
  
  for (cluster_level in cluster_levels) {
    cluster_indices <- which(cluster_labels == cluster_level)
    
    if (length(cluster_indices) < 2L) next
    
    cluster_coordinates <- coordinate_matrix[
      cluster_indices,
      ,
      drop = FALSE
    ]
    
    nearest_result <- RANN::nn2(
      data = cluster_coordinates,
      query = cluster_coordinates,
      k = 2L
    )
    
    nearest_same_cluster_distances[cluster_indices] <-
      nearest_result$nn.dists[, 2L]
  }
  
  valid_distances <- is.finite(nearest_same_cluster_distances)
  
  if (!any(valid_distances)) return(NA_real_)
  
  mean(nearest_same_cluster_distances[valid_distances])
}

# ==============================================================================
# 8. PAS
# ==============================================================================

precompute_pas_neighbours <- function(coordinate_matrix, k = 10L) {
  coordinate_matrix <- as.matrix(coordinate_matrix)
  
  if (nrow(coordinate_matrix) <= k) {
    stop("PAS requires more than ", k, " spots.", call. = FALSE)
  }
  
  nearest_result <- RANN::nn2(
    data = coordinate_matrix,
    query = coordinate_matrix,
    k = as.integer(k) + 1L
  )
  
  nearest_result$nn.idx[, -1L, drop = FALSE]
}

calculate_pas_score <- function(
  neighbour_indices,
  cluster_labels,
  minimum_different_neighbours = 6L
) {
  cluster_labels <- as.character(cluster_labels)
  
  if (nrow(neighbour_indices) != length(cluster_labels)) {
    stop(
      "PAS neighbour matrix and cluster-label lengths differ.",
      call. = FALSE
    )
  }
  
  different_neighbour_counts <- vapply(
    seq_len(nrow(neighbour_indices)),
    function(spot_index) {
      neighbour_labels <- cluster_labels[
        neighbour_indices[spot_index, , drop = TRUE]
      ]
      sum(neighbour_labels != cluster_labels[[spot_index]])
    },
    integer(1)
  )
  
  abnormal_spots <- different_neighbour_counts >=
    as.integer(minimum_different_neighbours)
  
  mean(abnormal_spots)
}

# ==============================================================================
# 9. Load data
# ==============================================================================

workflow_message("Loading ris3q29_integrate...")
ris3q29_integrate <- load_expected_object(
  rdata_file = integrated_rdata_file,
  expected_object_name = "ris3q29_integrate"
)

workflow_message("Loading ris3q29_st_data...")
ris3q29_st_data <- load_expected_object(
  rdata_file = spatial_rdata_file,
  expected_object_name = "ris3q29_st_data"
)

if (!inherits(ris3q29_integrate, "Seurat")) {
  stop("`ris3q29_integrate` is not a Seurat object.", call. = FALSE)
}

# ==============================================================================
# 10. Clustering table and sample order
# ==============================================================================

cluster_table <- as.data.frame(
  ris3q29_st_data$clusters,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cluster_columns <- c("sample", "barcode")
missing_required_cluster_columns <- setdiff(
  required_cluster_columns,
  colnames(cluster_table)
)

if (length(missing_required_cluster_columns) > 0L) {
  stop(
    "Missing required column(s) in ris3q29_st_data$clusters: ",
    paste(missing_required_cluster_columns, collapse = ", "),
    call. = FALSE
  )
}

cluster_table$sample <- as.character(cluster_table$sample)
cluster_table$barcode <- as.character(cluster_table$barcode)

sample_names <- unique(cluster_table$sample)

if (!is.null(ris3q29_st_data$samples)) {
  candidate_sample_names <- as.character(ris3q29_st_data$samples)
  candidate_sample_names <- candidate_sample_names[
    !is.na(candidate_sample_names) & candidate_sample_names != ""
  ]
  
  if (
    length(candidate_sample_names) > 0L &&
    all(candidate_sample_names %in% sample_names)
  ) {
    sample_names <- unique(candidate_sample_names)
  }
}

if (length(sample_names) == 0L) {
  stop("No sample names were detected.", call. = FALSE)
}

if (anyDuplicated(sample_names) > 0L) {
  stop("Duplicated sample names detected.", call. = FALSE)
}

workflow_message("Samples detected: ", length(sample_names))
workflow_message("Sample order: ", paste(sample_names, collapse = ", "))

if (length(sample_names) != 26L) {
  warning(
    "Expected 26 samples, but found ",
    length(sample_names),
    ". Analysis will continue with all available samples.",
    call. = FALSE
  )
}

# ==============================================================================
# 11. Clustering resolutions
# ==============================================================================

resolution_columns <- grep(
  pattern = "^cluster_resolution_",
  x = colnames(cluster_table),
  value = TRUE
)

if (length(resolution_columns) == 0L) {
  stop("No `cluster_resolution_*` columns were found.", call. = FALSE)
}

resolution_values <- as.numeric(
  sub(
    pattern = "^cluster_resolution_",
    replacement = "",
    x = resolution_columns
  )
)

if (anyNA(resolution_values)) {
  stop("Could not parse one or more clustering resolutions.", call. = FALSE)
}

resolution_order <- order(resolution_values)
resolution_columns <- resolution_columns[resolution_order]
resolution_values <- resolution_values[resolution_order]

workflow_message("Clustering resolutions detected: ", length(resolution_values))
workflow_message(
  "Resolution range: ",
  min(resolution_values),
  " - ",
  max(resolution_values)
)

for (cluster_column in resolution_columns) {
  if (anyNA(cluster_table[[cluster_column]])) {
    stop(
      "NA cluster assignments detected in column `",
      cluster_column,
      "`.",
      call. = FALSE
    )
  }
}

# ==============================================================================
# 12. Construct spot identifiers and align with Seurat object
# ==============================================================================

cluster_table$sample_barcode <- paste(
  cluster_table$sample,
  cluster_table$barcode,
  sep = "_"
)

if (anyDuplicated(cluster_table$sample_barcode) > 0L) {
  stop(
    "Duplicated sample_barcode identifiers detected in cluster table.",
    call. = FALSE
  )
}

seurat_spots <- colnames(ris3q29_integrate)
missing_seurat_spots <- setdiff(cluster_table$sample_barcode, seurat_spots)

if (length(missing_seurat_spots) > 0L) {
  stop(
    length(missing_seurat_spots),
    " cluster-table spots were not found in ris3q29_integrate.\n",
    "First missing IDs:\n",
    paste(head(missing_seurat_spots, 20L), collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# 13. Spatial coordinates
# ==============================================================================

coordinate_table <- as.data.frame(
  ris3q29_st_data$bcs_information,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_coordinate_columns <- c(
  "sample",
  "barcode",
  "imagerow",
  "imagecol"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  colnames(coordinate_table)
)

if (length(missing_coordinate_columns) > 0L) {
  stop(
    "Missing required spatial-coordinate column(s): ",
    paste(missing_coordinate_columns, collapse = ", "),
    call. = FALSE
  )
}

coordinate_table$sample <- as.character(coordinate_table$sample)
coordinate_table$barcode <- as.character(coordinate_table$barcode)
coordinate_table$sample_barcode <- paste(
  coordinate_table$sample,
  coordinate_table$barcode,
  sep = "_"
)

if (anyDuplicated(coordinate_table$sample_barcode) > 0L) {
  stop(
    "Duplicated sample_barcode identifiers detected in bcs_information.",
    call. = FALSE
  )
}

coordinate_match <- match(
  cluster_table$sample_barcode,
  coordinate_table$sample_barcode
)

if (anyNA(coordinate_match)) {
  missing_ids <- cluster_table$sample_barcode[is.na(coordinate_match)]
  stop(
    length(missing_ids),
    " clustered spots do not have spatial coordinates.\n",
    "First missing IDs:\n",
    paste(head(missing_ids, 20L), collapse = "\n"),
    call. = FALSE
  )
}

cluster_table$imagecol <- as.numeric(
  coordinate_table$imagecol[coordinate_match]
)
cluster_table$imagerow <- as.numeric(
  coordinate_table$imagerow[coordinate_match]
)

if (anyNA(cluster_table$imagecol) || anyNA(cluster_table$imagerow)) {
  stop("NA values found in spatial coordinates.", call. = FALSE)
}

# ==============================================================================
# 14. Transcriptomic embedding
# ==============================================================================

available_reductions <- SeuratObject::Reductions(ris3q29_integrate)

if (!transcriptomic_reduction %in% available_reductions) {
  stop(
    "Requested transcriptomic reduction `",
    transcriptomic_reduction,
    "` is absent.\nAvailable reductions: ",
    paste(available_reductions, collapse = ", "),
    call. = FALSE
  )
}

transcriptomic_embedding <- SeuratObject::Embeddings(
  ris3q29_integrate[[transcriptomic_reduction]]
)

if (max(transcriptomic_dims) > ncol(transcriptomic_embedding)) {
  stop(
    "Requested transcriptomic dimension ",
    max(transcriptomic_dims),
    " but only ",
    ncol(transcriptomic_embedding),
    " dimensions are available in `",
    transcriptomic_reduction,
    "`.",
    call. = FALSE
  )
}

missing_embedding_spots <- setdiff(
  cluster_table$sample_barcode,
  rownames(transcriptomic_embedding)
)

if (length(missing_embedding_spots) > 0L) {
  stop(
    length(missing_embedding_spots),
    " clustered spots were not found in the PCA embedding.",
    call. = FALSE
  )
}

workflow_message("Transcriptomic ASW reduction: ", transcriptomic_reduction)
workflow_message(
  "Transcriptomic dimensions: ",
  min(transcriptomic_dims),
  "-",
  max(transcriptomic_dims)
)

# ==============================================================================
# 15. Calculate metrics
# ==============================================================================

analysis_started <- Sys.time()

expected_result_rows <-
  length(sample_names) * length(resolution_values) * length(method_order)

result_rows <- vector("list", expected_result_rows)
result_index <- 1L

for (sample_index in seq_along(sample_names)) {
  sample_started <- Sys.time()
  sample_id <- sample_names[[sample_index]]
  
  workflow_message("============================================================")
  workflow_message(
    "SAMPLE ",
    sample_index,
    "/",
    length(sample_names),
    ": ",
    sample_id
  )
  
  sample_table <- cluster_table[
    cluster_table$sample == sample_id,
    ,
    drop = FALSE
  ]
  rownames(sample_table) <- NULL
  
  n_sample_spots <- nrow(sample_table)
  workflow_message("Spots: ", n_sample_spots)
  
  if (n_sample_spots <= pas_k) {
    stop(
      "Sample ",
      sample_id,
      " contains too few spots for PAS.",
      call. = FALSE
    )
  }
  
  spatial_coordinates <- as.matrix(
    sample_table[, c("imagecol", "imagerow"), drop = FALSE]
  )
  storage.mode(spatial_coordinates) <- "double"
  
  sample_embedding <- transcriptomic_embedding[
    sample_table$sample_barcode,
    transcriptomic_dims,
    drop = FALSE
  ]
  storage.mode(sample_embedding) <- "double"
  
  transcriptomic_indices <- select_common_validation_indices(
    cluster_table = sample_table,
    cluster_columns = resolution_columns,
    max_spots = transcriptomic_asw_max_spots_per_sample,
    minimum_spots_per_cluster =
      transcriptomic_asw_minimum_spots_per_cluster,
    seed = random_seed + sample_index
  )
  
  spatial_indices <- select_common_validation_indices(
    cluster_table = sample_table,
    cluster_columns = resolution_columns,
    max_spots = spatial_asw_max_spots_per_sample,
    minimum_spots_per_cluster = spatial_asw_minimum_spots_per_cluster,
    seed = random_seed + 1000L + sample_index
  )
  
  workflow_message(
    "Transcriptomic ASW spots: ",
    length(transcriptomic_indices)
  )
  workflow_message("Spatial ASW spots: ", length(spatial_indices))
  
  workflow_message("Precomputing transcriptomic distance matrix...")
  transcriptomic_distance <- stats::dist(
    sample_embedding[transcriptomic_indices, , drop = FALSE],
    method = "euclidean"
  )
  
  workflow_message("Precomputing spatial distance matrix...")
  spatial_distance <- stats::dist(
    spatial_coordinates[spatial_indices, , drop = FALSE],
    method = "euclidean"
  )
  
  spot_pitch <- estimate_spatial_spot_pitch(spatial_coordinates)
  workflow_message(
    "Median spatial spot pitch: ",
    format(spot_pitch, digits = 6)
  )
  
  workflow_message("Precomputing PAS spatial neighbours...")
  pas_neighbour_indices <- precompute_pas_neighbours(
    coordinate_matrix = spatial_coordinates,
    k = pas_k
  )
  
  for (resolution_index in seq_along(resolution_values)) {
    resolution_started <- Sys.time()
    resolution_value <- resolution_values[[resolution_index]]
    cluster_column <- resolution_columns[[resolution_index]]
    cluster_labels <- as.character(sample_table[[cluster_column]])
    
    n_clusters <- length(unique(cluster_labels))
    
    transcriptomic_asw <- calculate_asw_from_distance(
      distance_object = transcriptomic_distance,
      cluster_labels = cluster_labels[transcriptomic_indices]
    )
    
    spatial_asw <- calculate_asw_from_distance(
      distance_object = spatial_distance,
      cluster_labels = cluster_labels[spatial_indices]
    )
    
    chaos_value <- calculate_chaos_score(
      coordinate_matrix = spatial_coordinates,
      cluster_labels = cluster_labels,
      spot_pitch = spot_pitch,
      normalize_by_spot_pitch = normalize_chaos_by_spot_pitch
    )
    
    pas_value <- calculate_pas_score(
      neighbour_indices = pas_neighbour_indices,
      cluster_labels = cluster_labels,
      minimum_different_neighbours = pas_minimum_different_neighbours
    )
    
    current_values <- c(
      transcriptomic_ASW = transcriptomic_asw,
      spatial_ASW = spatial_asw,
      CHAOS = chaos_value,
      PAS = pas_value
    )
    
    for (method_name in method_order) {
      result_rows[[result_index]] <- data.frame(
        sample_ID = sample_id,
        resolution = resolution_value,
        method = method_name,
        value = as.numeric(current_values[[method_name]]),
        nSpots = n_sample_spots,
        nClusters = n_clusters,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      result_index <- result_index + 1L
    }
    
    resolution_finished <- Sys.time()
    
    workflow_message(
      "resolution=",
      formatC(resolution_value, format = "f", digits = 2L),
      " | clusters=",
      n_clusters,
      " | transcriptomic_ASW=",
      formatC(transcriptomic_asw, format = "f", digits = 4L),
      " | spatial_ASW=",
      formatC(spatial_asw, format = "f", digits = 4L),
      " | CHAOS=",
      formatC(chaos_value, format = "f", digits = 4L),
      " | PAS=",
      formatC(pas_value, format = "f", digits = 4L),
      " | elapsed=",
      format_elapsed_time(
        difftime(
          resolution_finished,
          resolution_started,
          units = "secs"
        )
      )
    )
  }
  
  rm(
    transcriptomic_distance,
    spatial_distance,
    sample_embedding,
    spatial_coordinates,
    pas_neighbour_indices
  )
  invisible(gc(verbose = FALSE))
  
  sample_finished <- Sys.time()
  workflow_message(
    "Finished sample ",
    sample_id,
    " | elapsed=",
    format_elapsed_time(
      difftime(sample_finished, sample_started, units = "secs")
    )
  )
}

result_rows <- result_rows[
  !vapply(result_rows, is.null, logical(1))
]

results_long <- do.call(rbind, result_rows)
rownames(results_long) <- NULL

expected_rows <-
  length(sample_names) * length(resolution_values) * length(method_order)

if (nrow(results_long) != expected_rows) {
  stop(
    "Unexpected number of result rows.\nExpected: ",
    expected_rows,
    "\nObserved: ",
    nrow(results_long),
    call. = FALSE
  )
}

# ==============================================================================
# 16. Wide result table
# ==============================================================================

wide_grid <- do.call(
  rbind,
  lapply(
    resolution_values,
    function(resolution_value) {
      data.frame(
        resolution = rep(resolution_value, length(method_order)),
        method = method_order,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)
rownames(wide_grid) <- NULL

results_wide <- wide_grid

result_key <- paste(
  results_long$resolution,
  results_long$method,
  results_long$sample_ID,
  sep = "\r"
)

for (sample_id in sample_names) {
  target_key <- paste(
    wide_grid$resolution,
    wide_grid$method,
    sample_id,
    sep = "\r"
  )
  
  matched_rows <- match(target_key, result_key)
  
  if (anyNA(matched_rows)) {
    stop(
      "Missing result values while constructing wide output for sample ",
      sample_id,
      ".",
      call. = FALSE
    )
  }
  
  results_wide[[sample_id]] <- results_long$value[matched_rows]
}

sample_value_matrix <- as.matrix(
  results_wide[, sample_names, drop = FALSE]
)

results_wide$median <- apply(sample_value_matrix, 1L, safe_median)
results_wide$mean <- apply(sample_value_matrix, 1L, safe_mean)
results_wide$SD <- apply(sample_value_matrix, 1L, safe_sd)

# ==============================================================================
# 17. TSV
# ==============================================================================

workflow_message("Writing TSV...")
write_tsv(results_wide, output_tsv)

# ==============================================================================
# 18. XLSX
# ==============================================================================

workflow_message("Writing XLSX...")

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, sheetName = "PAS_CHAOS_ASW")

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

numeric_style <- openxlsx::createStyle(numFmt = "0.00000")

openxlsx::writeData(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  x = results_wide,
  startRow = 1L,
  startCol = 1L,
  headerStyle = header_style,
  withFilter = TRUE
)

openxlsx::freezePane(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  firstRow = TRUE,
  firstActiveCol = 3L
)

openxlsx::addStyle(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  style = numeric_style,
  rows = 2L:(nrow(results_wide) + 1L),
  cols = c(1L, 3L:ncol(results_wide)),
  gridExpand = TRUE,
  stack = TRUE
)

openxlsx::setColWidths(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  cols = 1L,
  widths = 12
)
openxlsx::setColWidths(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  cols = 2L,
  widths = 22
)
openxlsx::setColWidths(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  cols = 3L:(2L + length(sample_names)),
  widths = 14
)
openxlsx::setColWidths(
  wb = workbook,
  sheet = "PAS_CHAOS_ASW",
  cols = (ncol(results_wide) - 1L):ncol(results_wide),
  widths = 12
)

openxlsx::saveWorkbook(
  wb = workbook,
  file = output_xlsx,
  overwrite = TRUE
)

# ==============================================================================
# 19. Plot data
# ==============================================================================

plot_summary <- results_wide[
  ,
  c(
    "resolution",
    "method",
    "median",
    "mean",
    "SD"
  ),
  drop = FALSE
]

# All tested resolutions are retained in the full PDF. To keep the x axis
# readable, tick labels are displayed every 0.10 when possible.
full_x_breaks <- resolution_values[
  abs(
    resolution_values * 10 -
      round(resolution_values * 10)
  ) < 1e-8
]

if (length(full_x_breaks) == 0L) {
  full_x_breaks <- resolution_values
}

# The lighter PDF contains only resolutions that are exact multiples of 0.10.
resolution_step_010_values <- resolution_values[
  abs(
    resolution_values * 10 -
      round(resolution_values * 10)
  ) < 1e-8
]

if (length(resolution_step_010_values) == 0L) {
  stop(
    "No clustering resolutions matching a 0.10 step were found.",
    call. = FALSE
  )
}

# ==============================================================================
# 20. Plot helper
# ==============================================================================

create_metric_summary_plot <- function(
  summary_table,
  metric_name,
  plot_title,
  y_axis_label,
  direction_label,
  x_breaks,
  y_limits = NULL
) {
  
  metric_table <- summary_table[
    summary_table$method == metric_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(metric_table) == 0L) {
    stop(
      "No rows available for metric: ",
      metric_name,
      call. = FALSE
    )
  }
  
  output_plot <- ggplot2::ggplot(
    metric_table,
    ggplot2::aes(
      x = resolution,
      y = mean
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7,
      colour = "black",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 2.6,
      colour = "black",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = paste0(
        direction_label,
        " | line and points = mean across samples"
      ),
      x = "Clustering resolution",
      y = y_axis_label
    ) +
    ggplot2::theme_bw(
      base_size = 13
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = ggplot2::element_text(
        size = 11
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(y_limits)) {
    output_plot <- output_plot +
      ggplot2::coord_cartesian(
        ylim = y_limits
      )
  }
  
  output_plot
}

# ==============================================================================
# 21. Four plots: all tested resolutions
# ==============================================================================

transcriptomic_asw_plot_all <- create_metric_summary_plot(
  summary_table = plot_summary,
  metric_name = "transcriptomic_ASW",
  plot_title = "Transcriptomic ASW by clustering resolution",
  y_axis_label = "Transcriptomic ASW",
  direction_label = "Higher = better",
  x_breaks = full_x_breaks,
  y_limits = c(0, 0.5)
)

spatial_asw_plot_all <- create_metric_summary_plot(
  summary_table = plot_summary,
  metric_name = "spatial_ASW",
  plot_title = "Spatial ASW by clustering resolution",
  y_axis_label = "Spatial ASW",
  direction_label = "Higher = better",
  x_breaks = full_x_breaks,
  y_limits = NULL
)

chaos_plot_all <- create_metric_summary_plot(
  summary_table = plot_summary,
  metric_name = "CHAOS",
  plot_title = "CHAOS by clustering resolution",
  y_axis_label = "CHAOS",
  direction_label = "Lower = better",
  x_breaks = full_x_breaks,
  y_limits = NULL
)

pas_plot_all <- create_metric_summary_plot(
  summary_table = plot_summary,
  metric_name = "PAS",
  plot_title = "PAS by clustering resolution",
  y_axis_label = "PAS",
  direction_label = "Lower = better",
  x_breaks = full_x_breaks,
  y_limits = c(0, 1)
)

# ==============================================================================
# 22. Four plots: resolutions in steps of 0.10 only
# ==============================================================================

plot_summary_step_010 <- plot_summary[
  plot_summary$resolution %in%
    resolution_step_010_values,
  ,
  drop = FALSE
]

transcriptomic_asw_plot_step_010 <- create_metric_summary_plot(
  summary_table = plot_summary_step_010,
  metric_name = "transcriptomic_ASW",
  plot_title = "Transcriptomic ASW by clustering resolution",
  y_axis_label = "Transcriptomic ASW",
  direction_label = "Higher = better",
  x_breaks = resolution_step_010_values,
  y_limits = c(0, 0.5)
)

spatial_asw_plot_step_010 <- create_metric_summary_plot(
  summary_table = plot_summary_step_010,
  metric_name = "spatial_ASW",
  plot_title = "Spatial ASW by clustering resolution",
  y_axis_label = "Spatial ASW",
  direction_label = "Higher = better",
  x_breaks = resolution_step_010_values,
  y_limits = NULL
)

chaos_plot_step_010 <- create_metric_summary_plot(
  summary_table = plot_summary_step_010,
  metric_name = "CHAOS",
  plot_title = "CHAOS by clustering resolution",
  y_axis_label = "CHAOS",
  direction_label = "Lower = better",
  x_breaks = resolution_step_010_values,
  y_limits = NULL
)

pas_plot_step_010 <- create_metric_summary_plot(
  summary_table = plot_summary_step_010,
  metric_name = "PAS",
  plot_title = "PAS by clustering resolution",
  y_axis_label = "PAS",
  direction_label = "Lower = better",
  x_breaks = resolution_step_010_values,
  y_limits = c(0, 1)
)

# ==============================================================================
# 23. Save PDF 1: all tested resolutions, four pages
# ==============================================================================

workflow_message(
  "Writing four-page PDF for all tested resolutions..."
)

grDevices::pdf(
  file = output_pdf_all_resolutions,
  width = 10,
  height = 7.5,
  onefile = TRUE
)

print(transcriptomic_asw_plot_all)
print(spatial_asw_plot_all)
print(chaos_plot_all)
print(pas_plot_all)

grDevices::dev.off()

# ==============================================================================
# 24. Save PDF 2: resolution step 0.10, four pages
# ==============================================================================

workflow_message(
  "Writing four-page PDF for resolutions in steps of 0.10..."
)

grDevices::pdf(
  file = output_pdf_resolution_step_010,
  width = 10,
  height = 7.5,
  onefile = TRUE
)

print(transcriptomic_asw_plot_step_010)
print(spatial_asw_plot_step_010)
print(chaos_plot_step_010)
print(pas_plot_step_010)

grDevices::dev.off()

# ==============================================================================
# 25. Final validation
# ==============================================================================

required_output_files <- c(
  output_tsv,
  output_xlsx,
  output_pdf_all_resolutions,
  output_pdf_resolution_step_010
)

missing_output_files <- required_output_files[
  !file.exists(required_output_files)
]

if (length(missing_output_files) > 0L) {
  stop(
    "One or more expected output files were not created:\n",
    paste(
      missing_output_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

analysis_finished <- Sys.time()

workflow_message("============================================================")
workflow_message("CLUSTERING VALIDATION FINISHED")
workflow_message("============================================================")
workflow_message("Samples: ", length(sample_names))
workflow_message("Resolutions: ", length(resolution_values))
workflow_message("Methods: ", paste(method_order, collapse = ", "))
workflow_message("Output rows: ", nrow(results_wide))
workflow_message(
  "Expected output rows: ",
  length(resolution_values) * length(method_order)
)
workflow_message(
  "Resolution-step-0.10 values: ",
  paste(
    formatC(
      resolution_step_010_values,
      format = "f",
      digits = 2L
    ),
    collapse = ", "
  )
)
workflow_message(
  "Elapsed: ",
  format_elapsed_time(
    difftime(
      analysis_finished,
      analysis_started,
      units = "secs"
    )
  )
)
workflow_message("")
workflow_message("TSV:  ", output_tsv)
workflow_message("XLSX: ", output_xlsx)
workflow_message("PDF all resolutions: ", output_pdf_all_resolutions)
workflow_message(
  "PDF resolution step 0.10: ",
  output_pdf_resolution_step_010
)
workflow_message("============================================================")
