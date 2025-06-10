# Debug-wersja funkcji obliczającej summary dla pojedynczego klastra
compute_data_summary_cluster_debug <- function(cluster,
                                               spatial_data,
                                               resolution = 0.4,
                                               trim = 0.05,
                                               control_samples,
                                               experiment_samples,
                                               data_type = "raw_data",
                                               min_number_spots = 20,
                                               metrics = c("mean", "median")) {
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  tryCatch({
    # 1) Wyciągnij wszystkie barcody z danego klastra
    all_barcodes <- spatial_data[[data_type]]$metadata %>%
      right_join(spatial_data$clusters, by = c("barcode", "sample")) %>%
      filter((!!sym(resolution_column)) == cluster) %>%
      mutate(
        sample_barcode = paste(sample, barcode, sep = "_"),
        group = ifelse(sample %in% control_samples, "control", "experiment")
      ) %>%
      select(group, sample_barcode)
    
    # 2) Podziel na control/experiment
    barcodes_list     <- split(all_barcodes$sample_barcode, all_barcodes$group)
    control_barcodes  <- barcodes_list[["control"]]
    experiment_barcodes <- barcodes_list[["experiment"]]
    
    message(">>> Debug cluster ", cluster, "@ res=", resolution)
    message("    # control barcodes:    ", length(control_barcodes))
    message("    # experiment barcodes: ", length(experiment_barcodes))
    
    # 3) Spróbuj wyciągnąć macierz kontrolną
    ctrl_mat <- spatial_data[[data_type]]$data[, control_barcodes, drop = FALSE]
    message("    dim(ctrl_mat): ", paste(dim(ctrl_mat), collapse = " x "))
    
    # 4) Spróbuj wyciągnąć macierz eksperymentalną
    exp_mat <- spatial_data[[data_type]]$data[, experiment_barcodes, drop = FALSE]
    message("    dim(exp_mat):  ", paste(dim(exp_mat), collapse = " x "))
    
    # 5) Dla pewności zerknij na kilka pierwszych wartości
    if (ncol(ctrl_mat) > 0) {
      message("    ctrl_mat[1:3,1:3]:")
      print(head(ctrl_mat[1:min(3,nrow(ctrl_mat)), 1:min(3,ncol(ctrl_mat))]))
    }
    if (ncol(exp_mat) > 0) {
      message("    exp_mat[1:3,1:3]:")
      print(head(exp_mat[1:min(3,nrow(exp_mat)), 1:min(3,ncol(exp_mat))]))
    }
    
    # --- (tutaj możesz wstawić dalsze kroki obliczeń i zwrócić dowolny obiekt)
    # Na potrzeby debugowania zwróć listę z wymiarami
    return(list(
      cluster = cluster,
      dim_control   = dim(ctrl_mat),
      dim_experiment = dim(exp_mat)
    ))
    
  }, error = function(e) {
    stop(sprintf("‼️ BŁĄD w cluster %s: %s", cluster, e$message))
  })
}

# Przykładowe wywołanie debug dla klastra 4:
result_debug <- compute_data_summary_cluster_debug(
  cluster            = 4,
  spatial_data       = ris3q29_st_data,
  resolution         = 0.4,
  trim               = 0.05,
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  data_type          = "raw_data",
  min_number_spots   = 20,
  metrics            = c("mean", "median")
)

print(result_debug)



library(dplyr)

# Debug-wersja funkcji testów dla pojedynczego klastra
perform_statistical_tests_cluster_debug <- function(cluster_data,
                                                    cluster_name,
                                                    spatial_data,
                                                    resolution = 0.4,
                                                    metric = "mean",
                                                    mean_threshold = 0,
                                                    control_samples,
                                                    experiment_samples,
                                                    quantile_normalization = FALSE) {
  resolution_column <- paste0("cluster_resolution_", resolution)
  
  tryCatch({
    message(">>> Debug tests for ", cluster_name, "@ res=", resolution, " metric=", metric)
    
    # 1) Wyciągnij macierze kontrolne/eksperymentalne
    control_expression  <- cluster_data$control[[metric]][, control_samples, drop=FALSE]
    experiment_expression <- cluster_data$experiment[[metric]][, experiment_samples, drop=FALSE]
    
    message("    dim(control_expression):    ", paste(dim(control_expression), collapse=" x "))
    message("    dim(experiment_expression): ", paste(dim(experiment_expression), collapse=" x "))
    
    # 2) Pokaż kilka pierwszych wartości
    if (nrow(control_expression)>0 && ncol(control_expression)>0) {
      message("    control_expression[1:3,1:3]:")
      print(head(control_expression[1:3,1:3]))
    }
    if (nrow(experiment_expression)>0 && ncol(experiment_expression)>0) {
      message("    experiment_expression[1:3,1:3]:")
      print(head(experiment_expression[1:3,1:3]))
    }
    
    # 3) Sprawdź ścieżki NA/low_mean
    na_control  <- all(is.na(control_expression))
    na_experiment <- all(is.na(experiment_expression))
    message("    all NA control? ", na_control, 
            " | all NA experiment? ", na_experiment)
    message("    mean(control): ", mean(control_expression, na.rm=TRUE),
            " | mean(experiment): ", mean(experiment_expression, na.rm=TRUE))
    
    # 4) Jeśli quantile_normalization, zrób krok normalizacji
    if (quantile_normalization) {
      message("    -> applying quantile normalization")
      expr_raw <- cbind(control_expression, experiment_expression)
      expr_norm <- normalize.quantiles(expr_raw)
      message("    dim(norm): ", paste(dim(expr_norm), collapse=" x "))
    }
    
    # 5) Spróbuj wykonać test dla pierwszego genu
    if (!na_control && !na_experiment && 
        mean(control_expression[1, ], na.rm=TRUE) >= mean_threshold &&
        mean(experiment_expression[1, ], na.rm=TRUE) >= mean_threshold) {
      t1 <- try(t.test(control_expression[1, ], experiment_expression[1, ], var.equal=TRUE), silent=TRUE)
      message("    t-test for first row: ",
              if (inherits(t1, "htest")) t1$p.value else "ERROR")
    } else {
      message("    Skipped t-test for first row (NA/low_mean)")
    }
    
    # Zwróć informację o statusie
    return(list(
      cluster = cluster_name,
      control_dim   = dim(control_expression),
      experiment_dim = dim(experiment_expression),
      na_control    = na_control,
      na_experiment = na_experiment
    ))
    
  }, error = function(e) {
    stop(sprintf("‼️ BŁĄD w perform_statistical_tests dla %s: %s", cluster_name, e$message))
  })
}

# Przykładowe wywołanie debug tylko dla cluster_4:
summary_data   <- compute_data_summary( # uprzednio obliczone
  spatial_data       = ris3q29_st_data,
  resolution         = 0.4,
  trim               = 0.05,
  num_cores          = 20,
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  data_type          = "raw_data",
  min_number_spots   = 20,
  metrics            = c("mean", "median")
)

cluster4_data <- summary_data[["cluster_4"]]

result_test_debug <- perform_statistical_tests_cluster_debug(
  cluster_data       = cluster4_data,
  cluster_name       = "cluster_4",
  spatial_data       = ris3q29_st_data,
  resolution         = 0.4,
  metric             = "mean",
  mean_threshold     = 0,
  control_samples    = samples_wt,
  experiment_samples = samples_del,
  quantile_normalization = TRUE
)

print(result_test_debug)
