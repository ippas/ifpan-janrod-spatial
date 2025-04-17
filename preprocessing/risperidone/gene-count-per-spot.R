
risperidone_summary_statistics_half$raw_data


risperidone_summary_statistics_half$raw_data$resolution_0.05$cluster_0$control$mean

risperidone_integrate_half$nFeature_RNA %>% mean()

risperidone_integrate_half@assays$RNA@data %>% dim



peak_gene_mapping <- info_peaks_risperidone %>% 
  select(peak_id, gene_name) %>% 
  mutate(peak_id = str_replace_all(peak_id, "_", "-"))

risperidone_integrate_half@assays$integrated@data %>% 
  .[1:10, 1:10] %>% apply(., 2, is.na)


info_peaks_risperidone %>% 
  select(peak_id, gene_name) %>% 
  head %>% 
  mutate(peak_id = str_replace_all(peak_id, "-_", "-"))

risperidone_integrate_half@assays$integrated@data %>% 
  .[1:10, 1:10] %>% 
  {colSums(. != 0)} %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Column_Name")

library(progressr)

risperidone_integrate_half@assays$RNA@data %>% dim


gene_sum_per_spot <- risperidone_integrate_half@assays$RNA@data %>%
  .[, 1:100] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "peak_id") %>%
  left_join(peak_gene_mapping, by = "peak_id") %>% 
  select(-peak_id) %>% 
  group_by(gene_name) %>%
  summarise(across(everything(), ~ sum(. != 0), .names = "sum_{.col}")) %>%
  summarise(across(starts_with("sum_"), sum)) # Sumowanie dla każdej kolumny

gene_sum_per_spot


################################################################################
# Define the number of columns processed in each iteration
chunk_size <- 100
num_cols <- ncol(risperidone_integrate_half@assays$RNA@data)
total_steps <- ceiling(num_cols / chunk_size) # Calculate the total number of iterations

# Initialize the progress bar
handlers(global = TRUE)
progress <- progressor(steps = total_steps)

# Create an empty list to store results
gene_sum_list <- list()

# Loop through the data in chunks
for (start_col in seq(1, num_cols, by = chunk_size)) {
  
  end_col <- min(start_col + chunk_size - 1, num_cols) # Determine the column range for this iteration
  
  print(start_col/num_cols)
  
  # Update the progress bar
  progress(sprintf("Processing columns %d-%d", start_col, end_col))
  
  # Extract and process a subset of columns
  gene_sum_per_spot <- risperidone_integrate_half@assays$RNA@data %>%
    .[, start_col:end_col] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "peak_id") %>%
    left_join(peak_gene_mapping, by = "peak_id") %>% 
    select(-peak_id) %>% 
    group_by(gene_name) %>%
    summarise(across(everything(), ~ sum(. != 0), .names = "sum_{.col}")) %>%
    summarise(across(starts_with("sum_"), sum)) # Summing across all columns
  
  # Store results in the list
  gene_sum_list[[length(gene_sum_list) + 1]] <- gene_sum_per_spot
  
  print(gene_sum_per_spot)
}

# Merge all results into a single dataframe
final_gene_sum <- reduce(gene_sum_list, full_join, by = "gene_name")

gene_sum_list %>% 
  unlist() %>% median

samples_risperidone
samples_saline

gene_sum_list %>%
  unlist() %>%
  tibble::enframe(name = "sample_id", value = "value") %>%
  mutate(sample_group = str_extract(sample_id, "S[0-9]+Nr[0-9]+")) %>%
  group_by(sample_group) %>%
  summarise(mean_value = mean(value))


gene_sum_list %>%
  unlist() %>%
  enframe(name = "sample_id", value = "value") %>%
  mutate(sample_group = str_extract(sample_id, "S[0-9]+Nr[0-9]+"),
         condition = case_when(
           sample_group %in% samples_risperidone ~ "risperidone",
           TRUE ~ "saline"
         )) %>%
  group_by(sample_group, condition) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    median_value = median(value),
    .groups = "drop"
  )
