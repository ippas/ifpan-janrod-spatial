metadata_ris3q29 

"S13839Nr3"

risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_2$experiment$sum %>% head
risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_2$control$sum %>% head

  
  
risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_0$control$sum %>%
  as.data.frame() %>% 
  rownames_to_column(var = "peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "count") %>%
  # filter(count < 300, count > 0) %>%  # usuń zera, bo log(0) jest undefined
  ggplot(aes(x = count)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ sample, scales = "fixed", nrow = 1) +
  theme_minimal() +
  labs(title = "Log-scaled histogram of peak counts per sample",
       x = "Count (log10 scale)", y = "Frequency (log10 scale)")


risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_0$control$sum %>% 
  apply(., 2, sum)

risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_0$control$sum %>% 
  apply(., 2, sum) %>%
  { . / max(.) }

library(forcats)

font_size <- 14

# Total transcript count per sample
apply(risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_0$control$sum, 2, sum) %>%
  enframe(name = "sample", value = "total_count") %>%
  arrange(desc(total_count)) %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  ggplot(aes(x = sample, y = total_count)) +
  geom_col(width = 0.7, fill = "tomato", color = "black") +
  theme_minimal(base_size = font_size) +
  labs(title = "Total transcript count per sample", x = "Sample", y = "Total count") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = font_size + 2, face = "bold")
  ) -> p1

# Mean ± SE per sample
risDelSalDel_summary_statistics$raw_data$resolution_0.4$cluster_0$control$sum %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  summarise(
    mean_count = mean(count),
    se_count = sd(count) / sqrt(n())
  ) %>%
  arrange(desc(mean_count)) %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  ggplot(aes(x = sample, y = mean_count)) +
  geom_col(width = 0.7, fill = "steelblue", color = "black") +
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  theme_minimal(base_size = font_size) +
  labs(title = "Mean transcript count per sample ± SE", x = "Sample", y = "Mean count") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = font_size + 2, face = "bold")
  ) -> p2

p1 + p2


risDelSalDel_summary_statistics$raw_data$resolution_0.4



paste0("cluster_", 0:19) %>%
set_names() %>%  # nazwy jako nazwy listy
imap_dfr(function(clust_data, clust_name) {
  risDelSalDel_summary_statistics$raw_data$resolution_0.4[[clust_name]]$control$sum %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcript") %>%
    pivot_longer(-transcript, names_to = "sample", values_to = "count") %>%
    mutate(cluster = clust_name)
}) %>% 
  group_by(sample) %>%  
  nest %>% 
  mutate(sum = map(data, ~ .x$count %>% sum(na.rm = TRUE))) %>% unnest(sum)


paste0("cluster_", 0:19) %>%
  set_names() %>%
  imap_dfr(function(clust_name, clust_label) {
    risDelSalDel_summary_statistics$raw_data$resolution_0.4[[clust_label]]$control$sum %>%
      as.data.frame() %>%
      rownames_to_column(var = "transcript") %>%
      pivot_longer(-transcript, names_to = "sample", values_to = "count") %>%
      mutate(cluster = clust_label)
  }) %>%
  group_by(transcript) %>% 
  filter(sample == "S13839Nr3") %>% 
  nest %>% 
  mutate(sum_transcript = map(data, ~.x$count %>% sum(na.rm = TRUE))) %>% 
  unnest(sum_transcript) %>% 
  arrange(desc(sum_transcript)) %>% 
  select(-data) %>% 
  as.data.frame() %>% head(30)
  
  
spatial_cluster_peak_barcode_plot <- function(spatial_data,
                                              summary_statistics,
                                              resolution = "resolution_0.4",
                                              cluster_id,
                                              peak_id,
                                              samples,
                                              size = 1.2,
                                              alpha = 1,
                                              tif_image = TRUE,
                                              normalization = TRUE,
                                              show_legend = TRUE,
                                              return_list = FALSE) {
  # 1. Ekspresja: pobieramy sumy dla peak_id w danym klastrze
  peak_vector <- summary_statistics$raw_data[[resolution]][[paste0("cluster_", cluster_id)]]$sum[peak_id, , drop = TRUE]
  
  peak_df <- tibble(sample = names(peak_vector), value = as.numeric(peak_vector)) %>%
    filter(sample %in% samples)
  
  # 2. Pobierz barcode'y należące do klastra z spatial_data
  barcode_in_cluster <- spatial_data$clusters[[resolution]] %>%
    filter(sample %in% samples, cluster == as.numeric(cluster_id)) %>%
    select(sample, barcode)
  
  # 3. Przypisz wartość ekspresji dla każdego barcode
  barcode_values <- barcode_in_cluster %>%
    left_join(peak_df, by = "sample")
  
  # 4. Normalizacja (jeśli trzeba)
  if (normalization) {
    barcode_values <- barcode_values %>%
      mutate(value = (value - min(value, na.rm = TRUE)) /
               (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)))
  }
  
  # 5. Dołącz współrzędne (tylko dla wybranych barcode'ów)
  plotting_df <- spatial_data$bcs_information %>%
    inner_join(barcode_values, by = c("sample", "barcode"))
  
  # 6. Rysuj
  plot_list <- list()
  
  for (sample_id in samples) {
    sample_plot <- plotting_df %>%
      filter(sample == sample_id) %>%
      ggplot(aes(x = imagecol, y = imagerow, fill = value)) +
      {if (tif_image) {
        geom_spatial(data = spatial_data$images_information[spatial_data$images_information$sample == sample_id, ],
                     aes(grob = grob),
                     x = 0.5, y = 0.5)
      }} +
      geom_point(shape = 21, colour = "black", size = size, stroke = 0, alpha = alpha) +
      coord_cartesian(expand = FALSE) +
      scale_fill_gradientn(colours = my_custom_palette(100),
                           limits = if (normalization) c(0, 1) else NULL) +
      xlim(0, max(plotting_df %>% filter(sample == sample_id) %>% pull(width))) +
      ylim(max(plotting_df %>% filter(sample == sample_id) %>% pull(height)), 0) +
      xlab("") + ylab("") +
      ggtitle(paste0(
        sample_id, ": ",
        spatial_data$sample_information %>%
          filter(sample_ID == sample_id) %>%
          transmute(label = paste0(treatment, ", ", mouse_genotype)) %>%
          pull(label)
      )) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    
    if (!show_legend) {
      sample_plot <- sample_plot + theme(legend.position = "none")
    }
    
    plot_list[[sample_id]] <- sample_plot
  }
  
  if (return_list) {
    return(plot_list)
  } else {
    return(wrap_plots(plot_list))
  }
}


spatial_cluster_peak_barcode_plot(
  spatial_data = ris3q29_st_data,
  summary_statistics = risDelSalDel_summary_statistics,
  resolution = "resolution_0.4",
  cluster_id = "2",
  peak_id = "peak-305-Xkr4",
  samples = c("S13839Nr3")
)

risDelSalDel_summary_statistics$raw_data$resolution_0.1$cluster_0$control$sum %>% head
# 

spatial_cluster_peak_barcode_plot(
  spatial_data = risDelSalDel_summary_statistics,
  summary_statistics = risDelSalDel_summary_statistics,
  barcode_cluster_table = barcode_clusters_df,
  resolution = "resolution_0.4",
  cluster_id = "2",
  peak_id = "peak-305-Xkr4",
  samples = c("S13839Nr2", "S13839Nr7")
)

spatial_cluster_peak_barcode_plot(
       spatial_data = ris3q29_st_data,
       summary_statistics = risDelSalDel_summary_statistics,
       resolution = "resolution_0.4",
       cluster_id = c("2"),
       peak_id = "peak-305-Xkr4",
       samples = c("S13839Nr3")
)
