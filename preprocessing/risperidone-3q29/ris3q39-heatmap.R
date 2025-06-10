ris3q29_integrate

ris3q29_st_data$sample_information %>% dplyr::select(c(sample_ID, treatment, mouse_genotype))

sample_info <- ris3q29_st_data$sample_information %>%
  select(sample_ID, treatment, mouse_genotype)


# preapre data
load("/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/wtDel_summary_statistics.RData")



tmp <- wtDel_summary_statistics$quantile_normalize$resolution_0.4

genes_del3q29 <- c(
  "Bdh1", "Dlg1", "Mfi2", "Pigz", "Ncbp2", "Senp5", "Pak2", "Pigx",
  "Cep19", "Nrros", "Bex6", "Fbxo45", "Wdr53", "Smco1", "Rnf168",
  "Ubxn7", "Tm4sf19", "Tctex1d2", "Pcyt1a", "Slc51a", "Zdhhc19", "Tfrc"
)

cluster_names <- names(tmp) %>%
  tibble(name = .) %>%
  mutate(index = as.integer(str_extract(name, "\\d+$"))) %>%
  arrange(index) %>%
  pull(name)

heatmap_df <- map_dfr(cluster_names, function(cluster_name) {
  cbind(tmp[[cluster_name]]$control$mean,
        tmp[[cluster_name]]$experiment$mean) %>%
    as.data.frame() %>%
    rownames_to_column("peak") %>%
    mutate(
      gene = str_extract(peak, "[^-]+$"),
      cluster = cluster_name
    ) %>%
    filter(gene %in% genes_del3q29)
}) %>%
  pivot_longer(
    cols = -c(cluster, peak, gene),
    names_to = "sample_ID",
    values_to = "expression"
  ) %>%
  left_join(sample_info, by = "sample_ID")

df_plot <- heatmap_df %>%
  filter(cluster %in% c("cluster_0", "cluster_1")) %>%
  group_by(cluster, peak) %>%
  mutate(scaled_expression = scale(expression)[, 1]) %>%
  ungroup() %>%
  mutate(
    sample_ID = forcats::fct_inorder(sample_ID),
    peak = forcats::fct_rev(peak),
    mouse_genotype = factor(mouse_genotype, levels = c("wtwt", "wtdel"))
  )

# Heatmapa tylko cluster_0 z usunięciem wierszy z samymi NA
df_cluster0 <- df_plot %>%
  filter(cluster == "cluster_0") %>%
  mutate(peak = forcats::fct_rev(factor(peak))) %>%
  group_by(peak) %>%
  filter(!all(is.na(scaled_expression))) %>%
  ungroup()

# Pasek genotype
df_genotype_bar <- df_cluster0 %>%
  distinct(sample_ID, mouse_genotype) %>%
  mutate(peak = "GENOTYPE")

# Pasek treatment
df_treatment_bar <- df_cluster0 %>%
  distinct(sample_ID, treatment) %>%
  mutate(
    treatment = factor(treatment, levels = c("saline", "risperidone")),
    peak = "TREATMENT"
  )

# Łączymy wszystko
df_combined <- bind_rows(df_cluster0, df_genotype_bar, df_treatment_bar)

# Wysokość etykiet nad górnym paskiem
label_y_pos <- max(as.numeric(df_cluster0$peak)) + 1.5

# Numeryczne pozycje próbek
df_labels <- df_cluster0 %>%
  distinct(sample_ID) %>%
  mutate(x = as.numeric(forcats::fct_inorder(sample_ID)))

df_cluster0 <- df_cluster0 %>%
  mutate(x = as.numeric(forcats::fct_inorder(sample_ID)))

df_genotype_bar <- df_genotype_bar %>%
  mutate(x = as.numeric(forcats::fct_inorder(sample_ID)))

# Główna część wykresu
# Zakładam, że masz już przygotowane:
# df_heat  (z kolumnami x, y, peak, scaled_expression)
# df_geno  (z kolumnami x, y = n_peaks+1, mouse_genotype)
# df_lbl   (z kolumnami x, y = n_peaks+2, sample_ID)
# oraz zmienne:
# peak_levels, n_peaks

# 1) nowy porządek kolumn
sample_order <- df_cluster0 %>% 
  distinct(sample_ID, mouse_genotype, treatment) %>%
  arrange(mouse_genotype, treatment) %>%
  pull(sample_ID)

# 2) przygotuj warstwy z x według tego porządku
df_heat <- df_cluster0 %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = as.numeric(factor(peak, levels = peak_levels))
  )

df_geno <- df_genotype_bar %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = n_peaks + 2
  )

df_treat <- df_treatment_bar %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = n_peaks + 1
  )

df_lbl <- df_heat %>%
  distinct(sample_ID, x) %>%
  mutate(
    y = n_peaks + 3
  )




#################################################################################33

# 0) Ustal oryginalne poziomy “peak” i liczbę wierszy heatmapy
peak_levels <- df_cluster0 %>% pull(peak) %>% unique()
n_peaks     <- length(peak_levels)

# 1) Ustal kolejność próbek: najpierw wg genotypu, potem wg leczenia
sample_order <- df_cluster0 %>%
  distinct(sample_ID, mouse_genotype, treatment) %>%
  arrange(mouse_genotype, treatment) %>%
  pull(sample_ID)

# 2) Przygotuj dane z numerycznymi współrzędnymi x/y
df_heat <- df_cluster0 %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = as.numeric(factor(peak, levels = peak_levels))
  )

df_geno <- df_genotype_bar %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = n_peaks + 2
  )

df_treat <- df_treatment_bar %>%
  mutate(
    x = as.numeric(factor(sample_ID, levels = sample_order)),
    y = n_peaks + 1
  )

df_lbl <- df_heat %>%
  distinct(sample_ID, x) %>%
  mutate(
    y = n_peaks + 3
  )

# 3) Rysuj
ggplot() +
  # A) sample_ID u góry
  geom_text(
    data = df_lbl,
    aes(x = x, y = y, label = sample_ID),
    angle = 90, hjust = 0, vjust = 0.5, size = 3
  ) +
  # B) pasek Genotype
  geom_tile(
    data = df_geno,
    aes(x = x, y = y, fill = mouse_genotype),
    height = 0.9
  ) +
  scale_fill_manual(
    values = c("wtwt" = "#1f77b4", "wtdel" = "#ff7f0e"),
    name   = "Genotype"
  ) +
  ggnewscale::new_scale_fill() +
  # C) pasek Treatment
  geom_tile(
    data = df_treat,
    aes(x = x, y = y, fill = treatment),
    height = 0.9
  ) +
  scale_fill_manual(
    values = c("saline"      = "#6baed6",
               "risperidone" = "#fd8d3c"),
    name   = "Treatment"
  ) +
  ggnewscale::new_scale_fill() +
  # D) heatmapa ekspresji
  geom_tile(
    data = df_heat,
    aes(x = x, y = y, fill = scaled_expression),
    color = "grey80"
  ) +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-4, 4),
    na.value = "grey70",
    name     = "Z-score"
  ) +
  # E) osie
  scale_x_continuous(
    breaks   = df_lbl$x,
    labels   = NULL,
    expand   = expansion(mult = c(0.01, 0.01)),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = seq_len(n_peaks + 3),
    labels = c(as.character(peak_levels), "Treatment", "Genotype", ""),
    expand = expansion(add = c(0, 0))
  ) +
  # F) wyłącz przycinanie i ustaw marginesy
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title         = element_blank(),
    axis.text.x.top    = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 8),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 10),
    panel.grid         = element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    plot.margin        = margin(t = 80, r = 10, b = 10, l = 10)
  )
