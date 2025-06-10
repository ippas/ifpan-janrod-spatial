# 0) Przygotowanie surowych danych
tmp <- wtDel_summary_statistics$quantile_normalize$resolution_0.4

genes_del3q29 <- c(
  "Bdh1", "Dlg1", "Mfi2", "Pigz", "Ncbp2", "Senp5", "Pak2", "Pigx",
  "Cep19", "Nrros", "Bex6", "Fbxo45", "Wdr53", "Smco1", "Rnf168",
  "Ubxn7", "Tm4sf19", "Tctex1d2", "Pcyt1a", "Slc51a", "Zdhhc19", "Tfrc"
)

cluster_names <- names(tmp) %>%
  tibble::tibble(name = .) %>%
  dplyr::mutate(index = as.integer(stringr::str_extract(name, "\\d+$"))) %>%
  dplyr::arrange(index) %>%
  dplyr::pull(name)

heatmap_df <- purrr::map_dfr(cluster_names, function(cl) {
  cbind(
    tmp[[cl]]$control$mean,
    tmp[[cl]]$experiment$mean
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(
      gene    = stringr::str_extract(peak, "[^-]+$"),
      cluster = cl
    ) %>%
    dplyr::filter(gene %in% genes_del3q29)
}) %>%
  tidyr::pivot_longer(
    cols      = -c(cluster, peak, gene),
    names_to  = "sample_ID",
    values_to = "expression"
  ) %>%
  dplyr::left_join(sample_info, by = "sample_ID")

# 1) Skalowanie i filtrowanie cluster_0
df_plot <- heatmap_df %>%
  dplyr::filter(cluster == "cluster_0") %>%
  dplyr::group_by(cluster, peak) %>%
  dplyr::mutate(scaled_expression = scale(expression)[,1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sample_ID      = forcats::fct_inorder(sample_ID),
    peak           = forcats::fct_rev(peak),
    mouse_genotype = factor(mouse_genotype, levels = c("wtwt", "wtdel"))
  )

df_cluster0 <- df_plot %>%
  dplyr::group_by(peak) %>%
  dplyr::filter(!all(is.na(scaled_expression))) %>%
  dplyr::ungroup()

df_genotype_bar <- df_cluster0 %>%
  dplyr::distinct(sample_ID, mouse_genotype) %>%
  dplyr::mutate(peak = "GENOTYPE")

df_treatment_bar <- df_cluster0 %>%
  dplyr::distinct(sample_ID, treatment) %>%
  dplyr::mutate(
    treatment = factor(treatment, levels = c("saline", "risperidone")),
    peak      = "TREATMENT"
  )

# 2) Wektory poziomów i porządek próbek
peak_levels  <- c("GENOTYPE", "TREATMENT", levels(df_cluster0$peak))
sample_order <- df_cluster0 %>%
  dplyr::distinct(sample_ID, mouse_genotype, treatment) %>%
  dplyr::arrange(mouse_genotype, treatment) %>%
  dplyr::pull(sample_ID)
y_levels     <- peak_levels %>% factor(levels = peak_levels)

# 3) Przypisanie współrzędnych
df_lbl <- df_cluster0 %>%
  dplyr::distinct(sample_ID) %>%
  dplyr::mutate(
    x = factor(sample_ID, levels = sample_order),
    y = factor("SAMPLE", levels = levels(y_levels))
  )

df_geno <- df_genotype_bar %>%
  dplyr::mutate(
    x = factor(sample_ID, levels = sample_order),
    y = factor("GENOTYPE", levels = levels(y_levels))
  )

df_treat <- df_treatment_bar %>%
  dplyr::mutate(
    x = factor(sample_ID, levels = sample_order),
    y = factor("TREATMENT", levels = levels(y_levels))
  )

df_heat <- df_cluster0 %>%
  dplyr::mutate(
    x = factor(sample_ID, levels = sample_order),
    y = peak  # already factor with levels(df_cluster0$peak)
  )

# 4) Rysowanie
ggplot2::ggplot() +
  # A) sample_ID
  ggplot2::geom_text(
    data = df_lbl,
    ggplot2::aes(x = x, y = y, label = sample_ID),
    angle = 90, hjust = 0, vjust = 0.5, size = 3
  ) +
  # B) Genotype bar
  ggplot2::geom_tile(
    data = df_geno,
    ggplot2::aes(x = x, y = y, fill = mouse_genotype),
    height = 0.9
  ) +
  ggplot2::scale_fill_manual(
    name   = "Genotype",
    values = c(wtwt = "#1f77b4", wtdel = "#ff7f0e")
  ) +
  ggnewscale::new_scale_fill() +
  # C) Treatment bar
  ggplot2::geom_tile(
    data = df_treat,
    ggplot2::aes(x = x, y = y, fill = treatment),
    height = 0.9
  ) +
  ggplot2::scale_fill_manual(
    name   = "Treatment",
    values = c(saline = "#6baed6", risperidone = "#fd8d3c")
  ) +
  ggnewscale::new_scale_fill() +
  # D) Heatmapa ekspresji
  ggplot2::geom_tile(
    data = df_heat,
    ggplot2::aes(x = x, y = y, fill = scaled_expression),
    color = "grey80"
  ) +
  ggplot2::scale_fill_gradient2(
    name     = "Z-score",
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-4, 4),
    na.value = "grey70"
  ) +
  # E) Skale dyskretne
  ggplot2::scale_x_discrete(position = "top", expand = ggplot2::expansion(mult = c(0.01,0.01))) +
  ggplot2::scale_y_discrete(limits = c("SAMPLE", peak_levels)) +
  # F) Marginesy i clip
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.title         = ggplot2::element_blank(),
    axis.text.x.top    = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, size = 8),
    axis.text.x.bottom = ggplot2::element_blank(),
    axis.ticks.x       = ggplot2::element_blank(),
    axis.text.y        = ggplot2::element_text(size = 10),
    panel.grid         = ggplot2::element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    plot.margin        = ggplot2::margin(t = 80, r = 10, b = 10, l = 10)
  )
