library(tidyverse)
library(patchwork)

# ===== 1. Wczytanie danych =====

file_path <- "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone/spatial-cluster-stability/risperidone_half_SDMbench_original_metrics.tsv"

df <- read_tsv(file_path)

# ===== 2. Obliczenie zmian między kolejnymi resolution =====

df_delta <- df %>%
  arrange(resolution) %>%
  mutate(
    delta_clusters = n_clusters - lag(n_clusters),
    delta_CHAOS = SDMBench_CHAOS - lag(SDMBench_CHAOS),
    delta_PAS = SDMBench_PAS - lag(SDMBench_PAS)
  ) %>%
  filter(!is.na(delta_PAS))

# ===== 3. Wykresy surowych wartości =====

p_clusters <- ggplot(df, aes(x = resolution, y = n_clusters)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "Number of clusters",
    x = "Resolution",
    y = "Number of clusters"
  ) +
  theme_classic()

p_chaos <- ggplot(df, aes(x = resolution, y = SDMBench_CHAOS)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "CHAOS",
    x = "Resolution",
    y = "CHAOS"
  ) +
  theme_classic()

p_pas <- ggplot(df, aes(x = resolution, y = SDMBench_PAS)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "PAS",
    x = "Resolution",
    y = "PAS"
  ) +
  theme_classic()

raw_plot <- p_clusters / p_chaos / p_pas

raw_plot


# ===== 4. Wykresy zmian / "pochodnych" =====

p_delta_clusters <- ggplot(df_delta, aes(x = resolution, y = delta_clusters)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "Change in number of clusters",
    x = "Resolution",
    y = expression(Delta~clusters)
  ) +
  theme_classic()

p_delta_chaos <- ggplot(df_delta, aes(x = resolution, y = delta_CHAOS)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "Change in CHAOS",
    x = "Resolution",
    y = expression(Delta~CHAOS)
  ) +
  theme_classic()

p_delta_pas <- ggplot(df_delta, aes(x = resolution, y = delta_PAS)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(
    title = "Change in PAS",
    x = "Resolution",
    y = expression(Delta~PAS)
  ) +
  theme_classic()

delta_plot <- p_delta_clusters / p_delta_chaos / p_delta_pas

delta_plot
