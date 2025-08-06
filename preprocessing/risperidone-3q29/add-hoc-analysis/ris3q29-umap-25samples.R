install.packages("ggrepel")

# ##############################################################################
# ---- load data ----
# ##############################################################################
load("results/risperidone-3q29/ris3q29_integrate.RData")


ris3q29_integrate <- FindClusters(ris3q29_integrate, resolution = 0.4)

DimPlot(ris3q29_integrate, reduction = "umap", label = TRUE)

spatial_obj <- RunUMAP(spatial_obj, dims = 1:30)

rm(ris3q29_integrate)


Embeddings(ris3q29_integrate, reduction = "umap") %>%
  as.data.frame()


# Wyciągnięcie koordynatów UMAP
umap_df <- Embeddings(ris3q29_integrate, reduction = "umap") %>%
  as.data.frame()

# Dodanie metadanych (wszystkie metadane)
umap_df <- umap_df %>%
  rownames_to_column(var = "spot") %>%
  left_join(
    ris3q29_integrate@meta.data %>% rownames_to_column(var = "spot"),
    by = "spot"
  )

umap_df$seurat_clusters <- as.character(umap_df$seurat_clusters)

umap_df$color <- palette_allen[umap_df$seurat_clusters]

centroids <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2)
  )


# Pobierz unikalne wartości i posortuj numerycznie
levels_num <- sort(as.numeric(unique(umap_df$seurat_clusters)))
levels_num <- as.character(levels_num)  # spowrotem w char

umap_df$seurat_clusters <- factor(
  umap_df$seurat_clusters,
  levels = levels_num
)

filename <- paste0(
  "results/risperidone-3q29/figures/umap/UMAP_ris3q29_25samples_res0.4_", 
  "25.06.2025", 
  ".svg"
)

svg(
  filename = filename,
  width = 10,
  height = 8
)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = palette_allen) +
  theme_classic() +
  labs(
    color = "clusters",
    title = "UMAP, ris3q29, 25 samples, res0.4",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 4)
    )
  ) +
  geom_text(
    data = centroids,
    aes(label = seurat_clusters),
    color = "black",
    size = 6,
    fontface = "bold"
  )

dev.off()


