cloz_env <- new.env()
load("results/clozapine/clozapine-half.RData", envir = cloz_env)

message("Loaded into separate environment: cloz_env")

cloz_env$clozapine_summary_statistics_half$quantile_normalize$resolution_0.4$cluster_0$control$mean
cloz_env$clozapine_summary_statistics_half$quantile_normalize$resolution_0.4$cluster_0$experiment$mean

data.frame(
  peak = cloz_env$clozapine_summary_statistics_half$quantile_normalize$resolution_0.4$cluster_0$peak,
  gene = cloz_env$clozapine_summary_statistics_half$quantile_normalize$resolution_0.4$cluster_0$gene
) %>% 
  filter(gene %in% c(
    "Bhlhe40",
    "Cacna1i",
    "Egr3",
    "Homer1",
    "Kalrn",
    "Olig2",
    "Phactr3",
    "Stum"
  )) -> supplementGenes_26.11.2025

cloz_env$clozapine_summary_statistics_half$quantile_normalize$resolution_0.4 %>% 
  lapply(., function(x){
    
    cbind(
      x$control$mean,
      x$experiment$mean 
    ) %>% 
      .[supplementGenes_26.11.2025$peak,] %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "peak") %>% 
      left_join(., supplementGenes_26.11.2025, by = "peak")
    
  }) %>% 
  bind_rows(.id = "cluster") %>% 
  select(c(peak, gene, cluster, everything())) %>% 
  filter(
    (gene == "Egr3"     & cluster %in% c("cluster_13", "cluster_5")) |
      (gene == "Bhlhe40"  & cluster %in% c("cluster_11")) |
      (gene == "Cacna1i"  & cluster %in% c("cluster_0", "cluster_11")) |
      (gene == "Homer1"   & cluster %in% c("cluster_8")) |
      (gene == "Kalrn"    & cluster %in% c("cluster_11")) |
      (gene == "Olig2"    & cluster %in% c("cluster_10")) |
      (gene == "Phactr3"  & cluster %in% c("cluster_9")) |
      (gene == "Stum"     & cluster %in% c("cluster_8"))
  ) %>% 
  write.table(
    x           = .,          # Twój data.frame
    file        = "results/clozapine/supplementGenesClozapine_26.11.2025.tsv",
    sep         = "\t",        # TAB = TSV
    quote       = FALSE,       # bez cudzysłowów
    row.names   = FALSE        # bez nazw wierszy
  )




cloz <- cloz_env$clozapine_integrate_half
cloz <- FindNeighbors(cloz, reduction = "pca", dims = 1:30)
cloz <- FindClusters(cloz, resolution = 0.4)

svg("results/clozapine/supplemntary_results_paper_03.12.2025/umap_clusters_clozapine_03.12.2025.svg", width = 7, height = 7)
DimPlot(cloz, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
dev.off()
