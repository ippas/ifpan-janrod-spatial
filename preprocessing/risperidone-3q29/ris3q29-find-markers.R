tmp <- FindClusters(ris3q29_integrate, resolution = 0.4, verbose = TRUE)


marker_list <- list()

for (clust in c(0:19)) {
  message(paste("Processing cluster", clust))
  marker_list[[as.character(clust)]] <- FindMarkers(tmp, ident.1 = clust)
}


marker_list$`0` %>% 
  arrange((avg_log2FC)) %>% 
  head(20)


all_markers <- FindAllMarkers(
  tmp,
  only.pos = TRUE,         # tylko nadekspresjonowane geny
  min.pct = 0.25,          # co najmniej 25% komórek w klastrze musi wyrażać gen
  logfc.threshold = 0.25   # tylko geny z FC > 0.25 (log2)
)

all_markers <- FindAllMarkers(
  tmp,
  # only.pos = TRUE,         # tylko nadekspresjonowane geny
  min.pct = 0.25,          # co najmniej 25% komórek w klastrze musi wyrażać gen
  logfc.threshold = 0.25   # tylko geny z FC > 0.25 (log2)
)

all_markers %>% 
  rename(gene="peak") %>% 
  mutate(gene = str_extract(peak, "[^-]+$")) %>% 
  { rownames(.) <- NULL; . } %>% 
  select(peak, gene, everything()) %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  group_by(cluster) %>%
  group_split(.keep = TRUE) %>%
  set_names(all_markers %>% mutate(cluster = paste0("cluster_", cluster)) %>% group_by(cluster) %>% group_keys() %>% pull(cluster)) %>% 
  writexl::write_xlsx(path = "/home/mateusz/projects/ifpan-janrod-spatial/results/risperidone-3q29/ris3q29-allSamples-markers-log2FC0.25-pct0.25.xlsx")


  

all_markers %>% 
  rename(gene="peak") %>% 
  mutate(gene = str_extract(peak, "[^-]+$")) %>% 
  { rownames(.) <- NULL; . } %>% 
  select(peak, gene, everything()) %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  mutate(gene = toupper(gene)) %>% 
  filter(gene %in%  c(
    "FKBP5", "PDK4", "ZBTB16", "SAA1", "TSC22D3", "FAM107A", "HIF3A", "LCN2", "GLUL", "PNMT",
    "CD163", "PER1", "MAOA", "ERRFI1", "SERPINE1", "KLF9", "TXNIP", "PHACTR3", "DUSP1", "ANGPTL4",
    "FLRT3", "PRG4", "CDKN1A", "MAP3K6", "ANKRD1", "CNR1", "HP", "CPM", "PTK2B", "CRYAB",
    "CCL20", "LEP", "RGS2", "AKR1B10", "SULT1E1", "HSPA1B", "KCNJ11", "APOD", "CRISPLD2", "IL1R2",
    "SLPI", "ADRA1B", "TFCP2L1", "LOX", "CXCL5", "FMO2", "OBP2A", "MMP25", "ART3", "FAM90A7",
    "FXYD4", "DCN", "SEC14L4", "GPX3", "IL17RA", "GPM6B", "VN1R5"
  ))
  
  