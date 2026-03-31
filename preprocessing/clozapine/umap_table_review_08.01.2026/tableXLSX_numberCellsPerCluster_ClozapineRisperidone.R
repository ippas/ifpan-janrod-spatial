

# utwórz tymczasowe środowisko
tmp_env <- new.env(parent = emptyenv())

# wczytaj RData do tego środowiska
load("results/risperidone/risperidone-half.RData", envir = tmp_env)

# zobacz co się wczytało
ls(tmp_env)

## ===============================
## RISPERIDONE – res 0.4
## ===============================
ris_res0.4 <- tmp_env$risperidone_st_data_half$clusters %>% 
  select(sample, barcode, cluster_resolution_0.4) %>% 
  mutate(cluster_resolution_0.4 = paste0("cluster_", cluster_resolution_0.4)) %>% 
  group_by(sample, cluster_resolution_0.4) %>% 
  summarise(n_spots = n(), .groups = "drop") %>% 
  pivot_wider(
    names_from  = sample,
    values_from = n_spots,
    values_fill = 0
  ) %>% 
  mutate(cluster_id = as.integer(str_remove(cluster_resolution_0.4, "cluster_"))) %>% 
  arrange(cluster_id) %>% 
  select(-cluster_id) %>% 
  select("cluster_resolution_0.4", 
         "S6269Nr1",  "S6269Nr3",  "S7788Nr1",  "S7788Nr2",  "S7788Nr3",  "S7788Nr11",
         "S6230Nr3",  "S6230Nr4",  "S6269Nr2",  "S6269Nr4",  "S7788Nr15", "S7788Nr16")

## ===============================
## CLOZAPINE – res 0.4
## ===============================
clo_res0.4 <- clozapine_st_data_half$clusters %>% 
  select(sample, barcode, cluster_resolution_0.4) %>% 
  mutate(cluster_resolution_0.4 = paste0("cluster_", cluster_resolution_0.4)) %>% 
  group_by(sample, cluster_resolution_0.4) %>% 
  summarise(n_spots = n(), .groups = "drop") %>% 
  pivot_wider(
    names_from  = sample,
    values_from = n_spots,
    values_fill = 0
  ) %>% 
  mutate(cluster_id = as.integer(str_remove(cluster_resolution_0.4, "cluster_"))) %>% 
  arrange(cluster_id) %>% 
  select(-cluster_id) %>% 
  select("cluster_resolution_0.4", 
         "S6269Nr1",  "S6269Nr3",  "S7788Nr1",  "S7788Nr2",  "S7788Nr3",  "S7788Nr11",
         "S6230Nr1",  "S6230Nr2",  "S7788Nr5",  "S7788Nr8",  "S7788Nr12", "S7788Nr13")


## ===============================
## CLOZAPINE – res 0.45
## ===============================
clo_res0.45 <- clozapine_st_data_half$clusters %>% 
  select(sample, barcode, cluster_resolution_0.45) %>% 
  mutate(cluster_resolution_0.45 = paste0("cluster_", cluster_resolution_0.45)) %>% 
  group_by(sample, cluster_resolution_0.45) %>% 
  summarise(n_spots = n(), .groups = "drop") %>% 
  pivot_wider(
    names_from  = sample,
    values_from = n_spots,
    values_fill = 0
  ) %>% 
  mutate(cluster_id = as.integer(str_remove(cluster_resolution_0.45, "cluster_"))) %>% 
  arrange(cluster_id) %>% 
  select(-cluster_id) %>% 
  select("cluster_resolution_0.45", 
         "S6269Nr1",  "S6269Nr3",  "S7788Nr1",  "S7788Nr2",  "S7788Nr3",  "S7788Nr11",
         "S6230Nr1",  "S6230Nr2",  "S7788Nr5",  "S7788Nr8",  "S7788Nr12", "S7788Nr13")




out_file <- "/home/mateusz/projects/ifpan-janrod-spatial/results/clozapine/umap_table_review_08.01.2026/numberSpots_perClusterSample_RisperidoneClozapine_09.01.2026.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "ris_res0.4")
addWorksheet(wb, "clo_res0.4")
addWorksheet(wb, "clo_res0.45")

writeData(wb, "ris_res0.4", ris_res0.4)
writeData(wb, "clo_res0.4", clo_res0.4)
writeData(wb, "clo_res0.45", clo_res0.45)

saveWorkbook(wb, out_file, overwrite = TRUE)


rm(
  tmp_env,
  ris_res0.4,
  clo_res0.4,
  clo_res0.45,
  wb,
  out_file
)
gc()

