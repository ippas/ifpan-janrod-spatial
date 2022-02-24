#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

lapply(samples_name, 
       function(x) paste("data/spaceranger-corrected/", 
                         x, 
                         "/outs/spatial/tissue_lowres_image.png",
                         sep = "") %>%
         read.bitmap(.)) -> images_cl


###############################
# Convert the Images to grobs #
###############################
lapply(images_cl, 
       function(x) rasterGrob(x, 
                              width = unit(1, "npc"),
                              height = unit(1, "npc"))) %>%
  tibble(sample=factor(samples_name), grob =.) %>%
  mutate(height = (rapply(images_cl, 
                          function(x) data.frame(height = nrow(x)))),
         width = (rapply(images_cl,
                        function(x) data.frame(weight = ncol(x))))) -> images_tibble


# create table bcs_merge containing information about barcode spots
# for all samples
lapply(samples_name,
       # load data about position tissue
       function(x, y = x) read.csv(
         paste("data/spaceranger-corrected/", 
               x, 
               "/outs/spatial/tissue_positions_list.csv",
               sep = ""),
         col.names=c("barcode","tissue","row","col","imagerow","imagecol"), 
         header = FALSE) %>% 
         
         mutate(
           imagerow = imagerow * (paste("data/spaceranger-corrected/", 
                                        x, 
                                        "/outs/spatial/scalefactors_json.json",
                                        sep = "") %>%
                                    rjson::fromJSON(file = .) %>% 
                                    .$tissue_lowres_scalef),
           imagecol = imagecol * (paste("data/spaceranger-corrected/", 
                                        x, 
                                        "/outs/spatial/scalefactors_json.json",
                                        sep = "") %>%
                                    rjson::fromJSON(file = .) %>% 
                                    .$tissue_lowres_scalef),
           tissue = as.factor(tissue),
           height = images_tibble[images_tibble$sample == x,]$height,
           width = images_tibble[images_tibble$sample == x,]$width
         )
       ) %>% 
  setNames(samples_name) %>% 
  bind_rows(., .id = "sample") -> bcs_merge


##### testing code to refresh how visualize data
resolution = 0.2
data_cluster <- FindClusters(integrated_analysis, resolution = resolution)

DimPlot(data_cluster, reduction = "umap", 
        split.by = "sample", ncol = 4)


data_cluster <- FindClusters(integrated_analysis, resolution = resolution) %>% 
  .$seurat_clusters %>%
  as.data.frame() %>% 
  rename(cluster = ".") %>% 
  rownames_to_column(var = "sample_barcode") %>%
  separate("sample_barcode", c("sample", "barcode"), sep = "_") %>% 
  left_join(., bcs_merge, by = c("barcode", "sample"))

plot_clusters(data_cluster)

colfilt_anno %>%
  filter(gene_name == "Fos") %>%
  select(gene_name, peak_id)


plot_feature(data_cluster = data_cluster,
             peak_id = "merged-samples-peak-94135",
             size = 1)







