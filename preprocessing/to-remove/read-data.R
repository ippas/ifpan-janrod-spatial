#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

meta_data <- read.table("data/metadata-antipsychotics.tsv",
                        header = TRUE,
                        sep = "\t") %>% filter(treatment %in% c("risperidone", "saline"))

# Require source files
info_peaks <- read.table("data/risperidone/gene-annotation/peaks-annotate-reduction.tsv", 
                         header = TRUE,
                         sep = "\t")
 

###############
# read images #
###############
lapply(samples_name, 
       function(x) paste(path_to_data, 
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
         paste(path_to_data, 
               x, 
               "/outs/spatial/tissue_positions_list.csv",
               sep = ""),
         col.names=c("barcode","tissue","row","col","imagerow","imagecol"), 
         header = FALSE) %>% 
         
         mutate(
           imagerow = imagerow * (paste(path_to_data, 
                                        x, 
                                        "/outs/spatial/scalefactors_json.json",
                                        sep = "") %>%
                                    rjson::fromJSON(file = .) %>% 
                                    .$tissue_lowres_scalef),
           imagecol = imagecol * (paste(path_to_data, 
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

