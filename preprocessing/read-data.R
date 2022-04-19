#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

meta_data <- read.table("data/samples-spatial-metadata.tsv",
                        header = TRUE,
                        sep = "\t") %>% filter(experiment == "tif-ldopa")

# Require source files
info_peaks <- read.table("data/ldopa/gene-annotation/peaks-annotate-reduction.tsv", 
                         header = TRUE,
                         sep = "\t")
                         # col.names = c('chr_peak',
                         #               'start_peak',
                         #               'end_peak',
                         #               'peak_id',
                         #               'score_int(-10*log10pvalue)',
                         #               'strand_coverage',
                         #               'fold_change_peak_summit',
                         #               '-log10pvalue_peak_summit',
                         #               '-log10qvalue_peak_summit',
                         #               'relative_summit_position_peak_start',
                         #               'type_peak',
                         #               'chr_gene',
                         #               'start_gene',
                         #               'end_gene',
                         #               'gene_id',
                         #               'gene_name',
                         #               'strand_gene')) 
                         # 

###############
# read images #
###############
lapply(samples_name, 
       function(x) paste("data/ldopa/spaceranger-corrected/", 
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
         paste("data/ldopa/spaceranger-corrected/", 
               x, 
               "/outs/spatial/tissue_positions_list.csv",
               sep = ""),
         col.names=c("barcode","tissue","row","col","imagerow","imagecol"), 
         header = FALSE) %>% 
         
         mutate(
           imagerow = imagerow * (paste("data/ldopa/spaceranger-corrected/", 
                                        x, 
                                        "/outs/spatial/scalefactors_json.json",
                                        sep = "") %>%
                                    rjson::fromJSON(file = .) %>% 
                                    .$tissue_lowres_scalef),
           imagecol = imagecol * (paste("data/ldopa/spaceranger-corrected/", 
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
