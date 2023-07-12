# read function
source("preprocessing/functions/functions-spatial-data.R")
source("preprocessing/functions/statistics-functions.R")
source("preprocessing/functions/visualization-functions.R")
source("preprocessing/functions/umi-per-spot.R")

# executes seurat analysis for risperidone
path_to_data <- "data/risperidone/spaceranger-corrected-half/"
nfeatures <- 2000
dims <- 1:30

# create sample_names vector contain id samples to analysis
sample_names <- list.files(path = path_to_data)

# read metadata for risperidone
metadata_risperidone <- read_metadata(file_path = "data/metadata-antipsychotics.tsv", 
                                      treatments = c("saline", "risperidone"))

# visualization test
samples_saline <- metadata_risperidone %>% 
  filter(treatment == "saline" & mouse_genotype == "wt") %>%
  .[, 1]

samples_risperidone <- metadata_risperidone %>% 
  filter(treatment == "risperidone" & mouse_genotype == "wt") %>%
  .[, 1]

# read peaks for risperidone
info_peaks_risperidone <- read_gene_annotation(file_path = "data/risperidone/gene-annotation/peaks-annotate-reduction.tsv")

# executes seurat analysis for risperidone
risperidone_integrate_half <- integrate_data_seurat(path_to_data, nfeatures, dims)

# read images to spatial transcriptoms for risperidone
images_risperidone_half <- create_images_tibble(path_to_data, sample_names)

# read barcode data for risperidone
barcode_risperidone_half <- create_barcode_data(
  path_to_data = path_to_data,
  sample_names = sample_names,
  images_tibble = images_risperidone_half
)

# poprawić do 
risperidone_st_data_half <- create_spatial_data(sample_names = sample_names,
                                           metadata = metadata_risperidone,
                                           barcode_info = barcode_risperidone_half,
                                           images_info = images_risperidone_half,
                                           integrated_data = risperidone_integrate_half,
                                           peaks_info = info_peaks_risperidone)

risperidone_st_data_half <- add_seurat_data(spatial_data = risperidone_st_data_half,
                                       integrated_data = risperidone_integrate_half)

risperidone_st_data_half <- add_filtered_data(spatial_data = risperidone_st_data_half,
                                         mean_expression_threshold = 0.5)

risperidone_st_data_half <- add_colfilt_data(spatial_data = risperidone_st_data_half,
                                        min_spot_threshold = 0,
                                        expression_threshold = 2)

risperidone_st_data_half <- add_range_normalize_data(spatial_data = risperidone_st_data_half,
                                                range = 1500,
                                                flatten = 1,
                                                threshold = 500)

risperidone_st_data_half <- add_clusters_data(spatial_data = risperidone_st_data_half,
                                         integrated_data = risperidone_integrate_half,
                                         resolution_start = 0.05,
                                         resolution_end = 2,
                                         resolution_step = 0.05)

risperidone_st_data_half <- evaluate_clustering_stability(spatial_data = risperidone_st_data_half,
                                                     seurat_object = risperidone_integrate_half,
                                                     resolution_start = 0.05,
                                                     resolution_end = 2,
                                                     resolution_step = 0.05)




for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)){
  risperidone_st_data_half <-
    add_quantile_norm_data(
      spatial_data = risperidone_st_data_half,
      resolution = {{resolution}},
      num_cores = 24,
      data_type = "raw_data"
    )
}


risperidone_st_data_half$stability_results %>%
  ggplot(aes(x = resolution, y = silhouette_score)) +
  geom_line()



data.frame(data_type = c(
  rep("raw_data", 3),
  rep("range_normalize", 3),
  c(
    "quantile_normalize_resolution_0.1",
    "quantile_normalize_resolution_0.4",
    "quantile_normalize_resolution_0.8"
  )
),
resolution = rep(c(0.1, 0.4, 0.8), 3),
data_type_name = c(
  rep("raw_data", 3),
  rep("range_normalize", 3),
  rep("quantile_normalize", 3)
)) -> data_params_df


data.frame(data_type = c(
  rep("raw_data", 6),
  rep("raw_data", 6),
  rep("range_normalize", 6),
  rep("seurat", 6),
  c(
    "quantile_normalize_resolution_0.05",
    "quantile_normalize_resolution_0.1",
    "quantile_normalize_resolution_0.15",
    "quantile_normalize_resolution_0.2",
    "quantile_normalize_resolution_0.4",
    "quantile_normalize_resolution_0.8"
  )
),
resolution = rep(c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8), 5),
data_type_name = c(
  rep("raw_data", 6),
  rep("quantile_metric", 6),
  rep("range_normalize", 6),
  rep("seurat", 6),
  rep("quantile_normalize", 6)
)) %>%
  mutate(quantile_normalization = ifelse(data_type_name == "quantile_metric", TRUE, FALSE)) -> data_params_df


summarize_and_test(spatial_data = risperidone_st_data_half,
                   trim = 0.05, 
                   num_cores = 24,
                   data_params_df = data_params_df,
                   control_samples = samples_saline,
                   experiment_samples = samples_risperidone,
                   mean_threshold = 0,
                   metrics = c("mean", "median", "skewness", "kurtosis")) -> risperidone_summary_statistics_half

save(samples_saline,
     samples_risperidone,
     risperidone_integrate_half,
     risperidone_st_data_half,
     risperidone_summary_statistics_half, 
     file = "results/risperidone/risperidone-half.RData")



# DimPlot(risperidone_integrate, reduction = "umap", 
#         split.by = "sample", ncol = 4)
# 
# DimPlot(risperidone_integrate, reduction = "umap", 
#         ncol = 4)
# 
# visualize clusters
# spatial_cluster(spatial_data = risperidone_st_data_half,
#                 resolution = 0.8,
#                 samples = c(samples_saline, samples_risperidone),
#                 palette = palette_allen, 
#                 size= 1.0, 
#                 ncol = 4)
# 
# spatial_interest_cluster(cluster = 3,
#                          # seurat_object = integrated_analysis,
#                          spatial_data = risperidone_st_data_half,
#                          resolution = 0.8,
#                          samples = c(samples_saline, samples_risperidone),
#                          size= 1,
#                          ncol = 4)
# 
# 
# spatial_gene_plot(spatial_data = risperidone_st_data_half,
#                   type_data = "quantile_normalize",
#                   gene = "Itpk1",
#                   samples =  c(samples_saline, samples_risperidone),
#                   min_percentile = 0.00,
#                   max_percentile = 1,
#                   size = 0.8,
#                   ncol = 4,
#                   normalization = T)
# 
# spatial_gene_plot(spatial_data = risperidone_st_data_half,
#                   type_data = "raw_data",
#                   gene = "Egr4",
#                   samples =  c(samples_saline, samples_risperidone),
#                   min_percentile = 0.00,
#                   max_percentile = 1,
#                   size = 0.8,
#                   ncol = 4,
#                   normalization = T)


# calculate_gene_expression_stats(spatial_data = risperidone_st_data,
#                                 stat_test = "t.test",
#                                 data_type = "raw_data",
#                                 expression_unit = "sample",
#                                 control_samples = samples_saline,
#                                 experiment_samples = samples_risperidone,
#                                 save_results = F
# ) -> stat_raw_risperidone_half



genes2mind_ris <- c(
  "Rps13",
  "Trove2",
  "Nrp1",
  "Prrt3",
  "Hes5",
  "Pkp2",
  "Ndufaf4",
  "Zfp866",
  "Dok3",
  "Dusp11",
  "Lrrtm1",
  "8430408G22Rik",
  "Cttnbp2",
  "Fos",
  "1810013L24Rik",
  "Tnfrsf25",
  "Polr3e",
  "2410131K14Rik",
  "Dclk1",
  "Cttnbp2",
  "Srsf5",
  "Clk1",
  "Ywhag",
  "Egr2",
  "Csrnp1",
  "Cd200",
  "Mest",
  "Plat",
  "Srsf6",
  "BC031353",
  "Arc",
  "Cebpb",
  "Angptl4",
  "Txnip",
  "Slco4a1",
  "Tra2a",
  "Sox18",
  "Dusp12",
  "Btrc",
  "Kif1b",
  "Iqcb1",
  "Tiparp",
  "Edn1",
  "Ddit4l",
  "Mrpl15",
  "Rhou",
  "Baiap2",
  "Cdkn1a",
  "Tsc22d3",
  "Pex13",
  "Spsb1",
  "Cfp",
  "Plin4",
  "Sgk1",
  "Arrdc2",
  "Eif2a",
  "Ccdc88c",
  "Epm2a",
  "Junb",
  "Prmt6",
  "Nfkbia",
  "Cirbp",
  "1200016B10Rik",
  "Slc38a5",
  "Sult1a1",
  "Cnp",
  "Clk4",
  "Ddit4",
  "Errfi1",
  "1190005F20Rik",
  "Cul3",
  "Dclk1",
  "Fosb",
  "Slc2a1",
  "Itgad",
  "Arrdc2",
  "Srfbp1",
  "H3f3b",
  "Dnase1l2",
  "Prg4",
  "Mett11d1",
  "Pmch",
  "Dusp1",
  "Arrdc3",
  "Cdkn1a",
  "Npas4",
  "Hcrt",
  "9830001H06Rik",
  "Grpel1",
  "Gucy1a3",
  "Zmym1",
  "Snrnp48",
  "Hbb-b2",
  "Eprs",
  "Dlst",
  "Fam120b",
  "Cnot3"
)
