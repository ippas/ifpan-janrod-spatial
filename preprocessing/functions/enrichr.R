SpatialDimPlot(risperidone_integrate_half, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)

tmp <- FindClusters(risperidone_integrate_half, resolution = 0.4)

markers <- FindMarkers(tmp, ident.1 = 3)

markers %>% 
  # head %>%
  rownames_to_column(var = "peak_id") %>%
  left_join(., 
            risperidone_st_data_half$raw_data$annotate[c("peak_id", "gene_name")], by = "peak_id") %>% 
  pull(gene_name) -> gene_vector1


perform_enrichment_analysis(genes = gene_vector1, database = "Mouse_Gene_Atlas", top_rows = 5) %>% head(10)


spatial_interest_cluster(cluster = 13,
                         # seurat_object = integrated_analysis,
                         spatial_data = risperidone_st_data_half,
                         resolution = 0.4,
                         samples = c(samples_saline, samples_risperidone),
                         size= 1,
                         ncol = 4)

risperidone_st_data_half$raw_data$annotate[c("peak_id", "gene_name")]


databases_enrichr <- listEnrichrDbs()
databases_enrichr %>% filter(libraryName == "Mouse_Gene_Atlas")


# Define the function
perform_enrichment_analysis <- function(genes, database, top_rows) {
  # Perform the enrichment analysis
  results <- enrichR::enrichr(genes, database = database)
  
  # If the analysis was successful, print the top rows
  if (!is.null(results)) {
    print(head(results, n = top_rows))
  } else {
    print("Enrichment analysis failed. Please check your inputs.")
  }
}


perform_enrichment_analysis(genes = gene_vector1, database = "Mouse_Gene_Atlas", top_rows = 5)


spatial_interest_cluster(cluster = 5,
                         # seurat_object = integrated_analysis,
                         spatial_data = risperidone_st_data_half,
                         resolution = 0.4,
                         samples = c(samples_saline, samples_risperidone),
                         size= 1,
                         ncol = 4)
