
cluster_names <- names(
  salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4
)

cluster_mean_expression_list <- list()
cluster_mean_long_df <- list()

for (cluster_name in cluster_names) {
  
  salWt_df <- salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$control$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  salDel_df <- salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$experiment$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  risWt_df <- risWtRisDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$control$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  risDel_df <- risWtRisDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$experiment$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  
  mean_df <- tibble(
    name_id = salWt_df$name_id,
    mean_salWt = rowMeans(salWt_df %>% select(-name_id), na.rm = TRUE),
    mean_salDel = rowMeans(salDel_df %>% select(-name_id), na.rm = TRUE),
    mean_risWt = rowMeans(risWt_df %>% select(-name_id), na.rm = TRUE),
    mean_risDel = rowMeans(risDel_df %>% select(-name_id), na.rm = TRUE)
  )
  
  
  expression_df <- bind_cols(
    salWt_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("salWt_", .x)),
    
    salDel_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("salDel_", .x)),
    
    risWt_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("risWt_", .x)),
    
    risDel_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("risDel_", .x))
  )
  
  
  cluster_mean_expression_list[[cluster_name]] <- bind_cols(
    mean_df,
    expression_df
  )
  
  
  cluster_mean_long_df[[cluster_name]] <- mean_df %>%
    mutate(cluster_id = cluster_name) %>%
    select(
      name_id,
      cluster_id,
      mean_salWt,
      mean_salDel,
      mean_risWt,
      mean_risDel
    )
}

cluster_mean_long_df <- bind_rows(cluster_mean_long_df)


library(readxl)

edger_df <- read_excel(
  "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_09.09.2025.xlsx"
)

edger_df 

edger_df_with_means <- edger_df %>%
  left_join(
    cluster_mean_long_df %>%
      select(
        name_id,
        cluster_id,
        mean_salWt,
        mean_salDel,
        mean_risWt,
        mean_risDel
      ),
    by = c(
      "name_id" = "name_id",
      "cluster" = "cluster_id"
    )
  )


edger_df_with_means

library(writexl)

write_xlsx(
  edger_df_with_means,
  "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_withMeans_09.06.2026.xlsx"
)

meanExpressionPerSampleMin20spots_quantileNormalizeData_df <- bind_rows(
  cluster_mean_expression_list,
  .id = "cluster"
)

meanExpressionPerSampleMin20spots_quantileNormalizeData_df %>%
  head()

write_xlsx(
  meanExpressionPerSampleMin20spots_quantileNormalizeData_df,
  "results/risperidone-3q29/edger_statistics/meanExpressionPerSampleMin20spots_quantileNormalizeData_09.06.2026.xlsx"
)

rm(
  cluster_names,
  cluster_name,
  salWt_df,
  salDel_df,
  risWt_df,
  risDel_df,
  mean_df,
  expression_df,
  cluster_mean_expression_list,
  cluster_mean_long_df,
  edger_df,
  edger_df_with_means,
  meanExpressionPerSampleMin20spots_quantileNormalizeData_df
)

################################################################################
library(tidyverse)
library(readxl)
library(writexl)

load("results/risperidone-3q29/salWtSalDel_summary_statistics.RData")
load("results/risperidone-3q29/risWtRisDel_summary_statistics.RData")

cluster_names <- names(
  salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4
)

cluster_mean_expression_list <- list()
cluster_mean_long_df <- list()

for (cluster_name in cluster_names) {
  
  salWt_df <- salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$control$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  salDel_df <- salWtSalDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$experiment$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  risWt_df <- risWtRisDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$control$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  risDel_df <- risWtRisDel_summary_statistics$quantile_normalize$resolution_0.4[[cluster_name]]$experiment$mean %>%
    as.data.frame() %>%
    rownames_to_column("name_id")
  
  
  mean_df <- tibble(
    name_id = salWt_df$name_id,
    
    mean_salWt = rowMeans(
      salWt_df %>% select(-name_id),
      na.rm = TRUE
    ),
    
    mean_salDel = rowMeans(
      salDel_df %>% select(-name_id),
      na.rm = TRUE
    ),
    
    mean_risWt = rowMeans(
      risWt_df %>% select(-name_id),
      na.rm = TRUE
    ),
    
    mean_risDel = rowMeans(
      risDel_df %>% select(-name_id),
      na.rm = TRUE
    ),
    
    mean_allSal = rowMeans(
      bind_cols(
        salWt_df %>% select(-name_id),
        salDel_df %>% select(-name_id)
      ),
      na.rm = TRUE
    ),
    
    mean_allRis = rowMeans(
      bind_cols(
        risWt_df %>% select(-name_id),
        risDel_df %>% select(-name_id)
      ),
      na.rm = TRUE
    ),
    
    mean_allWt = rowMeans(
      bind_cols(
        salWt_df %>% select(-name_id),
        risWt_df %>% select(-name_id)
      ),
      na.rm = TRUE
    ),
    
    mean_allDel = rowMeans(
      bind_cols(
        salDel_df %>% select(-name_id),
        risDel_df %>% select(-name_id)
      ),
      na.rm = TRUE
    )
  )
  
  
  expression_df <- bind_cols(
    salWt_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("salWt_", .x)),
    
    salDel_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("salDel_", .x)),
    
    risWt_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("risWt_", .x)),
    
    risDel_df %>% 
      select(-name_id) %>% 
      rename_with(~ paste0("risDel_", .x))
  )
  
  
  cluster_mean_expression_list[[cluster_name]] <- bind_cols(
    mean_df,
    expression_df
  )
  
  
  cluster_mean_long_df[[cluster_name]] <- mean_df %>%
    mutate(cluster_id = cluster_name) %>%
    select(
      name_id,
      cluster_id,
      mean_salWt,
      mean_salDel,
      mean_risWt,
      mean_risDel,
      mean_allSal,
      mean_allRis,
      mean_allWt,
      mean_allDel
    )
}


cluster_mean_long_df <- bind_rows(cluster_mean_long_df)


edger_df <- read_excel(
  "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_09.09.2025.xlsx"
)


edger_df_with_means <- edger_df %>%
  left_join(
    cluster_mean_long_df %>%
      select(
        name_id,
        cluster_id,
        mean_salWt,
        mean_salDel,
        mean_risWt,
        mean_risDel,
        mean_allSal,
        mean_allRis,
        mean_allWt,
        mean_allDel
      ),
    by = c(
      "name_id" = "name_id",
      "cluster" = "cluster_id"
    )
  )


write_xlsx(
  edger_df_with_means,
  "results/risperidone-3q29/edger_statistics/edgerMarginalProportional_combineMultiComparision_pInteraction0.05_withMeans_09.06.2026.xlsx"
)


meanExpressionPerSampleMin20spots_quantileNormalizeData_df <- bind_rows(
  cluster_mean_expression_list,
  .id = "cluster"
)


write_xlsx(
  meanExpressionPerSampleMin20spots_quantileNormalizeData_df,
  "results/risperidone-3q29/edger_statistics/meanExpressionPerSampleMin20spots_quantileNormalizeData_09.06.2026.xlsx"
)


rm(
  cluster_names,
  cluster_name,
  salWt_df,
  salDel_df,
  risWt_df,
  risDel_df,
  mean_df,
  expression_df,
  cluster_mean_expression_list,
  cluster_mean_long_df,
  edger_df,
  edger_df_with_means,
  meanExpressionPerSampleMin20spots_quantileNormalizeData_df
)
