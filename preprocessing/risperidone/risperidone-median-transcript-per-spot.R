risperidone_st_data_half$raw_data$data  %>% colSums() %>% names %>% str_split("_", simplify = TRUE) %>%
  .[, 1] %>% unique()
