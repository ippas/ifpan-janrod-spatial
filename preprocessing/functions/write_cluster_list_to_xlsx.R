write_cluster_list_to_xlsx <- function(cluster_list, file_path) {
  wb <- openxlsx::createWorkbook()
  
  sorted_names <- names(cluster_list)[order(as.integer(gsub(".*_(\\d+)$", "\\1", names(cluster_list))))]
  
  for (sheet_name in sorted_names) {
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = cluster_list[[sheet_name]])
    openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  }
  
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
}
