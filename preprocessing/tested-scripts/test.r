install.packages("here")

library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)

# Load the data
peaks_macs3_default <- read.csv("data/risperidone-3q29/gene-annotation/macs3/ris3q29-26s_peaks.narrowPeak", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(peaks_macs3_default) <- c(
  "chrom",
  "start",
  "end",
  "name",
  "score",
  "strand",
  "signalValue",
  "pValue_log10",
  "qValue_log10",
  "peak_offset"
) 
  
peaks_macs3_default %>% 
  dplyr::mutate(len_peak = end - start)  %>% 
  .$len_peak %>% sd()

# Load the data
peaks_macs3_ext100shift50q001 <- read.csv("data/risperidone-3q29/gene-annotation/rmacs3-ris3q29-26s-ext100shift50q001/ris3q29-26s_ext100shift50q001_peaks.narrowPeak", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(peaks_macs3_ext100shift50q001) <- c(
  "chrom",
  "start",
  "end",
  "name",
  "score",
  "strand",
  "signalValue",
  "pValue_log10",
  "qValue_log10",
  "peak_offset"
) 

peaks_macs3_ext100shift50q001 %>% 
  dplyr::mutate(len_peak = end - start)  %>% 
  .$len_peak %>% sd()

peaks_macs3_ext100shift50q001 %>% 
  dplyr::mutate(len_peak = end - start)  %>% 
  dim

peaks_macs3_default %>% 
  dplyr::filter(chrom == "chr7") %>% 
  dplyr::filter(start > 105556963) %>% 
  dplyr::filter(start < 105559597)
