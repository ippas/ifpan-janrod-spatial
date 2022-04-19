#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Feb 2022
#############################################################

######################################
# Install and require needed package #
######################################
install.packages("ggplot2")
install.packages("Matrix")
install.packages("rjson")
install.packages("cowplot")
install.packages("RColorBrewer")
install.packages("Seurat")
install.packages("grid")
install.packages("readbitmap")
install.packages("dplyr")
install.packages("data.table")
install.packages("doSNOW")
install.packages("hdf5r")
install.packages('remotes')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("rhdf5")

require(remotes)
remotes::install_github("satijalab/seurat-data")
require(SeuratData)
require(patchwork)
require(ggplot2)
require(pbapply)
require(Seurat)
require(stringr)
require(ggplot2)
require(Matrix)
require(rjson)
require(cowplot)
require(RColorBrewer)
require(grid)
require(readbitmap)
require(Seurat)
require(dplyr)
require(hdf5r)
require(data.table)
require(doSNOW)
require(tidyr)
require(tibble)
require(gridExtra)




# Require source files
# get list of sample
samples_name <- list.files(path = "data/ldopa/spaceranger-corrected/")

samples_name


for(sample_name in samples_name) {
  directory <- paste("data/ldopa/spaceranger-corrected/", 
                     sample_name, 
                     "/outs/filtered_feature_bc_matrix/", sep = "")
  print(directory)
  sample_data <- Read10X(
    data.dir = directory
  )
  
  sample <- CreateSeuratObject(counts = sample_data,
                               project = sample_name)
  sample <- AddMetaData(sample, metadata=sample_name, col.name = "sample")
  
  assign(sample_name, sample)
  
  rm(sample_name,
     directory,
     sample_data,
     sample)
}

###################
# Merging samples #
###################
merged_samples <- merge(
  get(samples_name[1]),
  sapply(samples_name[2:length(samples_name)], get),
  add.cell.ids = samples_name,
  project = "merged.sample"
)



####################
# integration data #f
####################
# split the dataset into a list
merged_samples_list <- SplitObject(merged_samples, split.by = "sample")

# normalize and identify variable features for each dataset independently
merged_samples_list <- lapply(X = merged_samples_list,
                              FUN = function(x) {
                                x <- NormalizeData(x, normalization.method = "LogNormalize")
                                x <- FindVariableFeatures(x, 
                                                          selection.method = "vst",
                                                          nfeatures = 2000)
                              })


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merged_samples_list, nfeatures = 2000)


merged_anchors <- FindIntegrationAnchors(object.list = merged_samples_list,
                                         anchor.features = features)

# this command creates an 'integrated' data assay
integrated_data <- IntegrateData(anchorset = merged_anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated_data) <- "integrated"


integrated_analysis <- ScaleData(integrated_data, verbose = FALSE)
integrated_analysis <- RunPCA(integrated_analysis, npcs = 30, verbose = FALSE)
integrated_analysis <- RunUMAP(integrated_analysis, reduction = "pca", dims = 1:30)
integrated_analysis <- FindNeighbors(integrated_analysis, reduction = "pca", dims = 1:30)




# remove varible 
rm(list = samples_name)
rm(integrated_data,
   merged_anchors,
   merged_samples,
   merged_samples_list)
