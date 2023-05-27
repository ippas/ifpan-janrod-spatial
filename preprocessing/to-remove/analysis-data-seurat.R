######################################
# Install and require needed packages #
######################################
# List of packages to install
packages <- c("ggplot2", "Matrix", "rjson", "cowplot", "RColorBrewer", "Seurat", 
              "grid", "readbitmap", "dplyr", "data.table", "doSNOW", "hdf5r",
              "remotes", "BiocManager", "tidyr", "tibble", "gridExtra", "parallel")

# Install packages if not already installed
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install specific versions and repositories
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.11")
BiocManager::install("rhdf5")

# Install from GitHub
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  require(remotes)
  remotes::install_github("satijalab/seurat-data")
}

# Load required packages
required_packages <- c("SeuratData", "patchwork", "pbapply", "stringr", "tidyr", "tibble", packages)
lapply(required_packages, require, character.only = TRUE)



####################
# Integration Data #
####################
path_to_data <- "data/risperidone/spaceranger-corrected//"


# Require source files
# get list of sample
samples_name <- list.files(path = path_to_data) #%>% .[-c(8,10,11,12)]

samples_name


for(sample_name in samples_name) {
  directory <- paste(path_to_data, 
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


# Split the dataset into a list
merged_samples_list <- SplitObject(merged_samples, split.by = "sample")

# Normalize and identify variable features for each dataset independently
merged_samples_list <- lapply(X = merged_samples_list, FUN = function(x) {
  x %>%
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merged_samples_list, nfeatures = 2000)

# Find integration anchors
merged_anchors <- FindIntegrationAnchors(object.list = merged_samples_list, anchor.features = features)

# Create an 'integrated' data assay
integrated_data <- IntegrateData(anchorset = merged_anchors)

# Set the default assay to the corrected (integrated) data
DefaultAssay(integrated_data) <- "integrated"

# Perform downstream analysis on the integrated data
integrated_analysis <- integrated_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)


path_to_data <- "data/risperidone/spaceranger-corrected//"
nfeatures <- 2000
dims <- 1:30

risperidone_integrate <- integrate_data(path_to_data, nfeatures, dims)

# Remove unnecessary variables
rm(list = samples_name)
rm(integrated_data, merged_anchors, merged_samples, merged_samples_list)
