if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


install.packages("preprocessCore")
BiocManager::install("edgeR")
require(DESeq2)


install.packages("enrichR")
require(enrichR)

install.packages("gplots")
# If necessary, install the e1071 package
if (!require(e1071)) {
  install.packages("e1071")
}

# Load the e1071 package
library(e1071)

install.packages('psych')
install.packages("DT")
require(DT)

library(shinydashboard)
install.packages("shinydashboard")