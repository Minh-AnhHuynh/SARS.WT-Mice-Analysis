
################################################################################
# Code for:
#   - packages required for the whole project
#   - subset function: my_subset(data, keep = "_", remove = "&")
#   - constants for the whole project (alpha, beta, dataset, dataset_no_outliers)
################################################################################

# Load packages
library(stringr) # string manipulation
library(ggplot2) # plots
library(ggrepel) # plots => labels

library(MASS) # MDS
library(matrixTests) # row_t_equalvar for deg
library(matrixStats) # rowMeans2 for deg
library(enrichR) # Enrichment analysis
library(gplots) # heatmap.2
library(dendextend) # Cut heatmap dendrogram
library(RVenn) # Venn diagramm + setmap

# Load the dataset
dataset <- read.delim("./02_donnees/Project-SARSMice.txt")

# Constants
alpha <- 0.01 # p-value threshold
beta <- 1.2   # fold-change threshold

my_subset <- function(data, keep = "_", remove = "&"){
  # Returns a subset of data columns
  # Containing any of the terms of "keep" but none of "remove" string vectors
  
  col_keep <- sapply(keep, function(x) grepl(x, colnames(data)))
  filter_keep <- apply(col_keep, 1, any)
  
  col_remove <- sapply(remove, function(x) grepl(x, colnames(data)))
  filter_remove <- apply(col_remove, 1, any)
  
  filter_ <- filter_keep & !filter_remove
  subset_ <- subset(data, select = filter_)
  
  return(subset_)
}

################################################################################
# Examples
################################################################################

.example <- function(){
  
  # The function my_subset() takes 3 arguments :
  # "data" : the dataset from which we want to obtain a subset
  # "keep" : string vector indicating column names that we want to keep from data
  # "remove" : string vector indicating column names that we want to remove from data
  
  # Create a dataset containing only mock
  dataset_mock <- my_subset(dataset, keep = "mock")
  
  # Create a dataset containing mock_d1 and WT.104_d1
  dataset_d1 <- my_subset(dataset, keep = c("mock_d1", "WT.104_d1"))
  
  # Create a dataset containing all WT except for d4
  dataset_no_d4 <- my_subset(dataset, keep = c("WT.104", "WT.105"), remove = "d4")
  
}

# Create a dataset without outliers
dataset_no_outliers <- my_subset(dataset, remove = c("mock_d4_4", "SARS.WT.104_d7_4", "SARS.WT.105_d2_1") )


