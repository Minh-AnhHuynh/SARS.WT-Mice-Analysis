
################################################################################
# Code for:
#   - function to generate a heatmap: my_heatmap(data, genes, show_keep)
#   - function to obtain genes by cutting row_dendrogram: get_lists(heatmap, k) 
################################################################################

my_heatmap <- function(data, genes, show_keep, show_remove = "&", divide_cond = "mock", 
                       dendrogram = 'both', dendrorow = 2, dendrocol = 2, name = ''){
  # Returns a fold change heatmap from "data"
  # The genes to cluster with are "genes", the samples to cluster are chosen with 
  # "show_keep" and "show_remove" as for subset() function.
  # Fold changes are obtained by dividing the dataframe with "divide_cond"
  # A dendrogram can be performed on "row", "column", "both", or "none"
  # Meta clusters a formed and colored, the numbers of metaclusters are
  # "dendrorow" and "dendrocol"
  # The heatmap is saved as "Heatmap_*name*.pdf"
  
  # Select heatmap columns
  data_show <- my_subset(data, keep = show_keep, remove = show_remove)
  data_divide <- my_subset(data, keep = divide_cond)
  
  # Normalize by the mean of data_divide
  data_divide_mean <- rowMeans2(as.matrix(data_divide))
  data_hm <- data_show / data_divide_mean
  # /!\ should be replaced with the following line but we did not have enough 
  # time to describe the new resulting plots in the report
  # data_hm <- log(data_show / data_divide_mean)
  
  # Filter selected genes
  data_hm <- data_hm[rownames(data_hm) %in% genes,]
  
  # Prepare heatmap
  data_scaled <- t(scale(t(data_hm), center = TRUE, scale = TRUE))
  size_palette <- max( abs(min(data_scaled)), abs(max(data_scaled)) )
  ppalette <- seq(-size_palette,size_palette,0.05)
  hmcols <- colorRampPalette(c("green","black","red"))(length(ppalette)-1)
  
  # Color dendrograms
  Rowv <- data_scaled %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", k = dendrorow) %>% set("branches_lwd", 4) %>%
    rotate_DendSer(ser_weight = dist(data_scaled))
  Colv  <- data_scaled %>% t %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", k = dendrocol) %>% set("branches_lwd", 4) %>%
    rotate_DendSer(ser_weight = dist(t(data_scaled)))
  
  # Generate heatmap and generate it again to save as pdf
  heatmap <- heatmap.2(as.matrix(data_scaled), margins=c(10,10), scale="none", breaks=ppalette, col=hmcols, 
                       symbreaks=FALSE, trace="none", density.info="none", key=TRUE, dendrogram = dendrogram,
                       Rowv = Rowv, Colv = Colv)
  
  pdf(file = paste0("./04_figures/Heatmap_", name, ".pdf"), width = 10)
  heatmap <- heatmap.2(as.matrix(data_scaled), margins=c(10,10), scale="none", breaks=ppalette, col=hmcols, 
                       symbreaks=FALSE, trace="none", density.info="none", key=TRUE, dendrogram = dendrogram,
                       Rowv = Rowv, Colv = Colv)
  dev.off()
  
  return(heatmap)
}

get_lists <- function(heatmap, k, name = ''){
  # Cut the row dendrogram into k metaclusters
  # Return a list of the k vectors of gene names and save them as "List_*name*i.txt"
  
  # Get the row dendrogram and cut it
  dend <- heatmap$rowDendrogram
  cut_dend <- cutree(dend, k = k)
  
  # Fill a list with the k vectors of gene names
  lists_genes <- list()
  for (i in 1:k){
    # Obtain the 'i'th gene list 
    lists_genes[[i]] <- labels(dend)[cut_dend == i]
    # Save it as .txt
    filename <- paste0("./03_results/List_", name, i, ".txt")
    write.table(lists_genes[[i]], filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  return(lists_genes)
}

################################################################################
# Examples
################################################################################

.example <- function(){
  
  # Generate a heatmap, using a dataset, a list of DEG, and a set of conditions for samples
  # to show, a condition for samples to divide with, and a name
  data <- dataset_no_outliers
  DEG_dataframe <- deg(data, "WT", "mock")
  genes <- DEG_dataframe[DEG_dataframe$significant == TRUE,"gene"]
  show_keep <- c("SARS.WT.104", "SARS.WT.105")
  show_remove <- "d1"
  divide_cond <- "mock"
  heatmap_1 <- my_heatmap(data = data, genes = genes, show_keep = show_keep, show_remove = show_remove, 
                        divide_cond = divide_cond, name = "example_1")
  
  # Make a dendrogram only for rows, with K=3 meta clusters
  heatmap_2 <- my_heatmap(data = data, genes = genes, show_keep = show_keep, show_remove = show_remove, 
                        divide_cond = divide_cond, dendrogram = 'row', dendrorow = 3, name = "example_2")

  # Extract the K=3 gene lists
  get_lists(heatmap_2, 3, name = 'example_')
  
}