
################################################################################
# Code for:
#   - function to make an MDS: make_MDS(data)
#   - function to plot an MDS: plot_MDS(MDS, param, name, outlier = TRUE)
################################################################################

make_MDS <- function(data){
  # Returns the result of the MDS algorithm applied on "data"
  # Adds columns used for plotting different categories of objects with the same color
  
  # Make data frame of points
  d_matrix <- dist(t(data))
  fit <- isoMDS(d_matrix,k = 2)
  MDS <- data.frame(fit$points)
  colnames(MDS) <- c("x","y")
  
  # Add columns useful for the plot
  MDS$sample <- colnames(data)
  MDS$strain_inoc <- as.vector(sapply(MDS$sample,function(x) strsplit(x,"\\_")[[1]][1]))
  MDS$strain <- as.vector(sapply(MDS$strain_inoc,function(x) gsub('.104|.105',replacement = '',x) ))
  MDS$day <- as.vector(sapply(MDS$sample,function(x) strsplit(x,"\\_")[[1]][2]))
  MDS$strain_inocXday <- str_c(MDS$strain_inoc, '_', MDS$day)
  MDS$strainXday <- str_c(MDS$strain, '_', MDS$day)
  
  return(MDS)
}

plot_MDS <- function(MDS, param, name, legend = '', outlier = TRUE){
  # Makes a plot of an "MDS" result and save it as *MDS_'name'.pdf*
  # "param" is the name of column of the MDS data frame, by which we can color the points
  # "legend" is the title to put on the legend
  # if "outlier" is TRUE, then labels pointing toward the outliers are added to the plot
  
  # Make the plot
  plot <- ggplot(data = MDS) +
    geom_point(aes(x, y, fill = MDS[,param]), shape = 21, color = "black", size = 2) +
    labs(fill = legend) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_blank(),
          legend.position="right",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
    ) +
    if (outlier == TRUE) {
      geom_label_repel(aes(x, y, label = ifelse(sample %in% c("mock_d4_4", "SARS.WT.104_d7_4", "SARS.WT.105_d2_1"), sample, "")), 
                       max.overlaps = Inf, box.padding = 1, min.segment.length = 0)
    }
  
  # Save as pdf
  pdf(file = paste0("./04_figures/MDS_",name,".pdf"))
  plot(plot)
  dev.off()
  
  return(plot)
}

################################################################################
# Examples
################################################################################

.example <- function(){
  
  # To perform an MDS on a dataset, we first must apply the function "make_MDS()" to the dataset
  # Then we must apply the function "plot_MDS" to the make_MDS result
  # The second argument "param" of "plot_MDS", which selects categories to color can be chosen 
  # from : "sample", "strain_inoc", "strain", "day", "strain_inocXday", "strainXday"

  # Perform an MDS with the complete dataset
  MDS1 <- make_MDS(dataset)
  # Color by strain, give a name for the legend, and name the file "example_complete_strain"
  plot_MDS(MDS1, "strain", legend = "Strain:", "example_complete_strain")
  
  # Color by day, and hide outliers
  plot_MDS(MDS1, "day", "example_complete_day_nooutliers", outlier = FALSE)
  
  # Perform an MDS on a sub-dataset
  data_subset <- my_subset(dataset, keep = "WT")
  MDS2 <- make_MDS(data_subset)
  plot_MDS(MDS2, "day", "example_WT_day")
  
}

