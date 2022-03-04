
################################################################################
# Example code (no function) for gene sets manipulation:
#   - generate gene sets with deg()
#   - put the gene sets together with Venn()
#   - plot the Venn object with setmap() or ggplot()
#   - manipulate gene sets with unite(), discern(), overlap()
################################################################################

.example <- function(){
  
  ### Preparation
  
  
  # D1 : WT vs mock
  cond1 <- c("SARS.WT.104_d1", "SARS.WT.105_d1")
  cond2 <- "mock_d1"
  DEG_d1 <- deg(dataset_no_outliers, cond1, cond2)
  genes1 <- DEG_d1[DEG_d1$significant == TRUE, "gene"]
  # D2 : WT vs mock
  cond1 <- c("SARS.WT.105_d2", "SARS.WT.105_d2")
  cond2 <- "mock_d2"
  DEG_d2 <- deg(dataset_no_outliers, cond1, cond2)
  genes2 <- DEG_d2[DEG_d2$significant == TRUE, "gene"]
  # D4 : WT vs mock
  cond1 <- c("SARS.WT.104_d4", "SARS.WT.105_d4")
  cond2 <- "mock_d4"
  DEG_d4 <- deg(dataset_no_outliers, cond1, cond2)
  genes4 <- DEG_d4[DEG_d4$significant == TRUE, "gene"]
  # D7 : WT vs mock
  cond1 <- c("SARS.WT.104_d7", "SARS.WT.105_d7")
  cond2 <- "mock_d7"
  DEG_d7 <- deg(dataset_no_outliers, cond1, cond2)
  genes7 <- DEG_d7[DEG_d7$significant == TRUE, "gene"]
  
  
  ### 4 gene sets (1 per day) to 3 gene sets (early, late, both)
  
  
  # We have 4 gene sets, 1 per day
  print(paste("Day 1 :", length(genes1), "genes"))
  print(paste("Day 2 :", length(genes2), "genes"))
  print(paste("Day 4 :", length(genes4), "genes"))
  print(paste("Day 7 :", length(genes7), "genes"))
  
  # Comparison of these 4 genes sets with a setmap
  quad_set <- Venn(list("day 1" = genes1, "day 2" = genes2, "day 4" = genes4, "day 7" = genes7))
  setmap(quad_set, element_fontsize = 0.1)
  dev.off()
  
  pdf(file = "./04_figures/Setmap_example.pdf")
  setmap(quad_set, element_fontsize = 0.1)
  dev.off()
  
  # Merged : "early" days = d1/d2, "late" days = d4/d7
  early <- unite(quad_set, slice = c("day 1", "day 2"))
  late <- unite(quad_set, slice = c("day 4", "day 7"))
  
  # Comparison of the two gene sets with a Venn diagramm
  duo_set <- Venn(list( "early" = early, "late" = late ))
  ggvenn(duo_set) + theme(panel.background = element_rect(fill = "white"))
  pdf(file = "./04_figures/Venn_example.pdf")
  ggvenn(duo_set) + theme(panel.background = element_rect(fill = "white"))
  dev.off()
  
  # Sets that are unique to early/late days or that are present in both
  inter_set <- overlap(duo_set)
  early_set <- discern(duo_set, "early")
  late_set <- discern(duo_set, "late")
  
  # We end up with 3 gene sets
  print(paste("Early :", length(early_set), "genes"))
  print(paste("Late :", length(late_set), "genes"))
  print(paste("Both :", length(inter_set), "genes"))
  
}

