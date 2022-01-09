
################################################################################
# Code for:
#   - function to find differentially expressed genes: deg(data, cond1, cond2)
#   - function to plot deg on a volcano plot: volcano_plot(DEG, cond1, cond2, ngene = 15)
################################################################################

deg <- function(data, cond1, cond2){
  # Returns a data frame with the columns:
  # $ gene        : names of all genes tested
  # $ padjust     : p-value for the gene when comparing conditions "cond1" and "cond2"
  # $ fc          : fold change for the gene (cond1/cond2 so >1 => cond1>cond2)
  # $ significant : logical value, equals TRUE when padjust < alpha (0.01) and |fc| > beta (1.2)  
  
  # Retrieve the 2 datasets corresponding to cond1 and cond2
  data1 <- my_subset(data, keep = cond1)
  data2 <- my_subset(data, keep = cond2)
  
  # Calculate a t-test cond1 vs cond2 for each gene
  test <- row_t_equalvar(data1, data2)
  
  # Add a column with the names of the gene, with adjusted p-values, fold changes, and change direction
  test$gene <- rownames(test)
  test$padjust <- p.adjust(test$pvalue, method = 'fdr')
  test$fc <- rowMeans2(as.matrix(data1) / rowMeans2(as.matrix(data2)))
  test$dir <- sign(test$fc-1)

  test <- test[,c("gene", "padjust", "fc", "dir")]
  
  # the row_t_equalvar test might create some NA's if gene have similar expression in both samples
  test[is.na(test$padjust),"padjust"] <- 1
  
  # Add a column indicating if the gene is differentially expressed
  test$significant <- test$padjust < alpha & (test$fc > beta | test$fc < 1/beta)
  
  # save differentially expressed genes to .txt
  name1 <- paste(cond1, collapse=".")
  name2 <- paste(cond2, collapse=".")
  filename <- paste0("./03_results/DEG_", name1, "vs", name2, ".txt")
  write.table(test[test$significant==TRUE,"gene"], filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(test)
}

volcano_plot <- function(DEG, cond1, cond2, ngene = 15){
  # Plots the result of "deg" as a volcano plot, and saves it as "Volcano_*cond1*vs*cond2*.pdf"
  # ngene indicates the number of genes to label on the plot
  
  # Create a column "col" indicating how to color each gene
  DEG$col <- "noDE"
  DEG$col[DEG$significant==TRUE & DEG$fc<(1/beta)] <- "down"
  DEG$col[DEG$significant==TRUE & DEG$fc>beta] <- "up"
  
  # Add a label for the first "ngene" genes with the best fc and padjust, NA for others labels
  DEG$sorted <- 1
  DEG$sorted[DEG$fc>beta | DEG$fc<(1/beta)] <- DEG$padjust[DEG$fc>beta | DEG$fc<(1/beta)]
  DEG$lab <- NA
  DEG$lab[order(DEG$sorted, decreasing = FALSE)][1:ngene] <- DEG$gene[order(DEG$sorted, decreasing = FALSE)][1:ngene]
  
  # Generate the volcano plot
  p <- ggplot(DEG, aes(x = log2(fc), y = -log10(padjust), fill = col, label = lab)) +
    geom_point(color = "black", shape = 21) +
    scale_fill_manual(values = c("down" = "green", "noDE" = "gray90", "up" = "red" ) ) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", size = 0.75) + 
    geom_vline(xintercept = c(log2(beta) , log2(1/beta) ), linetype = "dashed", size = 0.75) +
    labs(x = "log2(fold change)", y="-log10(p_value)", fill = "Legend:") +
    geom_text_repel(max.overlaps = ngene/2) +
    theme_minimal() +
    theme(text = element_text(size = 14)) +
    xlim(-1, 1)
  print(p)
  
  # Save the volcano plot to pdf
  name1 <- paste(cond1, collapse=".")
  name2 <- paste(cond2, collapse=".")
  filename <- paste0("./04_figures/Volcano_", name1, "vs", name2, ".pdf")
  pdf(filename)
  print(p)
  dev.off()
  
}

################################################################################
# Examples
################################################################################

.example <- function(){
  
  # To perform a differential gene expression analysis, we must first generate the list of
  # differentially expressed genes with "deg" and then plot it with "volcano_plot()"
  
  # DEG of WT vs mock
  cond1 <- "WT"
  cond2 <- "mock"
  DEG <- deg(dataset_no_outliers, cond1, cond2)
  
  # The list of deg genes can be accessed through its "significant" column
  DEG_list <- DEG[DEG$significatif == TRUE,]
  DEG_names <- DEG_list$gene
  
  # volcano_plot() is applied on the DEG result and saved as "Volcano_*example_cond1*vs*cond2*.pdf"
  volcano_plot(DEG, paste0("example_", cond1), cond2)
  
  # the number of labeled plotted is controlled with "ngene"
  volcano_plot(DEG, paste0("example_ngene100_", cond1), cond2, ngene = 100)
  
}
