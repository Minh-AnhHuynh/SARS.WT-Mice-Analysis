
################################################################################
# Code for:
#   - function to find enriched terms in a gene list: list_enriched(genes)
#   - function to plot the list of enriched terms: plot_func(enriched, name = '', nterms = 15)
################################################################################

list_enriched <- function(genes, reg){
  # Returns the list of significantly enriched terms in the given vector "genes"
  # Currently only looks for Gene Ontology 2021 databases
  
  setEnrichrSite("Enrichr")
  dbs_used <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021")
  
  res <- data.frame()
  for (d in c(-1, 1)){
    dir_genes <- genes[reg == d]
    if (length(dir_genes > 0)){
      enriched <- enrichr(dir_genes, dbs_used)
      enriched_sig <- lapply(enriched, function(x) x[x$Adjusted.P.value<alpha,])
      
      for (dbs in 1:length(dbs_used)){
        if (nrow(enriched_sig[[dbs]]) > 0){
          name <- dbs_used[[dbs]]
          enriched_sig[[dbs]]$database <- name
          enriched_sig[[dbs]]$dir <- ifelse(d==1, "red", "green")
          res <- rbind(res, enriched_sig[[dbs]])
        }
      }
    }
    
  }
  
  return(res)
}

plot_func <- function(enriched, name = '', nterms = 15, color = TRUE){
  # Plots a list of enriched terms obtained with list_enriched() and saves it as Enrich_*name*.pdf
  # Only plots the "nterms" most significant ones

  # If there is a non-empty list of enriched terms
  if (nrow(enriched) > 0) {
    
    # Obtain the nterms with the best Combined.Score
    if (nrow(enriched) > nterms) {
      enriched <- enriched[order(enriched$Combined.Score, decreasing = TRUE)[1:nterms],]
    }
    
    # Color the names of the terms with fold change direction
    colors <- enriched$dir[order(enriched$Combined.Score, decreasing = FALSE)]
    if (color == FALSE) { colors = "black" }
    
    # Generate the plot
    p <- ggplot(data=enriched, 
                aes(x = Combined.Score, y = reorder(Term, Combined.Score ), fill = database)
    ) +
      scale_fill_manual(values = c("GO_Biological_Process_2021" = "#F8766D",
                                   "GO_Cellular_Component_2021" = "#7CAE00",
                                   "GO_Molecular_Function_2021" = "#00BFC4")) +
      geom_col() +
      labs(y = "", x = "Combined score", fill = "Database") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid = element_line(color = "grey80"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = colors ),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.background = element_rect(color = "black"),
        legend.position = "right",
        plot.tag.position = c(0.75, 0.3),
      )
  
    # Save it as pdf
    filename <- paste0("./04_figures/Enrich_", name, ".pdf")
    print(p)
    pdf(filename)
    print(p)
    dev.off()
    
  } else {print("Nothing is enriched.")}
  
}

################################################################################
# Examples
################################################################################

.example <- function(){
  
  # The function list_enrich() takes 1 arguments :
  # "genes" : the list of genes into which we look for over-represented terms
  # The function plot_funct() takes 3 arguments :
  # "enriched" : a result from list_enrich() to be plotted
  # "name" : name for the created pdf file
  # "nterms" : max number of terms to be plotted
  
  # Create a gene set to analyze using deg()
  DEG <- deg(dataset, "WT", "mock")
  genes <- DEG[DEG$significant == TRUE,"gene"]
  reg <- DEG[DEG$significant == TRUE,"dir"]
  
  # Generate the list of enriched terms
  func_res <- list_enriched(genes, reg)
  
  # Plot the 10 most enriched terms, and save the file to pdf as "Enrich_example"
  plot_func(func_res, name = "example", nterms = 10)
  
}


