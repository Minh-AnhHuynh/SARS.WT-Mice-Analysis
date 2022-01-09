
################################################################################
# Comparison of all strains
################################################################################

### MDS with complete dataset, separated per strain => we observe outliers and discriminate WT vs rest

MDS <- make_MDS(dataset)
plot_MDS(MDS, "strain", legend = "Strain", "complete_strain", outlier = TRUE)

# Create a dataset without outliers
dataset_no_outliers <- my_subset(dataset, remove = c("mock_d4_4", "SARS.WT.104_d7_4", "SARS.WT.105_d2_1") )

### DEG virus strain vs mock => there are DEG only for WT vs mock

# SARS.BatSRBD vs mock : 0 gene
cond1 <- "SARS.BatSRBD"
cond2 <- "mock"
DEG_Bat <- deg(dataset_no_outliers, cond1, cond2)
list_genes <- DEG_Bat[DEG_Bat$significant == TRUE,"gene"]
print(paste("Differentially expressed genes for", cond1, "vs", cond2, ":", length(list_genes) ))

# SARS.icSARS vs mock : CXCL10
cond1 <- "SARS.icSARS"
cond2 <- "mock"
DEG_ic <- deg(dataset_no_outliers, cond1, cond2)
list_genes <- DEG_ic[DEG_ic$significant == TRUE,"gene"]
print(paste("Differentially expressed genes for", cond1, "vs", cond2, ":", length(list_genes) ))

# SARS.WT vs mock : 317 genes
cond1 <- "SARS.WT"
cond2 <- "mock"
DEG_WT <- deg(dataset_no_outliers, cond1, cond2)
DEG_WTvsmock <- DEG_WT[DEG_WT$significant == TRUE,"gene"]
print(paste("Differentially expressed genes for", cond1, "vs", cond2, ":", length(DEG_WTvsmock) ))
# Print the genes with a volcano plot
volcano_plot(DEG_WT, cond1, cond2, ngene = 15)

# EnrichR of DEG WT vs mock => activation of chemotaxis, mostly NK cells
enriched <- list_enriched(DEG_WT[DEG_WT$significant == TRUE,"gene"], DEG_WT[DEG_WT$significant == TRUE,"dir"])
plot_func(enriched, "WTvsmock", nterm = 15)

### Heatmap with DEG WT vs mock => we get two clusters of genes and two clusters of strains

# Heatmap using only the genes of WTvsmock
hm_res <- my_heatmap(data = dataset, name = "with_DEG_WTvsmock", DEG_WTvsmock, 
                     show_keep = "_", show_remove = "&", divide_cond = "mock", dendrogram = 'both',
                     dendrorow = 2, dendrocol = 2)

