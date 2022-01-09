
################################################################################
# Comparison of the response to SARS-WT strain at various time points
################################################################################

### MDS of SARS.WT strain => SARS.WT strain clusters well according to days but not to inoculum

# SARS.WT without outliers
data_SARS.WT <- my_subset(dataset_no_outliers, keep = "SARS.WT")
MDS <- make_MDS(data_SARS.WT)

# Separate with inoculum
plot_MDS(MDS, "strain_inoc", legend = "Inoculum", "SARS.WT_inoc")
# deg inoculum : 104 vs 105 => 0 gene
cond1 <- "SARS.WT.104"
cond2 <- "SARS.WT.105"
DEG_inoc <- deg(dataset_no_outliers, cond1, cond2)
genes_inoc <- DEG_inoc[DEG_inoc$significant == TRUE, "gene"]

# Separate with day the SARS.WT strain (all genes)
plot_MDS(MDS, "day", legend = "Day", "SARS.WT_day")
# MDS on mock using only the genes of WTvsmock => we should separate mock d1 to d7
data_mock <- my_subset(dataset_no_outliers, keep = "mock")
data_MDS <- data_mock[rownames(data_mock) %in% DEG_WTvsmock,]
MDS <- make_MDS(data_MDS)
# Separate with day the mock strain (only deg WT vs mock)
plot_MDS(MDS, "day", legend = "Day", "mock_day")

### Lists of deg for each day, and deg 104vs105

# D1 : WT vs mock
cond1 <- c("SARS.WT.104_d1", "SARS.WT.105_d1")
cond2 <- "mock_d1"
DEG_d1 <- deg(dataset_no_outliers, cond1, cond2)
genes1 <- DEG_d1[DEG_d1$significant == TRUE, "gene"]
# D2 : WT vs mock
cond1 <- c("SARS.WT.104_d2", "SARS.WT.105_d2")
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

### Comparison of these 4 genes sets with a setmap => early days and late days

quad_set <- Venn(list("day 1" = genes1, "day 2" = genes2, "day 4" = genes4, "day 7" = genes7))
setmap(quad_set, element_fontsize = 0.1)
dev.off()

pdf(file = "./04_figures/Setmap_4_days.pdf")
setmap(quad_set, element_fontsize = 0.1)
dev.off()

### Merged : "early" days = d1/d2, "late" days = d4/d7

# early : 404 genes
cond1 <- c("SARS.WT.104_d1", "SARS.WT.105_d1", "SARS.WT.104_d2", "SARS.WT.105_d2")
cond2 <- c("mock_d1", "mock_d2")
early <- deg(dataset_no_outliers, cond1, cond2)
early_genes <- early[early$significant == TRUE, "gene"]
early_reg <- early[early$significant == TRUE, "dir"]
# late : 446 genes
cond1 <- c("SARS.WT.104_d4", "SARS.WT.105_d4", "SARS.WT.104_d7", "SARS.WT.105_d7")
cond2 <- c("mock_d4", "mock_d7")
late <- deg(dataset_no_outliers, cond1, cond2)
late_genes <- late[late$significant == TRUE, "gene"]
late_reg <- late[late$significant == TRUE, "dir"]

### Comparison of early and late gene sets

# Comparison of the two gene sets with a Venn diagramm
duo_set <- Venn(list( "early" = early_genes, "late" = late_genes ))
ggvenn(duo_set) + theme(panel.background = element_rect(fill = "white"))
pdf(file = "./04_figures/Venn_earlyvslate.pdf")
ggvenn(duo_set) + theme(panel.background = element_rect(fill = "white"))
dev.off()

# Genes that are common or unique to early and late sets
inter_set <- overlap(duo_set)
u_early_genes <- early_genes[!early_genes %in% late_genes]
u_early_reg <- early_reg[!early_genes %in% late_genes]
u_late_genes <- late_genes[!late_genes %in% early_genes]
u_late_reg <- late_reg[!late_genes %in% early_genes]

### Enrichr of what these genes lists contain

# early genes => broad innate immune response
enriched <- list_enriched(u_early_genes, u_early_reg)
plot_func(enriched, "unique_early_set", nterm = 10)

# late genes => mitosis blockage
enriched <- list_enriched(u_late_genes, u_late_reg)
plot_func(enriched, "unique_late_set", nterm = 10)

### Fold change heatmap on inter_set (overlapping genes) to know when these 
### genes are up or down regulated

# days cluster well according to these genes
heatmap <- my_heatmap(data = dataset_no_outliers,
                      name = "inter_set",
                      genes = inter_set,
                      show_keep = "SARS.WT",
                      show_remove = "&",
                      divide_cond = "mock",
                      dendrogram = 'both',
                      dendrorow = 4,
                      dendrocol = 4)
# There are 4 major types of gene expression
lists <- get_lists(heatmap, 4)

# 85 genes are upregulated in d1/d2, downregulated in d4/d7
list_1 <- lists[[1]]
# they correspond to positive regulation of NK cell chemotaxis
enriched <- list_enriched(list_1, rep(1, length(list_1)))
plot_func(enriched, "dore", 10, color = FALSE)

# 52 genes are upregulated in d2 only
list_2 <- lists[[2]]
# they correspond to CXCR3 binding / T cell migration
enriched <- list_enriched(list_2, rep(1, length(list_2)))
plot_func(enriched, "rose", 10, color = FALSE)

# 13 genes are upregulated in d4/d7, downregulated in d1/d2
list_3 <- lists[[3]]
# they correspond to "Nothing is enriched."
enriched <- list_enriched(list_3, rep(1, length(list_3)))
plot_func(enriched, "violet", 10, color = FALSE)

# 24 genes are upregulated in d1 only
list_4 <- lists[[4]]
# they correspond to lymphocyte chemotaxis
enriched <- list_enriched(list_4, rep(1, length(list_4)))
plot_func(enriched, "turquoise", 10, color = FALSE)


