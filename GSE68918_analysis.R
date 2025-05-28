library(data.table)

anno <- fread("GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt", skip = 8)
anno_preclean <- anno[anno$Obsolete_Probe_Id == ""]
anno_clean <- anno_preclean[anno_preclean$Source != "ILMN_Controls"] #Final filtered data 

anno_premap <- anno_clean[, .(Probe_Id, Symbol, Entrez_Gene_ID)]
anno_map <- anno_premap[anno_premap$Symbol != ""]


data2 <- fread("GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt", skip = 8)

# Step 1: Find which row contains "ID_REF"
lines <- readLines("GSE68918_series_matrix.txt")
header_line <- grep("^\"ID_REF\"", lines)

# Step 2: Use fread() to skip everything before that row
expr <- fread("GSE68918_series_matrix.txt", skip = header_line - 1)

setnames(expr, "ID_REF", "Probe_Id")
expr_annotated <- merge(anno_map, expr, by = "Probe_Id")
expr_by_gene <- expr_annotated[, lapply(.SD, mean, na.rm = TRUE), by = Symbol, .SDcols = patterns("^GSM")]

# Step 1: Set the gene symbols as row names
gene_symbols <- expr_by_gene$Symbol
rownames(expr_by_gene) <- expr_by_gene$Symbol

# Step 2: Remove the Symbol column
expr_by_gene$Symbol <- NULL

# Step 3: Convert to a numeric matrix
expr_matrix <- as.matrix(expr_by_gene)
mode(expr_matrix) <- "numeric"
rownames(expr_matrix) <- gene_symbols

metadata <- data.frame(
  Sample = colnames(expr_matrix),
  Condition = c("siSCR", "siKDM3A-B", "siSCR", "siKDM3A-B", "siSCR", "siKDM3A-B",
                "siSCR", "siKDM3A-B", "siSCR", "siKDM3A-B", "siSCR")
)

install.packages("ggfortify")  # Only once
install.packages("pheatmap")  
install.packages("BiocManager")
BiocManager::install("limma")

# Now to experiment with the data: 
# PCA plot
pca <- prcomp(t(expr_by_gene), scale. = TRUE)
library(ggplot2)
library(ggfortify)
pca_plot <- autoplot(pca, data = metadata, colour = 'Condition') + ggtitle("PCA of Gene Expression (siSCR vs siKDM3A-B)") + theme_minimal()
#Save the figure for my github
ggsave("figures/pca_plot.png", plot = pca_plot, width = 6, height = 5, dpi = 300)


# Heatmap of top variable genes
library(pheatmap)
top_var_genes <- head(order(apply(expr_by_gene, 1, var), decreasing = TRUE), 50)
pheat_map <- pheatmap(
  expr_by_gene[top_var_genes, ],
  fontsize_row = 6,          # smaller row font
  fontsize_col = 8,          # smaller column font
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,         # removes ugly grid lines
  cellwidth = 15,            # adjust if squished
  cellheight = 7            # space out rows a bit
)

#Finding statistically different genes 
library(limma)
metadata$Condition <- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "siSCR")

design <- model.matrix(~ Condition, data = metadata)

#Fit the linear model
rownames(expr_matrix)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)


results <- topTable(fit, coef = "ConditionsiKDM3A-B", number = Inf)

# logFC tells you how much the expression of a gene has changed between your treatment and control conditions
pre_sig_hits <- results[abs(results$logFC) > 1, , drop = FALSE]
pre_sig_hits <- pre_sig_hits[order(pre_sig_hits$adj.P.Val), , drop = FALSE]  # optional: sort by significance

sig_hits <- pre_sig_hits[pre_sig_hits$adj.P.Val < 0.05, , drop = FALSE]
rownames(sig_hits) <- rownames(pre_sig_hits)[pre_sig_hits$adj.P.Val < 0.05]


#Heatmap
# Normalize expression matrix by row (Z-score)
expr_norm <- t(scale(t(expr_matrix)))

# Gather the top 50 significant genes 
top_genes_for_plot <- rownames(sig_hits)[1:50]
genes_to_plot <- intersect(top_genes_for_plot, rownames(expr_norm))

# Create sample annotation
annotation_col <- data.frame(Condition = metadata$Condition)
rownames(annotation_col) <- metadata$Sample

# Heatmap
pheatmap(
  expr_norm[genes_to_plot, ],
  annotation_col = annotation_col,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  cellheight = 7,
  cellwidth = 15
)
# Save heatmap to a PNG file
png("figures/sig_genes_heatmap.png", width = 1200, height = 1000)
pheatmap(
  expr_norm[genes_to_plot, ],
  annotation_col = annotation_col,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  cellheight = 7,
  cellwidth = 15
)
dev.off()


