library(DESeq2)
library(ComplexHeatmap)
library(ggplot2)

setwd(<Working_dir>)

load("expression_data/VST_data_object_replicates.RData")
vsd

# Get samples of interest
colData(vsd)
vsd_select <- 
  vsd[,colData(vsd)$Group == "BT12_20" | colData(vsd)$Group == "BT12_DMSO"]

# Get expression data
exp <- assay(vsd_select)
exp <- as.data.frame(exp)
head(exp)

# Get results data frame
setwd("comparisons/BT12/BT12_20_vs_BT12_DMSO/")
res <- read.table("BT12_20_vs_BT12_DMSO_noiseq_results.txt", sep = "\t",
                  header = T, quote = "\"", fill = T)
head(res)

# How many differentially expressed genes are present?
# We will use prob > 0.9 as a threshold following developer's
# recommendation
res_sig <- subset(res, prob > 0.9)
dim(res_sig)
res_annot <- merge(res_sig, exp, by.x = "Row.names", by.y = 0)

# Get top 50 genes by ranking
top50 <- res_annot[order(abs(res_annot$ranking), decreasing = T),]
top50 <- top50[1:50,]
# Sort by fold change
top50 <- top50[order(top50$M, decreasing = T),]
mat <- top50[,grep("Megan", names(top50))]
row.names(mat) <- top50$hgnc_symbol                 
head(mat)
mat <- as.matrix(mat)
heat <- t(scale(t(mat)))
head(heat)

# Set colour scheme and choose breaks
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

# Create bar chart to place at the right side of the heatmap
# with log2 fold changes
fc_df <- data.frame(symbol = top50$hgnc_symbol, 
                    FC = top50$M)
rownames(fc_df) <- fc_df$symbol
fc_df <- fc_df[,-1]
row_ann <- rowAnnotation(bar = anno_barplot(fc_df), 
                         show_annotation_name = F)
colors <- list(Megan2="lightgreen", Megan4="dodgerblue",
               Megan2_1="lightgreen", Megan4_1="dodgerblue")

treatments = c("DMSO", "Ext20", "DMSO", "Ext20")

set.seed(1)
col_ann <- HeatmapAnnotation(Treatment = treatments, 
                             show_legend = F,
                             col = colors)

# Create heatmap
pdf("BT12_20_vs_DMSO_complex_heatmap.pdf", width = 7, height = 8)
Heatmap(heat, cluster_rows = F, clustering_distance_columns = "euclidean",
        column_km = 2, rect_gp=gpar(col="white",lwd=1), 
        right_annotation = row_ann, show_column_dend = F,
        top_annotation = col_ann, column_title = c("Ext20", "DMSO"),
        heatmap_legend_param = list(title = NULL),
        show_column_names = F)
dev.off()
