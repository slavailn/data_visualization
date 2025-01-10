library(clusterProfiler)
library(ggplot2)
# This sample scripts demonstrates over-representation analysis relative
# to KEGG pathways followed by visualization of the results as bubble plot
# using ggplot2 package

# target_ids - list of genes of interest
# universe_ids - list of genes in the universe, for example all 
# of the genes expressed in the experiment

# Run KEGG over-representation analysis with clusterProfiler
kk <- enrichKEGG(
  gene = target_ids,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_ids,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.1,
  use_internal_data = FALSE
)

head(kk@result)
write.csv(kk@result, file = "KEGG_ORA_results.csv")

# Plot enriched terms
head(kk@result)
# Extract results data frame
ora <- kk@result
# We need to calculate gene ratios
ratio_list <- strsplit(ora$GeneRatio, split = '/')
target_genes <- lapply(ratio_list, '[', 1)
universe_genes <- lapply(ratio_list, '[', 2)
geneRatio <- as.numeric(target_genes) / as.numeric(universe_genes)
head(geneRatio)

# Trim long pathway names, here I simply removed redundant expression 
# expression present in every pathway name
# Trim ' -  Mus musculus (house mouse)' expression
pathName <- paste(ora$ID, ora$Description, sep=":")
pathName <- gsub(" - Mus musculus \\(house mouse\\)", "", pathName)
data <- data.frame(pathName = pathName, geneRatio = geneRatio, 
                   geneCount = ora$Count, pValue = ora$pvalue, 
                   adjPval = ora$p.adjust)
head(data)

# Order by p-value
data <- data[order(data$pValue, decreasing = F),]

# Get top 20 categories
data <- head(data, n = 20)

# Order by count
data <- data[order(data$geneCount, decreasing = T),]
data$Significant <- data$adjPval < 0.2

ggplot(data, aes(x = geneRatio, y = reorder(pathName, geneCount),
                 size = geneCount, color=Significant)) + 
  geom_point(alpha = 0.7) +  
  labs(y="Pathway", x = "Gene ratio") +
  guides(color = guide_legend(title = "adj.pvalue < 0.2")) + 
  guides(size = guide_legend(title = "Gene count"))
ggsave("KEGG_ORA_bubble_plot.pdf", device = "pdf", units = 'in',
       width = 10, height = 8)
