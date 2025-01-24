library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(enrichplot)

# First run enrichment analysis with clusterProfiler
# Use the resulting object to calculate pairwise distances
# Next, cluster and create treeplot
sim_ora <- pairwise_termsim(ego) # ego - clusterProfiler output
print(sim_ora)

# Checked, how many terms have adjusted p-value below 0.05
sum(ego@result$p.adjust < 0.05)

# Cluster and create treeplot in one step
treeplot(sim_ora, showCategory = 50,
         color = "p.adjust", nWords = 4, nCluster = 5,
         hclust_method = "ward.D")
ggsave("treeplot.pdf", device = "pdf", units = 'in',
       width = 12, height = 10)
