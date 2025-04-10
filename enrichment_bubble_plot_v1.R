setwd(<working_dir>)
list.files()

# This script is based on GOmeth output from MissMethyl package
kegg_res <- read.csv("KEGG_missMethyl.csv", header = T)
head(kegg_res)

# Calculate fold enrichment 
# Enrichment = (DEG in pathway/total DEG)/(Genes in Pathway/Total Genes)
head(kegg_res)[,1:5]

# Fold enrichment
# Let's set total DEG to 1000
total_deg <- 1000 # Substitute with actual number of DEG in the comparison
total_genes <- 20000 # Substitute with actual number of genes compared 
fe <- (kegg_res$DE/total_deg)/(kegg_res$N/total_genes)
head(fe)

# Add enrichment fold column
kegg_res$FE <- fe
head(kegg_res)

# Keep significant KEGG pathways
sig <- kegg_res[kegg_res$FDR < 0.05,]
dim(sig)

# Bubble plot
ggplot(sig, aes(x = FE, y = reorder(Description, -FE))) +
  geom_point(aes(size = DE, color = FDR)) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR") +
  scale_size(range = c(0.5, 5), name = "Number of DM genes") +
  theme_minimal() +
  labs(
    x = "Fold Enrichment",
    y = "Pathway ID",
    title = "Bubble Plot of Pathway Enrichment"
  )
ggsave("enrichment_bubble_plot.pdf", device = "pdf", width = 8, height = 9)
