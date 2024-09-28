library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

# ComplexHeatmap tutorial by Kevin Blighe
# Link: https://github.com/kevinblighe/E-MTAB-6141

tmpfile <- tempfile()
download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/mat.tsv',
              tmpfile, method = 'auto')
mat <- read.table(tmpfile, sep = '\t', row.names = 1,
                  header = TRUE, stringsAsFactors = FALSE)

tmpfile <- tempfile()
download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/metadata.tsv',
              tmpfile, method = 'auto')
metadata <- read.table(tmpfile, sep = '\t', row.names = 1,
                       header = TRUE, stringsAsFactors = FALSE)

tmpfile <- tempfile()
download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/sig_genes.list',
              tmpfile, method = 'auto')
sig_genes <- read.table(tmpfile, sep = '\t',
                        header = FALSE, stringsAsFactors = FALSE)[,1]



digest::digest(mat, algo = 'md5')
digest::digest(metadata, algo = 'md5')
digest::digest(sig_genes, algo = 'md5')

# first 5 rows; first 5 columns
mat[1:5,1:5]

# take a peek at the metadata
head(metadata)

# take a peek at the genes identified as statistically significant
head(sig_genes)

# dimensions of expression data and metadata, and length of sig_genes
dim(mat)

# verify integrity of metadata and expression matrix:
# --check that both objects are aligned by name
all(rownames(metadata) == colnames(mat))

# Subset the expression matrix for the statistically significant genes
mat <- mat[sig_genes,]

## Generate the heatmap
# First scale the data to Z-scores by row
heat <- t(scale(t(mat)))
head(heat)

# Set colour scheme and choose breaks
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

## Generate color mappings for the metadata histology scores
# CD3
cd3 <- metadata$CD3
cd3 <- cd3[!is.na(cd3)] # remove missing values - we don't want to include these in the mapping
pick.col <- brewer.pal(9, 'Greens')
col.cd3 <- colorRampPalette(pick.col)(length(unique(cd3)))

# CD20
cd20 <- metadata$CD20
cd20 <- cd20[!is.na(cd20)]
pick.col <- brewer.pal(9, 'Blues')
col.cd20 <- colorRampPalette(pick.col)(length(unique(cd20)))

# CD68L
cd68L <- metadata$CD68L
cd68L <- cd68L[!is.na(cd68L)]
pick.col <- brewer.pal(9, 'Reds')
col.cd68L <- colorRampPalette(pick.col)(length(unique(cd68L)))

# CD68SL
cd68SL <- metadata$CD68SL
cd68SL <- cd68L[!is.na(cd68L)]
pick.col <- brewer.pal(9, 'Oranges')
col.cd68SL <- colorRampPalette(pick.col)(length(unique(cd68SL)))

# CD138
cd138 <- metadata$CD138
cd138 <- cd138[!is.na(cd138)]
pick.col <- brewer.pal(9, 'Purples')
col.cd138 <- colorRampPalette(pick.col)(length(unique(cd68SL)))

# Creat inital annotation data frame for the heatmap
ann <-data.frame(
  Pathotype = metadata$Pathotype,
  CD3 = metadata$CD3,
  CD20 = metadata$CD20,
  CD68L = metadata$CD68L,
  CD68SL = metadata$CD68SL,
  CD138 = metadata$CD138,
  stringsAsFactors = FALSE) 

# create the colour mapping
colours <- list(
  Pathotype = c('Lymphoid' = 'blue', 'Myeloid' = 'red', 
                'Fibroid' = 'green3', 'Ungraded' = 'grey'),
  CD3 = c('0' = '#F7FCF5', '1' = '#C7E9C0', '2' = '#74C476', 
          '3' = '#238B45', '4' = '#00441B'),
  CD20 = c('0' = '#F7FBFF', '1' = '#C6DBEF', '2' = '#6BAED6', 
           '3' = '#2171B5', '4' = '#08306B'),
  CD68L = c('0' = '#FFF5F0', '1' = '#FCBBA1', '2' = '#FB6A4A', 
            '3' = '#CB181D', '4' = '#67000D'),
  CD68SL = c('0' = '#FFF5EB', '1' = '#FDD0A2', '2' = '#FD8D3C', 
             '3' = '#D94801', '4' = '#7F2704'),
  CD138 = c('0' = '#FCFBFD', '1' = '#DADAEB', '2' = '#9E9AC8', 
            '3' = '#6A51A3', '4' = '#3F007D'))
colours

# Create ComplexHeatmap Annotation object
colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Pathotype = list(
      nrow = 4, # number of rows across which the legend will be arranged
      title = 'Pathotype',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    CD3 = list(
      nrow = 5,
      title = 'CD3',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    CD20 = list(
      nrow = 5,
      title = 'CD20',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    CD68L = list(
      nrow = 5,
      title = 'CD68L',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    CD68SL = list(
      nrow = 5,
      title = 'CD68SL',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    CD138 = list(
      nrow = 5,
      title = 'CD138',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold'))))

# Create annotation box-and-whisker plots
boxplotCol <- HeatmapAnnotation(
  boxplot = anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'left')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'col')

boxplotRow <- HeatmapAnnotation(
  boxplot = row_anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'top')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'row')

# Create gene labels
# Here, we can ‘step through’ the 
# variables / genes and choose to only label a select few.
# The number of rows (genes) in our object is: 2772
# In this code snippet, we ‘step through’ the rownames and only 
# retain every 40th successive label.
genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = seq(1, nrow(heat), 40),
    labels = rownames(heat)[seq(1, nrow(heat), 40)],
    labels_gp = gpar(fontsize = 10, fontface = 'bold'),
    padding = 0.75),
  width = unit(2.0, 'cm') +
    
    max_text_width(
      rownames(heat)[seq(1, nrow(heat), 40)],
      gp = gpar(fontsize = 10,  fontface = 'bold')))

# Perform partitioning around medoids to identify internal
# 'structure' in the data that may relate to biologically 
# meaningful pathways
pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)

# fix order of the clusters to have 1 to 4, top to bottom
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 
                                            'Cluster 3', 'Cluster 4'))

# Create the actual heatmap object
hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE,
                
                name = 'Gene\nZ-\nscore',
                
                col = colorRamp2(myBreaks, myCol),
                
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                #row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn,
                bottom_annotation = boxplotCol)

# Draw the heatmap
draw(hmap + genelabels,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')


# Change colour scheme, breaks, and do extra clustering on columns
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

hmap1 <- Heatmap(heat,
                 
                 name = 'Gene Z-score',
                 
                 col = colorRamp2(myBreaks, myCol),
                 
                 heatmap_legend_param = list(
                   color_bar = 'continuous',
                   legend_direction = 'horizontal',
                   legend_width = unit(8, 'cm'),
                   legend_height = unit(5.0, 'cm'),
                   title_position = 'topcenter',
                   title_gp=gpar(fontsize = 30, fontface = 'bold'),
                   labels_gp=gpar(fontsize = 24, fontface = 'bold')),
                 
                 cluster_rows = TRUE,
                 show_row_dend = TRUE,
                 row_title = 'Statistically significant genes',
                 row_title_side = 'left',
                 row_title_gp = gpar(fontsize = 30,  fontface = 'bold'),
                 row_title_rot = 90,
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                 row_names_side = 'left',
                 row_dend_width = unit(25,'mm'),
                 
                 cluster_columns = TRUE,
                 show_column_dend = TRUE,
                 column_title = 'Samples',
                 column_title_side = 'bottom',
                 column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
                 column_title_rot = 0,
                 show_column_names = FALSE,
                 column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                 column_names_max_height = unit(10, 'cm'),
                 column_dend_height = unit(25,'mm'),
                 
                 clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                 clustering_method_columns = 'ward.D2',
                 clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                 clustering_method_rows = 'ward.D2')


myCol <- colorRampPalette(c('forestgreen', 'black', 'purple'))(100)
myBreaks <- seq(-2, 2, length.out = 100)

hmap2 <- Heatmap(heat,
                 
                 split = pamClusters$clustering,
                 cluster_row_slices = FALSE,
                 
                 column_km = 6,
                 
                 name = 'Gene Z-score',
                 col = colorRamp2(myBreaks, myCol),
                 
                 heatmap_legend_param = list(
                   color_bar = 'continuous',
                   legend_direction = 'horizontal',
                   legend_width = unit(8, 'cm'),
                   legend_height = unit(5.0, 'cm'),
                   title_position = 'topcenter',
                   title_gp=gpar(fontsize = 30, fontface = 'bold'),
                   labels_gp=gpar(fontsize = 24, fontface = 'bold')),
                 
                 cluster_rows = TRUE,
                 show_row_dend = FALSE,
                 #row_title = 'Statistically significant genes',
                 row_title_side = 'right',
                 row_title_gp = gpar(fontsize = 30,  fontface = 'bold'),
                 row_title_rot = 90,
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 12, fontface = 'bold'),
                 row_names_side = 'left',
                 row_dend_width = unit(25,'mm'),
                 
                 cluster_columns = TRUE,
                 show_column_dend = TRUE,
                 column_title = 'Samples',
                 column_title_side = 'bottom',
                 column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
                 column_title_rot = 0,
                 show_column_names = FALSE,
                 column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                 column_names_max_height = unit(10, 'cm'),
                 column_dend_height = unit(25,'mm'),
                 
                 clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                 clustering_method_columns = 'ward.D2',
                 clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                 clustering_method_rows = 'ward.D2')

pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(hmap1,
     heatmap_legend_side = 'top',
     row_sub_title_side = 'left',
     newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(hmap2,
     heatmap_legend_side = 'top',
     row_sub_title_side = 'right',
     newpage = FALSE)
popViewport()
popViewport()
