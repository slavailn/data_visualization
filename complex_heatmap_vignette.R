library(ComplexHeatmap)
library(circlize)
library(gridtext)

## Create a single heatmap
# First, let's generate random matrix where there are three groups in the
# columns and three groups in the rows
set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))
head(mat)

# Draw a heatmap
Heatmap(mat)

# You need to use draw() finction explictly to draw heatmaps
# inside of functions
#for(...) {
#  ht = Heatmap(mat)
#  draw(ht)
#}

# Set color pallete and re-draw a map
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(mat, name = "mat", col = col_fun)

# colorRamp2() makes colors in multiple heatmaps comparable
# if they are set with a same color mapping function
Heatmap(mat, name = "mat", col = col_fun, column_title = "mat")
Heatmap(mat/4, name = "mat", col = col_fun, column_title = "mat/4")
Heatmap(abs(mat), name = "mat", col = col_fun, column_title = "abs(mat)")

# If the matrix is continuous, you can also simply provide a vector of 
# colors and colors will be linearly interpolated. But remember 
# this method is not robust to outliers because the mapping starts 
# from the minimal value in the matrix and ends with the maximal value.
Heatmap(mat, name = "mat", col = rev(rainbow(10)), 
        column_title = "set a color vector for a continuous matrix")

# In the following example, we set colors for a discrete numeric matrix. 
# You don’t need to convert it to a character matrix, 
# just set the numbers as “names” of the color vector.
discrete_mat = matrix(sample(1:4, 100, replace = TRUE), 10, 10)
colors = structure(1:4, names = c("1", "2", "3", "4")) # black, red, green, blue
Heatmap(discrete_mat, name = "mat", col = colors,
        column_title = "a discrete numeric matrix")

# Discrete character matrix
discrete_mat = matrix(sample(letters[1:4], 100, replace = TRUE), 10, 10)
colors = structure(1:4, names = letters[1:4])
Heatmap(discrete_mat, name = "mat", col = colors,
        column_title = "a discrete character matrix")

# Heatmap allows NA values
mat_with_na = mat
na_index = sample(c(TRUE, FALSE), nrow(mat)*ncol(mat), replace = TRUE, prob = c(1, 9))
mat_with_na[na_index] = NA
Heatmap(mat_with_na, name = "mat", na_col = "black",
        column_title = "a matrix with NA values")

# Specify color space, for example LAB or RGB
f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"))
f2 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"), 
                space = "RGB")
Heatmap(mat, name = "mat1", col = f1, column_title = "LAB color space")
Heatmap(mat, name = "mat2", col = f2, column_title = "RGB color space")

# Set heatmap borders
Heatmap(mat, name = "mat", border_gp = gpar(col = "black", lty = 2),
        column_title = "set heatmap borders")

# Set cell borders
Heatmap(mat, name = "mat", rect_gp = gpar(col = "white", lwd = 2),
        column_title = "set cell borders")

# Set row and column titles
Heatmap(mat, name = "mat", column_title = "I am a column title", 
        row_title = "I am a row title")
Heatmap(mat, name = "mat", column_title = "I am a column title at the bottom", 
        column_title_side = "bottom")

# Change title fontsize
Heatmap(mat, name = "mat", column_title = "I am a big column title", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))

# Rotate titles
Heatmap(mat, name = "mat", row_title = "row title", row_title_rot = 0)

# Set background colors for the titles
Heatmap(mat, name = "mat", column_title = "I am a column title", 
        column_title_gp = gpar(fill = "red", col = "white", border = "blue"))

# There is not enough space at the top of the column title when the background
# is colored. We can add padding to fix this problem.
ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
Heatmap(mat, name = "mat", column_title = "I am a column title", 
        column_title_gp = gpar(fill = "red", col = "white", border = "blue"))

# Title can be set a mathematical formula
Heatmap(mat, name = "mat", 
        column_title = expression(hat(beta) == (X^t * X)^{-1} * X^t * y)) 

# More complicated text can be drawn using gridtext package
Heatmap(mat, name = "mat",
        column_title = gt_render(
          paste0("Some <span style='color:blue'>blue text **in bold.**</span><br>",
                 "And *italics text.*<br>And some ",
                 "<span style='font-size:18pt; color:black'>large</span> text."), 
          r = unit(2, "pt"), 
          padding = unit(c(2, 2, 2, 2), "pt")
        )
)

############################## CLUSTERING ################################
# Turn off clustering
Heatmap(mat, name = "mat", cluster_rows = FALSE)

# Cluster rows and columns
Heatmap(mat, name = "mat", row_dend_side = "right", 
        column_dend_side = "bottom")

# Change position and appearance of dendrograms
Heatmap(mat, name = "mat", column_dend_height = unit(4, "cm"), 
        row_dend_width = unit(4, "cm"))

# Clustering based on Pearson correlation
Heatmap(mat, name = "mat", clustering_distance_rows = "pearson",
        column_title = "pre-defined distance method (1 - pearson)")

# Clustering based on pair-wise distances, robust to outliers
mat_with_outliers = mat
for(i in  1:10) mat_with_outliers[i, i] = 1000
robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

Heatmap(mat_with_outliers, name = "mat", 
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        clustering_distance_rows = robust_dist,
        clustering_distance_columns = robust_dist,
        column_title = "robust_dist")

# Specify clustering method, for example "single"
Heatmap(mat, name = "mat", clustering_method_rows = "single")

# Note that if you are specifying clustering method, you need to 
# transpose the matrix
library(cluster)
Heatmap(mat, name = "mat", cluster_rows = diana(mat),
        cluster_columns = agnes(t(mat)), column_title = "clustering objects")

# if cluster_columns is set as a function, you don't need to transpose the matrix
Heatmap(mat, name = "mat", cluster_rows = diana,
        cluster_columns = agnes, column_title = "clustering functions")

# Render dendrograms with dendextend method
library(dendextend)
row_dend = as.dendrogram(hclust(dist(mat)))
row_dend = color_branches(row_dend, k = 2) # `color_branches()` 
# returns a dendrogram object
Heatmap(mat, name = "mat", cluster_rows = row_dend)

# Add graphics to the nodes of the dendrogram
row_dend = dendrapply(row_dend, function(d) {
  attr(d, "nodePar") = list(cex = 0.8, pch = sample(20, 1), col = rand_color(1))
  return(d)
})
Heatmap(mat, name = "mat", cluster_rows = row_dend, row_dend_width = unit(2, "cm"))

# Changing the way dendrograms are reordered, use dendsort package
# Let's compare default reordering and dendsort method
Heatmap(mat, name = "mat", column_title = "default reordering")

library(dendsort)
row_dend = dendsort(hclust(dist(mat)))
col_dend = dendsort(hclust(dist(t(mat))))
Heatmap(mat, name = "mat", cluster_rows = row_dend, cluster_columns = col_dend,
        column_title = "reorder by dendsort")

# Apply seriation to matrix ordering
library(seriation)
o = seriate(max(mat) - mat, method = "BEA_TSP")
Heatmap(max(mat) - mat, name = "mat", 
        row_order = get_order(o, 1), column_order = get_order(o, 2),
        column_title = "seriation by BEA_TSP method")

# Dimension labels
Heatmap(mat, name = "mat", row_names_side = "left", row_dend_side = "right", 
        column_names_side = "top", column_dend_side = "bottom")

# Remove row names
Heatmap(mat, name = "mat", show_row_names = FALSE)

# Change fontsize
Heatmap(mat, name = "mat", row_names_gp = gpar(fontsize = 20))

# Change colors of the text labels
Heatmap(mat, name = "mat", 
        row_names_gp = gpar(col = c(rep("red", 10), rep("blue", 8))))

# Rotate labels
Heatmap(mat, name = "mat", column_names_rot = 45)
Heatmap(mat, name = "mat", column_names_rot = 45, column_names_side = "top",
        column_dend_side = "bottom")

# Set maximum space for labels that are too long
mat2 = mat
rownames(mat2)[1] = paste(c(letters, LETTERS), collapse = "")
Heatmap(mat2, name = "mat", row_title = "default row_names_max_width")
Heatmap(mat2, name = "mat", row_title = "row_names_max_width as length of a*",
        row_names_max_width = max_text_width(
          rownames(mat2), 
          gp = gpar(fontsize = 12)
        ))

############################ HEATMAP SPLIT ##################################
## Split by k means clustering
Heatmap(mat, name = "mat", row_km = 2)
Heatmap(mat, name = "mat", column_km = 3)
Heatmap(mat, name = "mat", row_km = 2, column_km = 3)

# Run consensus k-means clustering 
Heatmap(mat, name = "mat", 
        row_km = 2, row_km_repeats = 100,
        column_km = 3, column_km_repeats = 100)

# Split by categorical variables
Heatmap(mat, name = "mat", 
        row_split = rep(c("A", "B"), 9), column_split = rep(c("C", "D"), 12))

# Split by data frame
Heatmap(mat, name = "mat", 
        row_split = data.frame(rep(c("A", "B"), 9), rep(c("C", "D"), 
                                                        each = 9)))

# Split on both dimensions
Heatmap(mat, name = "mat", row_split = factor(rep(c("A", "B"), 9)),
        column_split = factor(rep(c("C", "D"), 12)))

# Split using different method
pa = cluster::pam(mat, k = 3)
Heatmap(mat, name = "mat", row_split = paste0("pam", pa$clustering))

# Split by dendrogram
Heatmap(mat, name = "mat", row_split = 2, column_split = 3)

# Hclust dendrogram
dend = as.dendrogram(hclust(dist(mat)))
dend = color_branches(dend, k = 2)
Heatmap(mat, name = "mat", cluster_rows = dend, row_split = 2)

# Split heatmap annotations
Heatmap(mat, name = "mat", row_km = 2, column_km = 3,
        top_annotation = HeatmapAnnotation(foo1 = 1:24, bar1 = anno_points(runif(24))),
        right_annotation = rowAnnotation(foo2 = 18:1, bar2 = anno_barplot(runif(18)))
)

########################## Heatmap annotations ###############################
set.seed(123)
mat = matrix(rnorm(100), 10)
rownames(mat) = paste0("R", 1:10)
colnames(mat) = paste0("C", 1:10)
column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)

# Assign bottom and left annotations
Heatmap(mat, name = "mat", bottom_annotation = column_ha, 
        left_annotation = row_ha)

# Simple annotation
ha <- HeatmapAnnotation(foo = 1:10)

# Discrete annotation
ha <-
  
  HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE))
