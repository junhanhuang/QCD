library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)

results_base_dir <- file.path(getwd(), "data")

datasets <- c(
  "bifurcating_four_clusters_1",
  "linear_three_clusters_1",
  "multifurcating_four_clusters_1",
  "tree_five_clusters_1",
  "GSE74767"
)

for (name in datasets) {
  cat("Processing:", name, "\n")
  
  dir <- file.path(results_base_dir, name)
  umap_file <- file.path(dir, "umap_coords.csv")
  label_file <- file.path(dir, "leiden_labels.csv")
  
  umap <- as.matrix(read.csv(umap_file, header = FALSE))
  labels <- as.factor(read.csv(label_file, header = FALSE)[, 1])
  n_cells <- nrow(umap)
  n_clusters <- length(levels(labels))
  
  dummy_expr <- matrix(0, nrow = 1, ncol = n_cells)
  sce <- SingleCellExperiment(assays = list(counts = dummy_expr))
  reducedDims(sce)$UMAP <- umap
  colData(sce)$cluster <- labels
  
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim = "UMAP",
    smoother = "smooth.spline"
  )

  colorset <- if (n_clusters <= 9) {
    RColorBrewer::brewer.pal(n_clusters, "Set1")
  } else {
    viridis::viridis(n_clusters)
  }
  cell_colors <- colorset[as.integer(labels)]
  
  out_img <- file.path(dir, "slingshot_result.png")
  png(out_img, width = 1200, height = 1000, res = 150)
  
  plot(
    reducedDims(sce)$UMAP,
    col = adjustcolor(cell_colors, alpha.f = 0.8),
    pch = 19,
    cex = 2,
    asp = 1,
    xlab = "UMAP1", ylab = "UMAP2",
    main = paste0(name, " | Slingshot Trajectory"),
    bty = "l"
  )
  lines(SlingshotDataSet(sce), lwd = 5, col = "#222222")
  legend(
    "topright",
    legend = levels(labels),
    col = colorset[1:n_clusters],
    pch = 19,
    pt.cex = 2,
    bty = "n"
  )

  dev.off()
  
  cat("Saved:", out_img, "\n\n")
}
