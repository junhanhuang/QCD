suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(matrixStats)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(irlba)
  library(princurve)
})

results_base_dir <- file.path(getwd(), "data")
dataset_names <- c(
  "linear_three_clusters_1",
  "multifurcating_four_clusters_1",
  "tree_five_clusters_1",
  "GSE74767",
  "bifurcating_four_clusters_1"
)

start_clusters <- list()

end_clusters <- list()

for (dataset in dataset_names) {
  cat(" 正在处理数据集：", dataset, "\n")
  
  dataset_dir  <- file.path(results_base_dir, dataset)
  expr_file    <- file.path(dataset_dir, "expression.csv")
  cluster_file <- file.path(dataset_dir, "leiden_labels_align.csv")
  
  log_mat <- as.matrix(read.csv(expr_file, row.names = 1))  # cells × genes
  ncells  <- nrow(log_mat)
  ngenes  <- ncol(log_mat)
  cat(sprintf("  表达矩阵维度：细胞数 = %d, 基因数 = %d\n", ncells, ngenes))
  
  cluster_labels <- factor(read.csv(cluster_file, header = FALSE)[,1])
  cluster_labels <- factor(cluster_labels, levels = sort(unique(cluster_labels)))
  
  ndim <- 5
  if (ncol(log_mat) <= ndim) {
    message(sprintf("   PCA维度设置为 %d，但当前特征数为 %d，跳过PCA", ndim, ncol(log_mat)))
    pca_coords <- as.matrix(log_mat)
  } else {
    pca <- irlba::prcomp_irlba(log_mat, n = ndim)
    
    if (ndim > 3) {
      x <- 1:ndim
      opt1 <- which.min(sapply(2:10, function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(pca$sdev[1:ndim] ~ x + x2)$residuals^2 * rep(1:2, each = 10))
      }))
      
      xmat <- cbind(1:ndim, pca$sdev[1:ndim])
      line <- xmat[c(1, nrow(xmat)), ]
      proj <- project_to_curve(xmat, line)
      opt2 <- which.max(proj$dist_ind) - 1
      
      optpoint <- max(c(min(c(opt1, opt2)), 3))
    } else {
      optpoint <- ndim
    }
    
    pca_coords <- pca$x[, seq_len(optpoint)]
    rownames(pca_coords) <- rownames(log_mat)
  }
  
  cat(sprintf("  使用 PCA 主成分数：%d\n", ncol(pca_coords)))
  
  sce <- SingleCellExperiment(
    assays = list(norm = t(log_mat))  
  )
  reducedDims(sce)$PCA <- pca_coords
  colData(sce)$cluster <- cluster_labels
  
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim    = "PCA",
    start.clus    = start_clusters[[dataset]],
    end.clus      = end_clusters[[dataset]],
    dist.method   = "simple",
    reweight      = FALSE,
    reassign      = FALSE,
    smoother      = "smooth.spline",
    shrink.method = "cosine",
    thresh        = 0.0001,
    stretch       = 1,
    allow.breaks  = FALSE
  )
  cat("  Slingshot轨迹推断完成\n")
  
  pal  <- if (length(levels(cluster_labels)) <= 9) {
    brewer.pal(length(levels(cluster_labels)), "Set1")
  } else {
    viridis(length(levels(cluster_labels)))
  }
  cols <- pal[as.integer(cluster_labels)]
  
  png(file.path(dataset_dir, paste0(dataset, "_lineage_structure_leiden.png")),
      width = 2700, height = 2250, res = 300)
  
  plot(
    reducedDims(sce)$PCA[, 1:2],
    col = alpha(cols, 0.7),
    pch = 19, cex = 2.2, asp = 1,
    xlab = "PC1", ylab = "PC2",
    main = paste0(dataset, " | lineage_structure | leiden ")
  )
  lines(SlingshotDataSet(sce), type = "lineages", lwd = 2.5, col = "black")
  legend("topright", legend = levels(cluster_labels), col = pal, pch = 19, pt.cex = 2, bty = "n")
  dev.off()
  
  saveRDS(sce, file = file.path(dataset_dir, "slingshot_sce_leiden.rds"))
  cat("  已保存 SCE 对象和轨迹图：slingshot_sce_leiden.rds\n\n")
}

cat(" 所有数据集处理完成！\n")
