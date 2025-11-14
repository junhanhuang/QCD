suppressPackageStartupMessages({
  library(slingshot)
  library(dynwrap)
  library(dynplot)
  library(tidyverse)
  library(SingleCellExperiment)
  library(purrr)
  library(dyneval)
  library(RColorBrewer)
  library(igraph)
})

method_name="spectral"
label_name="spectral_labels_align.csv" 
results_base_dir  <- file.path(getwd(), "data")
result_output_dir <- file.path(getwd(), "result")
dir.create(result_output_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- c(
  "bifurcating_four_clusters_1",
  "linear_three_clusters_1",
  "multifurcating_four_clusters_1",
  "tree_five_clusters_1",
  "GSE74767"
)

create_trajectory <- function(cells, sce) {
  pseudotime <- slingPseudotime(sce)[cells, , drop = FALSE]
  weights    <- slingCurveWeights(sce)[cells, , drop = FALSE]
  lineages   <- slingLineages(sce)
  all_milestones <- sort(unique(unlist(lineages)))
  milestone_ids  <- paste0("M", all_milestones)
  
  milestone_network <- do.call(rbind, lapply(lineages, function(l) {
    data.frame(
      from     = paste0("M", head(l, -1)),
      to       = paste0("M", tail(l, -1)),
      length   = 1,
      directed = FALSE
    )
  })) %>% distinct() %>% as_tibble()
  
  
  cluster <- slingClusterLabels(sce)[cells, , drop = FALSE]
  lin_assign <- apply(weights, 1, which.max)
  
  progressions <- map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    if (!any(ind)) return(tibble())  
    
    lin <- lineages[[l]]
    pst.full <- pseudotime[, l]
    pst <- pst.full[ind]
    

    means <- sapply(lin, function(clID){
      cluster_weights <- cluster[, clID]
      if (sum(cluster_weights) == 0) return(NA)
      stats::weighted.mean(pst.full, cluster_weights)
    })
    
    non_ends <- means[-c(1, length(means))]
    non_ends <- non_ends[!is.na(non_ends)]
    
    if (length(non_ends) > 0) {
      edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    } else {
      edgeID.l <- rep(1, length(pst))
    }
    
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]
    
    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1
    pct[is.na(pct)] <- 0
    
    tibble(
      cell_id = cells[ind], 
      from = paste0("M", from.l), 
      to = paste0("M", to.l), 
      percentage = pct
    )
  })
  
  # 如果没有进展数据，创建默认进展
  if (nrow(progressions) == 0) {
    progressions <- tibble(
      cell_id = cells[1],
      from = milestone_ids[1],
      to = milestone_ids[min(2, length(milestone_ids))],
      percentage = 0
    )
  }
  
  wrap_data(cell_ids = cells) %>%
    add_trajectory(
      milestone_ids     = milestone_ids,
      milestone_network = milestone_network,
      progressions      = progressions
    ) %>%
    add_cell_waypoints(num_cells_selected = length(cells))
}

results_df <- tibble(
  dataset         = character(),
  F1_branches     = double(),
  F1_milestones   = double(),
  HIM             = double(),
  Correlation     = double(),
  featureimp_wcor = double()
)

for (name in datasets) {
  cat(sprintf("\n处理数据集: %s\n", name))
  dir <- file.path(results_base_dir, name)
  
  tryCatch({
    sce       <- readRDS(file.path(dir, paste0("slingshot_sce_",method_name,".rds")))
    truth_sce <- readRDS(file.path(dir, "slingshot_sce_truth.rds"))
    expr      <- assay(sce, "norm")
    
    grouping <- as.factor(read.csv(file.path(dir, label_name), header = FALSE)[, 1])
    names(grouping) <- colnames(sce)
    
    common_cells <- intersect(colnames(sce), colnames(truth_sce))
    cat(sprintf("  共细胞数量: %d\n", length(common_cells)))
    grouping <- grouping[common_cells]
    
    pred_data <- create_trajectory(common_cells, sce)
    true_data <- create_trajectory(common_cells, truth_sce)
    
    reducedDim <- reducedDims(sce)$PCA[common_cells, 1:2]
    pred_data  <- add_dimred(pred_data, reducedDim)
    true_data  <- add_dimred(true_data, reducedDim)
    
    n_color <- max(3, length(unique(grouping)))
    # palette <- setNames(brewer.pal(n_color, "Set1"), sort(unique(grouping)))
    if (n_color <= 9) {
      palette <- setNames(brewer.pal(n_color, "Set1"), sort(unique(grouping)))
    } else {
      palette <- setNames(viridis(n_color), sort(unique(grouping)))
    }
    plot_obj <- dynplot::plot_dimred(
      pred_data,
      grouping = grouping,
      color_cells = "grouping",
      size_cells = 3,
      border_radius_percentage = 0,
      groups = palette
    ) +
      theme(legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)) 
   #   ggtitle(paste0(name, " | Dynplot (",method_name,"聚类分组)"))
    
    ggsave(
      filename = file.path(result_output_dir,method_name, paste0(name, "_dynplot_",method_name,"_grouping.png")),
      plot = plot_obj,
      width = 6, height = 5, dpi = 300
    )
    
    metrics <- dyneval::calculate_metrics(
      dataset = true_data,
      model   = pred_data,
      metrics = c("correlation", "F1_branches", "F1_milestones", "him", "featureimp_wcor"),
      expression_source = t(expr)
    )
    
    results_df <- results_df %>% add_row(
      dataset         = name,
      F1_branches     = metrics$F1_branches,
      F1_milestones   = metrics$F1_milestones,
      HIM             = metrics$him,
      Correlation     = metrics$correlation,
      featureimp_wcor = metrics$featureimp_wcor
    )
    
    cat(sprintf(
      "  评估指标: F1_branch=%.4f | F1_milestone=%.4f | HIM=%.4f | Corr=%.4f | wcor=%.4f\n",
      metrics$F1_branches, metrics$F1_milestones, metrics$him,
      metrics$correlation, metrics$featureimp_wcor
    ))
    
  }, error = function(e) {
    warning(sprintf("处理数据集 %s 出错: %s", name, conditionMessage(e)))
  })
}

output_file <- file.path(result_output_dir, paste0("trajectory_metrics_summary_",method_name,".csv"))
readr::write_csv(results_df, output_file)
cat("\n有指标已保存至：", output_file, "\n")
