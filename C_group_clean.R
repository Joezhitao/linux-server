#############################################################
# 单细胞RNA-seq数据特定细胞簇优化处理脚本
# 功能：删除指定细胞类型中远离主体的小簇
# 日期：2023-11-15
#############################################################

# 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stats) # 用于kmeans聚类
})

# 打印处理进度信息
print_progress <- function(message) {
  cat(paste0("\n", message, "\n"))
  cat(paste0(rep("-", nchar(message)), collapse = ""), "\n")
}

# 自定义主题
custom_theme <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
}

# 主函数：处理指定细胞类型并删除远离主体的小簇 ---------------------------
remove_distant_subclusters <- function(seurat_obj, 
                                       target_clusters,
                                       reduction_type = "umap",
                                       output_dir = "./") {
  
  print_progress("1. 初始化处理")
  
  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 确认目标簇存在
  all_clusters <- levels(seurat_obj@meta.data$seurat_clusters)
  valid_targets <- target_clusters[target_clusters %in% all_clusters]
  
  if (length(valid_targets) == 0) {
    stop("未找到指定的目标簇!")
  }
  
  cat(sprintf("将处理 %d 个指定细胞簇: %s\n", length(valid_targets), paste(valid_targets, collapse = ", ")))
  
  # 设置颜色
  color_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  
  print_progress("2. 处理指定细胞簇")
  
  # 存储每个簇的处理结果
  cluster_results <- list()
  cells_to_keep <- rownames(seurat_obj@meta.data)
  
  # 处理每个指定的细胞簇
  for (cluster in valid_targets) {
    cat(sprintf("\n处理细胞簇: %s\n", cluster))
    
    # 提取目标簇的细胞坐标
    target_indices <- which(seurat_obj@meta.data$seurat_clusters == cluster)
    target_coords <- seurat_obj@reductions[[reduction_type]]@cell.embeddings[target_indices, 1:2]
    cat(sprintf("  该簇细胞总数: %d\n", nrow(target_coords)))
    
    # 使用k-means将目标簇分成两个子簇
    set.seed(42) # 设置随机种子以确保结果可重复
    kmeans_result <- kmeans(target_coords, centers = 2)
    
    # 计算每个子簇的细胞数量
    subcluster_sizes <- table(kmeans_result$cluster)
    cat(sprintf("  子簇1细胞数: %d\n", subcluster_sizes[1]))
    cat(sprintf("  子簇2细胞数: %d\n", subcluster_sizes[2]))
    
    # 确定哪个是较小的子簇
    smaller_subcluster <- which.min(subcluster_sizes)
    
    # 计算两个子簇中心之间的距离
    centers_distance <- dist(kmeans_result$centers)[1]
    cat(sprintf("  两个子簇中心之间的距离: %.2f\n", centers_distance))
    
    # 找出属于较小子簇的细胞
    small_subcluster_cells <- rownames(target_coords)[kmeans_result$cluster == smaller_subcluster]
    
    # 存储结果
    cluster_results[[cluster]] <- list(
      main_subcluster_cells = rownames(target_coords)[kmeans_result$cluster != smaller_subcluster],
      small_subcluster_cells = small_subcluster_cells
    )
    
    # 更新要保留的细胞列表（移除较小子簇的细胞）
    cells_to_keep <- setdiff(cells_to_keep, small_subcluster_cells)
    
    # 输出统计信息
    cat(sprintf("  保留的主体簇细胞: %d (%.1f%%)\n", 
                length(cluster_results[[cluster]]$main_subcluster_cells), 
                100*length(cluster_results[[cluster]]$main_subcluster_cells)/nrow(target_coords)))
    
    cat(sprintf("  移除的远离主体的小簇细胞: %d (%.1f%%)\n", 
                length(small_subcluster_cells), 
                100*length(small_subcluster_cells)/nrow(target_coords)))
  }
  
  print_progress("3. 创建优化后的Seurat对象")
  
  # 创建优化后的Seurat对象
  seurat_filtered <- subset(seurat_obj, cells = cells_to_keep)
  
  # 可视化结果 ----------------------------------------------------------
  print_progress("4. 可视化结果")
  
  # 降维可视化
  p_orig <- DimPlot(seurat_obj, reduction = reduction_type, group.by = "seurat_clusters", pt.size = 1,
                    label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle(paste0("原始数据 (", toupper(reduction_type), ")")) + 
    custom_theme()
  
  p_filtered <- DimPlot(seurat_filtered, reduction = reduction_type, group.by = "seurat_clusters", pt.size = 1,
                        label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle(paste0("优化后: 移除远离主体的小簇")) + 
    custom_theme()
  
  # 组合图表
  combined_plot <- p_orig + p_filtered + 
    plot_layout(ncol = 2, guides = "collect")
  
  print(combined_plot)
  
  # 保存图表
  ggsave(file.path(output_dir, paste0("removed_distant_subclusters_", reduction_type, ".pdf")), 
         combined_plot, width = 12, height = 6, dpi = 300)
  
  # 标记过滤状态 --------------------------------------------------------
  print_progress("5. 标记过滤状态")
  
  # 添加过滤状态到元数据
  seurat_obj@meta.data$filter_status <- "保留"
  
  # 标记移除的细胞
  for (cluster in valid_targets) {
    removed_cells <- cluster_results[[cluster]]$small_subcluster_cells
    seurat_obj@meta.data$filter_status[rownames(seurat_obj@meta.data) %in% removed_cells] <- 
      paste0("移除的", cluster, "远离主体小簇")
  }
  
  # 可视化过滤状态
  filter_colors <- c("保留" = "grey")
  for (i in seq_along(valid_targets)) {
    filter_colors[paste0("移除的", valid_targets[i], "远离主体小簇")] <- color_palette[i]
  }
  
  filter_plot <- DimPlot(seurat_obj, reduction = reduction_type, group.by = "filter_status") + 
    scale_color_manual(values = filter_colors) +
    ggtitle("细胞过滤状态") +
    custom_theme()
  
  print(filter_plot)
  
  # 保存过滤状态图表
  ggsave(file.path(output_dir, paste0("filter_status_", reduction_type, ".pdf")), 
         filter_plot, width = 10, height = 8, dpi = 300)
  
  # 保存结果 ------------------------------------------------------------
  print_progress("6. 保存结果")
  
  # 保存过滤后的Seurat对象
  saveRDS(seurat_filtered, file.path(output_dir, paste0("filtered_removed_distant_subclusters.rds")))
  
  # 返回结果
  return(list(
    original = seurat_obj,
    filtered = seurat_filtered,
    cluster_results = cluster_results
  ))
}

# 使用示例 --------------------------------------------------------------
# 加载数据
print_progress("加载数据")
pbmc <- readRDS("/home/lin/c_group/ident.rds")

# 检查数据
print_progress("数据概况")
cat(sprintf("细胞总数: %d\n", ncol(pbmc)))
cat(sprintf("基因总数: %d\n", nrow(pbmc)))
print(table(pbmc@meta.data$seurat_clusters))

# 设置参数
target_clusters <- c("Hepatic Progenitor Cells", "Neutrophils")  # 指定要处理的细胞类型

# 运行优化处理（删除远离主体的小簇）
print_progress("运行优化处理")
results <- remove_distant_subclusters(
  seurat_obj = pbmc,
  target_clusters = target_clusters,
  reduction_type = "umap",
  output_dir = "/home/lin/c_group/optimization_results/distant_subclusters"
)

# 查看结果统计
print_progress("处理完成!")
cat(sprintf("原始细胞数: %d\n", ncol(results$original)))
cat(sprintf("过滤后细胞数: %d (保留%.1f%%)\n", 
            ncol(results$filtered), 
            100*ncol(results$filtered)/ncol(results$original)))

# 查看每个处理簇的过滤情况
print_progress("过滤统计")
for (cluster in names(results$cluster_results)) {
  total_cells <- length(results$cluster_results[[cluster]]$main_subcluster_cells) + 
    length(results$cluster_results[[cluster]]$small_subcluster_cells)
  
  cat(sprintf("\n簇 %s:\n", cluster))
  cat(sprintf("  总细胞数: %d\n", total_cells))
  cat(sprintf("  保留的主体簇细胞: %d (%.1f%%)\n", 
              length(results$cluster_results[[cluster]]$main_subcluster_cells),
              100*length(results$cluster_results[[cluster]]$main_subcluster_cells)/total_cells))
  cat(sprintf("  移除的远离主体小簇细胞: %d (%.1f%%)\n", 
              length(results$cluster_results[[cluster]]$small_subcluster_cells),
              100*length(results$cluster_results[[cluster]]$small_subcluster_cells)/total_cells))
}
