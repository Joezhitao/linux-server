#############################################################
# 单细胞RNA-seq数据细胞簇优化处理脚本 - 模式一
# 功能：循环处理每个seurat_clusters簇，去除散在细胞
# 日期：2023-11-15
#############################################################

# 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(MASS)
  library(ggplot2)
  library(patchwork)
})

# 设置工作目录
setwd("/home/lin/c_group")

# 定义辅助函数 ----------------------------------------------------------
# 判断一个点是否在核心区域内
in_core_region <- function(x, y, kde, threshold) {
  ix <- findInterval(x, kde$x)
  iy <- findInterval(y, kde$y)
  
  if (ix > 0 && ix < length(kde$x) && iy > 0 && iy < length(kde$y)) {
    return(kde$z[ix, iy] > threshold * max(kde$z))
  } else {
    return(FALSE)
  }
}

# 打印处理进度信息
print_progress <- function(message) {
  cat(paste0("\n", message, "\n"))
  cat(paste0(rep("-", nchar(message)), collapse = ""), "\n")
}

# 自定义主题，移除所有坐标轴元素
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

# 主函数：循环处理每个簇并去除散在细胞 ----------------------------------
optimize_all_clusters_mode1 <- function(seurat_obj, 
                                        density_threshold = 0.1,
                                        reduction_type = "umap",
                                        output_dir = "./") {
  
  print_progress("1. 初始化处理")
  
  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 获取所有簇
  all_clusters <- levels(seurat_obj@meta.data$seurat_clusters)
  cat(sprintf("发现 %d 个细胞簇: %s\n", length(all_clusters), paste(all_clusters, collapse = ", ")))
  
  # 设置颜色
  n_colors <- length(all_clusters)
  color_palette <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                                      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
                                      "#bcbd22", "#17becf"))(n_colors)
  
  print_progress("2. 处理每个细胞簇")
  
  # 存储每个簇的处理结果
  cluster_results <- list()
  all_cells_to_keep <- character(0)
  
  # 处理每个细胞簇
  for (cluster in all_clusters) {
    cat(sprintf("\n处理细胞簇: %s\n", cluster))
    
    # 提取降维坐标和簇信息
    df <- seurat_obj@reductions[[reduction_type]]@cell.embeddings %>%
      as.data.frame() %>%
      cbind(cluster = seurat_obj@meta.data$seurat_clusters)
    
    # 提取目标簇的细胞坐标
    dim_cols <- colnames(df)[1:2]  # 获取前两个维度的列名
    target_cells <- df[df$cluster == cluster, dim_cols]
    cat(sprintf("  该簇细胞总数: %d\n", nrow(target_cells)))
    
    # 计算核密度估计
    kde <- kde2d(target_cells[[1]], target_cells[[2]], n = 100)
    
    # 识别核心区域内的目标簇细胞
    target_indices <- which(seurat_obj@meta.data$seurat_clusters == cluster)
    target_coords <- seurat_obj@reductions[[reduction_type]]@cell.embeddings[target_indices, 1:2]
    
    in_core <- sapply(1:nrow(target_coords), function(i) {
      in_core_region(target_coords[i, 1], target_coords[i, 2], kde, density_threshold)
    })
    
    core_target_cells <- rownames(target_coords)[in_core]
    removed_target_cells <- rownames(target_coords)[!in_core]
    
    # 存储结果
    cluster_results[[cluster]] <- list(
      core_cells = core_target_cells,
      removed_cells = removed_target_cells
    )
    
    # 添加到要保留的细胞列表
    all_cells_to_keep <- c(all_cells_to_keep, core_target_cells)
    
    # 输出统计信息
    cat(sprintf("  保留的核心细胞: %d (%.1f%%)\n", 
                length(core_target_cells), 
                100*length(core_target_cells)/nrow(target_cells)))
    
    cat(sprintf("  移除的散在细胞: %d (%.1f%%)\n", 
                length(removed_target_cells), 
                100*length(removed_target_cells)/nrow(target_cells)))
  }
  
  print_progress("3. 创建优化后的Seurat对象")
  
  # 创建优化后的Seurat对象
  seurat_filtered <- subset(seurat_obj, cells = all_cells_to_keep)
  
  # 可视化结果 ----------------------------------------------------------
  print_progress("4. 可视化结果")
  
  # 降维可视化 (UMAP或tSNE)
  p_orig <- DimPlot(seurat_obj, reduction = reduction_type, group.by = "seurat_clusters", pt.size = 1,
                    label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle(paste0("原始数据 (", toupper(reduction_type), ")")) + 
    scale_color_manual(values = color_palette) +
    custom_theme()
  
  p_filtered <- DimPlot(seurat_filtered, reduction = reduction_type, group.by = "seurat_clusters", pt.size = 1,
                        label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle(paste0("模式一：移除散在细胞 (", toupper(reduction_type), ")")) + 
    scale_color_manual(values = color_palette) +
    custom_theme()
  
  # 组合图表
  combined_plot <- p_orig + p_filtered + 
    plot_layout(ncol = 2, guides = "collect")
  
  print(combined_plot)
  
  # 保存图表
  ggsave(file.path(output_dir, paste0("original_vs_filtered_", reduction_type, ".pdf")), 
         combined_plot, width = 12, height = 6, dpi = 300)
  
  # 标记过滤状态 --------------------------------------------------------
  print_progress("5. 标记过滤状态")
  
  # 添加过滤状态到元数据
  seurat_obj@meta.data$filter_status <- "保留"
  
  # 标记移除的细胞
  for (cluster in all_clusters) {
    removed_cells <- cluster_results[[cluster]]$removed_cells
    seurat_obj@meta.data$filter_status[rownames(seurat_obj@meta.data) %in% removed_cells] <- 
      paste0("移除的簇", cluster, "散在细胞")
  }
  
  # 可视化不同的移除类别
  filter_colors <- c("保留" = "grey")
  for (i in seq_along(all_clusters)) {
    cluster <- all_clusters[i]
    filter_colors[paste0("移除的簇", cluster, "散在细胞")] <- color_palette[i]
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
  saveRDS(seurat_filtered, file.path(output_dir, paste0("filtered_all_clusters_threshold", density_threshold, ".rds")))
  
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
pbmc <- readRDS("/home/lin/c_group/C_group_harmony.rds")

# 检查数据
print_progress("数据概况")
cat(sprintf("细胞总数: %d\n", ncol(pbmc)))
cat(sprintf("基因总数: %d\n", nrow(pbmc)))
cat(sprintf("细胞簇数量: %d\n", length(levels(pbmc@meta.data$seurat_clusters))))
print(table(pbmc@meta.data$seurat_clusters))

# 设置密度阈值
density_threshold <- 0.2

# 运行优化处理（循环处理所有簇）- UMAP
print_progress("运行UMAP优化处理")
results_umap <- optimize_all_clusters_mode1(
  seurat_obj = pbmc,
  density_threshold = density_threshold,
  reduction_type = "umap",
  output_dir = "/home/lin/c_group/optimization_results/umap"
)

# 运行优化处理（循环处理所有簇）- tSNE
print_progress("运行tSNE优化处理")
results_tsne <- optimize_all_clusters_mode1(
  seurat_obj = pbmc,
  density_threshold = density_threshold,
  reduction_type = "tsne",
  output_dir = "/home/lin/c_group/optimization_results/tsne"
)

# 查看结果统计
print_progress("处理完成!")
cat(sprintf("原始细胞数: %d\n", ncol(results_umap$original)))
cat(sprintf("UMAP过滤后细胞数: %d (保留%.1f%%)\n", 
            ncol(results_umap$filtered), 
            100*ncol(results_umap$filtered)/ncol(results_umap$original)))
cat(sprintf("tSNE过滤后细胞数: %d (保留%.1f%%)\n", 
            ncol(results_tsne$filtered), 
            100*ncol(results_tsne$filtered)/ncol(results_tsne$original)))

# 查看每个簇的过滤情况 (UMAP)
print_progress("UMAP过滤统计")
for (cluster in names(results_umap$cluster_results)) {
  total_cells <- length(results_umap$cluster_results[[cluster]]$core_cells) + 
    length(results_umap$cluster_results[[cluster]]$removed_cells)
  
  cat(sprintf("\n簇 %s:\n", cluster))
  cat(sprintf("  总细胞数: %d\n", total_cells))
  cat(sprintf("  保留细胞数: %d (%.1f%%)\n", 
              length(results_umap$cluster_results[[cluster]]$core_cells),
              100*length(results_umap$cluster_results[[cluster]]$core_cells)/total_cells))
  cat(sprintf("  移除细胞数: %d (%.1f%%)\n", 
              length(results_umap$cluster_results[[cluster]]$removed_cells),
              100*length(results_umap$cluster_results[[cluster]]$removed_cells)/total_cells))
}
