#############################################################
# 单细胞RNA-seq数据细胞簇优化处理脚本 (精简版)
# 功能：去除散在细胞，优化聚类结果
# 日期：2025-10-21
#############################################################
rm(list = ls())

# 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(MASS)
  library(ggsci)
  library(ggplot2)
  library(patchwork)
})

# 判断细胞是否在核心区域内
in_core_region <- function(x, y, kde, threshold) {
  x_idx <- which.min(abs(kde$x - x))
  y_idx <- which.min(abs(kde$y - y))
  density_value <- kde$z[x_idx, y_idx]
  return(density_value >= threshold * max(kde$z))
}

# 主函数：细胞簇优化 ----------------------------------------------------
optimize_clusters <- function(seurat_obj, 
                              target_clusters = NULL,
                              cluster_column = "seurat_clusters",
                              reduction = "umap",
                              density_threshold = 0.1,
                              filter_mode = "both",
                              output_dir = "./cluster_opt") {
  
  # 验证参数
  if (!reduction %in% names(seurat_obj@reductions))
    stop(paste0("错误: 对象中不存在'", reduction, "'降维。可用的降维有: ", 
                paste(names(seurat_obj@reductions), collapse=", ")))
  
  if (!cluster_column %in% colnames(seurat_obj@meta.data))
    stop(paste0("错误: 元数据中不存在'", cluster_column, "'列。"))
  
  # 如果未指定目标簇，则使用所有簇
  if (is.null(target_clusters))
    target_clusters <- unique(seurat_obj@meta.data[[cluster_column]])
  
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 提取降维坐标和细胞类型信息
  coords <- seurat_obj@reductions[[reduction]]@cell.embeddings %>%
    as.data.frame() %>%
    setNames(c("dim1", "dim2"))
  
  # 设置颜色
  colors <- c(pal_d3("category20")(20), pal_d3("category20b")(20), 
              pal_d3("category20c")(20), pal_d3("category10")(10))
  
  # 存储结果
  cluster_results <- list()
  core_cells_all <- character(0)
  scattered_cells_all <- character(0)
  infiltrating_cells_all <- character(0)
  
  # 创建信息摘要
  summary_info <- c(
    "# 单细胞RNA-seq数据细胞簇优化处理结果摘要",
    "=================================================",
    paste("处理日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("使用降维方法:", reduction),
    paste("聚类列:", cluster_column),
    paste("密度阈值:", density_threshold),
    paste("过滤模式:", filter_mode),
    paste("目标簇:", paste(target_clusters, collapse=", ")),
    "=================================================",
    ""
  )
  
  # 处理每个目标细胞簇
  for (cluster in target_clusters) {
    cat(sprintf("处理簇: %s\n", cluster))
    
    # 提取目标簇的细胞坐标
    cluster_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cluster_column]] == cluster]
    if (length(cluster_cells) < 10) {
      cat(sprintf("  簇 '%s' 细胞数量太少 (%d)，跳过\n", cluster, length(cluster_cells)))
      summary_info <- c(summary_info, 
                        paste("簇", cluster, "细胞数量太少 (", length(cluster_cells), ")，已跳过"))
      next
    }
    
    cluster_coords <- coords[cluster_cells, ]
    
    # 计算核密度估计
    kde <- try(kde2d(cluster_coords$dim1, cluster_coords$dim2, n = 100), silent = TRUE)
    if (inherits(kde, "try-error")) {
      cat(sprintf("  无法为簇 '%s' 计算核密度估计，跳过\n", cluster))
      summary_info <- c(summary_info, 
                        paste("簇", cluster, "无法计算核密度估计，已跳过"))
      next
    }
    
    # 识别核心区域内的目标簇细胞
    in_core <- sapply(1:nrow(cluster_coords), function(i) {
      in_core_region(cluster_coords$dim1[i], cluster_coords$dim2[i], kde, density_threshold)
    })
    
    core_cells <- cluster_cells[in_core]
    scattered_cells <- cluster_cells[!in_core]
    core_cells_all <- c(core_cells_all, core_cells)
    scattered_cells_all <- c(scattered_cells_all, scattered_cells)
    
    # 识别渗入的非目标簇细胞
    infiltrating_cells <- character(0)
    if (filter_mode %in% c("infiltrating", "both")) {
      non_cluster_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cluster_column]] != cluster]
      non_cluster_coords <- coords[non_cluster_cells, ]
      
      in_region <- sapply(1:nrow(non_cluster_coords), function(i) {
        in_core_region(non_cluster_coords$dim1[i], non_cluster_coords$dim2[i], kde, density_threshold)
      })
      
      infiltrating_cells <- non_cluster_cells[in_region]
      infiltrating_cells_all <- c(infiltrating_cells_all, infiltrating_cells)
    }
    
    # 存储结果
    cluster_results[[cluster]] <- list(
      core_cells = core_cells,
      scattered_cells = scattered_cells,
      infiltrating_cells = infiltrating_cells
    )
    
    # 输出统计信息
    scattered_pct <- 100 * length(scattered_cells) / length(cluster_cells)
    infiltrating_pct <- if(length(infiltrating_cells) > 0) {
      100 * length(infiltrating_cells) / length(non_cluster_cells)
    } else {
      0
    }
    
    cat(sprintf("  核心细胞: %d (%.1f%%)\n", 
                length(core_cells), 
                100 - scattered_pct))
    
    cat(sprintf("  散在细胞: %d (%.1f%%)\n", 
                length(scattered_cells), scattered_pct))
    
    cat(sprintf("  渗入细胞: %d\n", length(infiltrating_cells)))
    
    # 添加到摘要
    summary_info <- c(summary_info,
                      paste("簇:", cluster),
                      paste("  总细胞数:", length(cluster_cells)),
                      paste("  核心细胞:", length(core_cells), sprintf("(%.1f%%)", 100 - scattered_pct)),
                      paste("  散在细胞:", length(scattered_cells), sprintf("(%.1f%%)", scattered_pct)),
                      paste("  渗入细胞:", length(infiltrating_cells)),
                      "")
  }
  
  # 创建过滤后的对象
  cells_to_keep <- switch(
    filter_mode,
    "scattered" = {
      # 仅去除散在细胞
      non_target_cells <- rownames(seurat_obj@meta.data)[
        !seurat_obj@meta.data[[cluster_column]] %in% target_clusters
      ]
      c(core_cells_all, non_target_cells)
    },
    "infiltrating" = {
      # 仅去除渗入细胞
      setdiff(colnames(seurat_obj), infiltrating_cells_all)
    },
    "both" = {
      # 同时去除散在细胞和渗入细胞
      non_target_cells <- rownames(seurat_obj@meta.data)[
        !seurat_obj@meta.data[[cluster_column]] %in% target_clusters
      ]
      non_infiltrating <- setdiff(non_target_cells, infiltrating_cells_all)
      c(core_cells_all, non_infiltrating)
    }
  )
  
  seurat_filtered <- subset(seurat_obj, cells = cells_to_keep)
  
  # 添加总体过滤结果到摘要
  cells_removed <- setdiff(colnames(seurat_obj), cells_to_keep)
  summary_info <- c(summary_info,
                    "总体过滤结果",
                    "=================================================",
                    paste("原始细胞总数:", ncol(seurat_obj)),
                    paste("过滤后细胞数:", ncol(seurat_filtered)),
                    paste("移除细胞总数:", length(cells_removed), 
                          sprintf("(%.1f%%)", 100*length(cells_removed)/ncol(seurat_obj))),
                    paste("  其中散在细胞:", length(scattered_cells_all)),
                    paste("  其中渗入细胞:", length(infiltrating_cells_all)),
                    "=================================================")
  
  # 可视化结果
  p_orig <- DimPlot(seurat_obj, reduction = reduction, group.by = cluster_column, 
                    pt.size = 1, label = TRUE, repel = TRUE) +
    ggtitle("原始数据") + 
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p_filtered <- DimPlot(seurat_filtered, reduction = reduction, group.by = cluster_column, 
                        pt.size = 1, label = TRUE, repel = TRUE) +
    ggtitle(paste0("过滤后 (", filter_mode, ", 阈值:", density_threshold, ")")) + 
    scale_color_manual(values = colors) +
    theme_minimal()
  
  combined_plot <- p_orig + p_filtered
  print(combined_plot)
  
  # 创建过滤状态列
  seurat_obj$filter_status <- "保留"
  
  # 根据过滤模式标记细胞
  if (filter_mode %in% c("scattered", "both")) {
    seurat_obj$filter_status[rownames(seurat_obj@meta.data) %in% scattered_cells_all] <- "散在细胞"
  }
  
  if (filter_mode %in% c("infiltrating", "both")) {
    seurat_obj$filter_status[rownames(seurat_obj@meta.data) %in% infiltrating_cells_all] <- "渗入细胞"
  }
  
  # 可视化过滤状态
  filter_colors <- c("保留" = "grey", "散在细胞" = "red", "渗入细胞" = "blue")
  
  p_filter_status <- DimPlot(seurat_obj, reduction = reduction, group.by = "filter_status", 
                             pt.size = 1) +
    ggtitle("细胞过滤状态") + 
    scale_color_manual(values = filter_colors) +
    theme_minimal()
  
  print(p_filter_status)
  
  # 保存结果
  ggsave(file.path(output_dir, paste0("filtered_", filter_mode, "_", 
                                      reduction, "_threshold", density_threshold, ".pdf")), 
         combined_plot, width = 12, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, paste0("filter_status_", filter_mode, "_", 
                                      reduction, "_threshold", density_threshold, ".pdf")), 
         p_filter_status, width = 8, height = 7, dpi = 300)
  
  saveRDS(seurat_filtered, file.path(output_dir, paste0("filtered_", filter_mode, 
                                                        "_threshold", density_threshold, ".rds")))
  
  # 保存信息摘要到文本文件
  writeLines(summary_info, file.path(output_dir, paste0("filter_summary_", filter_mode, 
                                                        "_threshold", density_threshold, ".txt")))
  
  # 返回结果
  return(list(
    filtered = seurat_filtered,
    original = seurat_obj,
    core_cells = core_cells_all,
    scattered_cells = scattered_cells_all,
    infiltrating_cells = infiltrating_cells_all,
    filter_mode = filter_mode,
    summary = summary_info
  ))
}

# 交互式参数设置函数
run_interactive <- function() {
  cat("\n========== 单细胞RNA-seq数据细胞簇优化 ==========\n\n")
  
  # 数据文件路径
  cat("请输入Seurat对象RDS文件路径: ")
  data_path <- readline()
  seurat_obj <- readRDS(data_path)
  
  # 显示可用的元数据列
  cat("\n可用的元数据列:\n")
  meta_cols <- colnames(seurat_obj@meta.data)
  for (i in seq_along(meta_cols)) cat(sprintf("%d. %s\n", i, meta_cols[i]))
  
  # 选择聚类列
  cat("\n请选择聚类列编号: ")
  col_idx <- as.integer(readline())
  cluster_column <- meta_cols[col_idx]
  
  # 显示可用的降维
  cat("\n可用的降维方法:\n")
  reductions <- names(seurat_obj@reductions)
  for (i in seq_along(reductions)) cat(sprintf("%d. %s\n", i, reductions[i]))
  
  # 选择降维方法
  cat("\n请选择降维方法编号: ")
  red_idx <- as.integer(readline())
  reduction <- reductions[red_idx]
  
  # 显示簇信息
  cat(sprintf("\n'%s'列中的簇:\n", cluster_column))
  clusters <- unique(seurat_obj@meta.data[[cluster_column]])
  for (i in seq_along(clusters)) cat(sprintf("%d. %s\n", i, clusters[i]))
  
  # 选择目标簇
  cat("\n请输入目标簇编号(多个用逗号分隔，留空处理所有簇): ")
  cluster_input <- readline()
  target_clusters <- if (cluster_input == "") NULL else {
    indices <- as.integer(unlist(strsplit(cluster_input, ",")))
    clusters[indices]
  }
  
  # 设置密度阈值
  cat("\n请输入密度阈值(0-1之间): ")
  density_threshold <- as.numeric(readline())
  
  # 选择过滤模式
  cat("\n请选择过滤模式:\n")
  cat("1. 仅去除散在细胞 (scattered)\n")
  cat("2. 仅去除渗入细胞 (infiltrating)\n")
  cat("3. 同时去除散在和渗入细胞 (both)\n")
  mode_idx <- as.integer(readline())
  filter_mode <- c("scattered", "infiltrating", "both")[mode_idx]
  
  # 设置输出目录
  cat("\n请输入结果输出目录: ")
  output_dir <- readline()
  
  # 运行优化
  result <- optimize_clusters(
    seurat_obj = seurat_obj,
    target_clusters = target_clusters,
    cluster_column = cluster_column,
    reduction = reduction,
    density_threshold = density_threshold,
    filter_mode = filter_mode,
    output_dir = output_dir
  )
  
  # 显示处理结果摘要
  cat("\n========== 处理完成! ==========\n")
  cat(sprintf("原始细胞数: %d\n", ncol(result$original)))
  cat(sprintf("过滤后细胞数: %d (保留%.1f%%)\n", 
              ncol(result$filtered), 
              100*ncol(result$filtered)/ncol(result$original)))
  cat(sprintf("移除的散在细胞数: %d\n", length(result$scattered_cells)))
  cat(sprintf("移除的渗入细胞数: %d\n", length(result$infiltrating_cells)))
  cat(sprintf("结果已保存至: %s\n", output_dir))
  
  return(result)
}

# 使用示例
# 方式1: 交互式运行
# result <- run_interactive()

# 方式2: 直接调用函数
seurat_obj <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_ident_11m26d.rds")
levels(seurat_obj@meta.data$celltype)
result <- optimize_clusters(
  seurat_obj = seurat_obj,
  target_clusters = c(),  # 可选，留空处理所有簇
  cluster_column = "celltype",           # 聚类列名
  reduction = "tsne",                           # 降维方法: umap或tsne
  density_threshold = 0.2,                      # 密度阈值
  filter_mode = "both",                         # 过滤模式: scattered/infiltrating/both
  output_dir = "./cluster_optimization"         # 输出目录
)
