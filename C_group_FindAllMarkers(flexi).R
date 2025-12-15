#############################################################
# 单细胞差异基因分析(Markers)参数化脚本 (v2: 修复Seurat V5报错版)
# 修复内容：调整 JoinLayers 顺序，解决 Assay5 subset 后结构失效问题
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 核心功能函数定义 ------------------------------------------------------

run_marker_analysis <- function(
    seurat_obj,                  
    group_by = "seurat_clusters", 
    subset_mode = FALSE,          
    subset_col = NULL,            
    subset_idents = NULL,         
    logfc_threshold = 0.25,       
    min_pct = 0.1,                
    output_dir = "./marker_results", 
    file_prefix = "Markers_Analysis",
    plot_width = 12,
    plot_height = 6
) {
  
  # --- 0. 初始化与环境设置 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("=================================================\n")
  cat(sprintf("开始差异基因(Markers)分析流程 | %s\n", file_prefix))
  
  final_obj <- seurat_obj
  
  # --- [关键修改] 1. 先进行 JoinLayers ---
  # 在 subset 之前合并图层，防止 subset 后出现空图层导致 Assay5 结构报错
  # 只有当对象是 V5 Assay 且有多层时，这一步才关键
  cat("  [INFO] 正在执行 JoinLayers (预处理)...\n")
  tryCatch({
    final_obj <- JoinLayers(final_obj)
  }, error = function(e) {
    cat("  [WARN] JoinLayers 遇到问题或无需 Join，尝试继续...\n")
  })
  
  # --- 2. 亚群筛选逻辑 ---
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("错误: 开启亚群模式时，必须设置 subset_col 和 subset_idents！")
    
    cat(sprintf("  [INFO] 模式: 亚群分析 (筛选 %s == %s)\n", subset_col, paste(subset_idents, collapse=",")))
    
    # 检查列是否存在
    if (!subset_col %in% colnames(final_obj@meta.data)) stop(paste("错误: 亚群筛选列不存在:", subset_col))
    
    Idents(final_obj) <- subset_col
    final_obj <- subset(final_obj, idents = subset_idents)
    
    cat(sprintf("  [INFO] 亚群筛选后剩余细胞数: %d\n", ncol(final_obj)))
    
    # 筛选后如果细胞数为0，直接报错停止
    if (ncol(final_obj) == 0) stop("错误: 筛选后没有剩余细胞，请检查 subset_idents 是否正确！")
  } else {
    cat("  [INFO] 模式: 全局分析 (使用所有细胞)\n")
  }
  
  # --- 3. 设置分组 ---
  if (!group_by %in% colnames(final_obj@meta.data)) stop(paste("错误: 分组列不存在:", group_by))
  
  Idents(final_obj) <- group_by
  # 移除未使用的 level，防止绘图或分析时出现空组
  final_obj@meta.data[[group_by]] <- droplevels(as.factor(final_obj@meta.data[[group_by]]))
  Idents(final_obj) <- final_obj@meta.data[[group_by]]
  
  cat(sprintf("  [INFO] 当前分析分组依据: %s (共 %d 个组)\n", group_by, length(unique(Idents(final_obj)))))
  
  # --- 4. 运行 FindAllMarkers ---
  cat(sprintf("  [INFO] 正在运行 FindAllMarkers (logfc=%.2f, min.pct=%.2f)...\n", logfc_threshold, min_pct))
  
  # 再次确保 DefaultAssay 正确，通常是 RNA
  # 如果你有 integrated assay，这里可能需要根据情况调整，一般找 marker 建议用 RNA count
  DefaultAssay(final_obj) <- "RNA" 
  
  combined_markers <- FindAllMarkers(
    object = final_obj, 
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = FALSE
  )
  
  # 检查是否找到了 marker
  if (nrow(combined_markers) == 0) {
    warning("  [WARN] 未找到符合阈值的 Marker！将跳过导出和绘图。")
    return(NULL)
  }
  
  all_markers <- combined_markers %>% group_by(cluster)
  
  # --- 5. 导出 CSV 结果 ---
  out_csv <- file.path(output_dir, paste0(file_prefix, "_All_Markers.csv"))
  cat(sprintf("  [INFO] 正在导出 Marker 表: %s\n", basename(out_csv)))
  
  write.csv(all_markers, 
            file = out_csv, 
            quote = FALSE, 
            row.names = FALSE)
  
  # --- 6. 绘制并保存 UMAP 图 ---
  cat("  [INFO] 正在绘制 UMAP...\n")
  
  plot_title <- paste0("UMAP by ", group_by, " (", file_prefix, ")")
  
  # 检查是否有 umap reduction
  if (!"umap" %in% names(final_obj@reductions)) {
    cat("  [WARN] 对象中没有 'umap' 降维信息，尝试使用 'tsne' 或跳过绘图...\n")
    red <- ifelse("tsne" %in% names(final_obj@reductions), "tsne", NULL)
  } else {
    red <- "umap"
  }
  
  if (!is.null(red)) {
    p1 <- DimPlot(final_obj, reduction = red, label = TRUE, group.by = group_by, pt.size = 0.4) + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5))
    
    out_pdf <- file.path(output_dir, paste0(file_prefix, "_UMAP.pdf"))
    ggsave(out_pdf, p1, width = plot_width, height = plot_height, dpi = 300)
    cat(sprintf("  [INFO] UMAP 图已保存: %s\n", basename(out_pdf)))
  } else {
    cat("  [WARN] 无法绘图：未找到 UMAP 或 tSNE 坐标。\n")
  }
  
  # --- 7. 结束 ---
  cat("=================================================\n")
  cat(sprintf("分析完成！结果保存在: %s\n", output_dir))
  
  return(all_markers)
}


# 3. 使用示例 ---------------------------------------------------------------

# 读取数据
seurat_obj <- readRDS("/home/lin/c_group/hep.rds")

# ==============================================
# 场景 A: 全局分析 (原有代码逻辑)
# ==============================================
# 对所有细胞，按 seurat_clusters 分组找差异基因
res_global <- run_marker_analysis(
  seurat_obj = seurat_obj,
  group_by = "seurat_clusters",       # 按照 cluster 分组
  subset_mode = FALSE,                # 关闭亚群筛选
  logfc_threshold = 0.25,
  output_dir = "/home/lin/c_group/marker_results/global",
  file_prefix = "Global_Clusters"
)

# ==============================================
# 场景 B: 亚群分析 (示例)
# ==============================================
# 例如：只取 seurat_clusters 为 0, 1, 2 的细胞，
# 然后在这些细胞中，按 "sample_type" (假设有这一列) 找差异基因
# 或者：只取 "CellType" 为 "T_cell" 的细胞，按 "stimulation" 分组找差异

res_subset <- run_marker_analysis(
  seurat_obj = seurat_obj,
  group_by = "seurat_clusters",       # 在筛选后的细胞中，依然按 cluster 找差异(或者换成其他列)
  subset_mode = TRUE,                 # 开启亚群筛选
  subset_col = "cellcluster",     # 筛选依据列
  subset_idents = c("Core-Metabolic-Hep"),         # 只保留 0, 1, 2 簇
  logfc_threshold = 0.25,
  output_dir = "/home/lin/c_group/marker_results/subset_012",
  file_prefix = "Subset_Clusters_Core"
)
levels(seurat_obj@meta.data$cellcluster)
