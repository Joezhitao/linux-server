#############################################################
# 单细胞比例冲击图 (Alluvial Plot) 参数化脚本 - 简化版
# 功能：读取Seurat对象，计算细胞比例，支持亚群筛选，自定义颜色
# 特点：去除强制排序参数，默认使用Seurat对象中的因子顺序
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载包与环境设置 ------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggalluvial)
  library(tidyverse)
  library(dplyr)
  library(readr)
})

# 2. 核心绘图函数 ----------------------------------------------------------

run_alluvial_plot <- function(
    rds_path,                    # Seurat对象的路径
    output_dir = "./alluvial_results", # 输出目录
    file_prefix = "Alluvial",    # 输出文件前缀
    
    # 核心列设置
    group_col = "orig.ident",    # x轴的分组列名 (如样本名)
    celltype_col = "cellcluster", # 细胞类型列名 (用于绘图和统计的最终列)
    
    # 颜色设置
    custom_colors = NULL,        # 命名向量 c("TypeA"="#xx", ...)，NULL则自动生成
    
    # 亚群筛选设置
    subset_mode = FALSE,         # 是否开启亚群筛选模式
    subset_col = NULL,           # 筛选所用的列
    subset_idents = NULL         # 保留哪些ID
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat(sprintf("开始冲击图分析: %s\n", file_prefix))
  
  # --- 1. 读取数据 ---
  cat(sprintf("  [INFO] 读取数据: %s\n", rds_path))
  if (!file.exists(rds_path)) stop("RDS文件不存在！")
  pbmc <- readRDS(rds_path)
  
  # --- 2. 亚群筛选 (可选) ---
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("亚群模式下需指定 subset_col 和 subset_idents")
    cat(sprintf("  [INFO] 亚群模式: 筛选 %s 属于 {%s}\n", subset_col, paste(subset_idents, collapse=",")))
    Idents(pbmc) <- subset_col
    pbmc <- subset(pbmc, idents = subset_idents)
  }
  
  # --- 3. 提取元数据 ---
  meta <- pbmc@meta.data
  
  # 检查列是否存在
  if (!celltype_col %in% colnames(meta)) stop(paste("列名不存在:", celltype_col))
  if (!group_col %in% colnames(meta)) stop(paste("列名不存在:", group_col))
  
  # --- 4. 计算比例 ---
  cat("  [INFO] 计算细胞比例...\n")
  # 统计数量
  cell_count <- table(meta[[celltype_col]], meta[[group_col]])
  
  # 计算比例
  cell_ratio <- prop.table(cell_count, margin = 2) %>%
    as.data.frame() %>%
    set_names(c("celltype", "group", "ratio"))
  
  # 导出CSV
  csv_path <- file.path(output_dir, paste0(file_prefix, "_ratio.csv"))
  write_csv(cell_ratio, csv_path)
  cat(sprintf("  [INFO] 比例表已保存: %s\n", csv_path))
  
  # --- 5. 颜色处理 ---
  # 如果提供了自定义颜色，我们需要确保绘图数据的因子顺序包含了这些颜色的Key
  if (!is.null(custom_colors)) {
    # 尝试将细胞类型设置为因子，顺序优先匹配颜色的顺序
    # 这样图例的顺序就会和颜色定义的顺序一致
    color_levels <- names(custom_colors)
    existing_levels <- as.character(unique(cell_ratio$celltype))
    
    # 取交集，防止颜色定义里有多余的，或者数据里有没定义的
    final_levels <- intersect(color_levels, existing_levels)
    # 把数据里有但颜色里没定义的补在后面
    remaining_levels <- setdiff(existing_levels, color_levels)
    final_levels <- c(final_levels, remaining_levels)
    
    cell_ratio$celltype <- factor(cell_ratio$celltype, levels = final_levels)
  }
  
  # --- 6. 绘图 (Alluvial Plot) ---
  cat("  [INFO] 正在绘图...\n")
  
  pp <- ggplot(cell_ratio, aes(x = group, y = ratio, fill = celltype,
                               stratum = celltype, alluvium = celltype)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) + 
    theme_classic() +
    labs(x = 'Group', y = 'Cell Ratio', fill = 'Cell Type', title = file_prefix) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    )
  
  # 应用颜色
  if (!is.null(custom_colors)) {
    pp <- pp + scale_fill_manual(values = custom_colors)
  }
  
  # 添加冲击图图层
  p_final <- pp +
    geom_col(width = 0.6, color = NA, linewidth = 0.5) +
    geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0, color = 'white', linewidth = 0.5) +
    geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0, fill = NA, color = 'white', linewidth = 0.5)
  
  # --- 7. 保存图片 ---
  pdf_path <- file.path(output_dir, paste0(file_prefix, "_Plot.pdf"))
  pdf(pdf_path, width = 8, height = 6)
  print(p_final)
  dev.off()
  
  cat(sprintf("  [SUCCESS] 绘图完成: %s\n", pdf_path))
  cat("=================================================\n")
}


# 3. 参数配置与运行示例 ----------------------------------------------------

# --- 常用参数准备 ---
rds_file <- "/home/lin/c_group/hep.rds"

# 自定义颜色 (可选)
# 只要这里的名字和 celltype_col 列里的内容对应即可
my_colors <- c(
  "Regenerative-Unit-Hep1" = "#D62728FF", 
  "Regenerative-Unit-Hep2" = "#9467BDFF",
  "Regenerative-Unit-Hep3" = "#8C564BFF",
  "Regenerative-Unit-Hep4" = "#E377C2FF"
)

# ==============================================
# 场景 A: 全局分析 (使用现有的 cellcluster 列)
# ==============================================
run_alluvial_plot(
  rds_path = rds_file,
  output_dir = "/home/lin/c_group/alluvial_results",
  file_prefix = "Global_Analysis",
  
  group_col = "orig.ident",   # X轴用什么分组
  celltype_col = "cellcluster", # 细胞类型用哪一列
  
  custom_colors = my_colors,  # 使用自定义颜色
  
  subset_mode = FALSE
)

# ==============================================
# 场景 B: 亚群分析 (筛选特定样本)
# ==============================================
run_alluvial_plot(
  rds_path = rds_file,
  output_dir = "/home/lin/c_group/alluvial_results",
  file_prefix = "Subgroup_Regenerative",
  
  group_col = "orig.ident",
  celltype_col = "seurat_clusters",
  
  custom_colors = my_colors,       # 演示：留空使用默认颜色
  
  subset_mode = TRUE,
  subset_col = "cellcluster",
  subset_idents = c("Regenerative-Unit-Hep")
)
color
seu_obj <- readRDS("/home/lin/c_group/hep.rds")
levels(seu_obj@meta.data$cellcluster)
