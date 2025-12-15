#############################################################
# Seurat 数据灵活合并工具 (v2.0)
# 功能：
#   1. 支持“全集+全集”、“子集+全集”、“子集+子集”任意组合
#   2. 自动处理 Seurat V5 的 Layer 分裂问题 (JoinLayers)
#   3. 自动添加新标签列，支持自动标准化
#   4. 兼容传入 RDS 文件路径或内存中的 Seurat 对象
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载包与环境设置 ------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# 设置内存限制，防止合并大文件时报错
options(future.globals.maxSize = 100 * 1024^3)

# 2. 定义通用合并函数 ------------------------------------------------------

merge_seurat_custom <- function(
    # --- 数据源 1 参数 ---
  obj1,                       # 输入：可以是 Seurat对象，也可以是 .rds 文件路径
  subset_col1 = NULL,         # 筛选列名 (如 "seurat_clusters")
  subset_idents1 = NULL,      # 筛选标签 (如 "Treg")，设为 NULL 则保留该数据所有细胞
  label1 = "Group1",          # 合并后的新标签名 (如 "CD4+ Treg")
  
  # --- 数据源 2 参数 ---
  obj2,                       # 输入：同上
  subset_col2 = NULL,         # 筛选列名
  subset_idents2 = NULL,      # 筛选标签，设为 NULL 则保留该数据所有细胞
  label2 = "Group2",          # 合并后的新标签名 (如 "CD8+ T")
  
  # --- 通用参数 ---
  output_col = "celltype",    # 结果元数据中存放新标签的列名
  project_name = "Combined",  # Seurat Project Name
  do_normalize = TRUE,        # 是否在合并后立即执行 NormalizeData (推荐)
  output_file = NULL          # 保存路径，NULL 则不保存
) {
  
  cat("\n=================================================\n")
  cat("启动 Seurat 数据合并流程\n")
  
  # --- 内部函数：读取与筛选单个对象 ---
  process_single_obj <- function(input_obj, col, idents, label_name, source_name) {
    # 1. 智能读取
    if (is.character(input_obj)) {
      if (!file.exists(input_obj)) stop(sprintf("错误: 文件不存在 -> %s", input_obj))
      cat(sprintf("  [LOAD] %s: 从文件读取 %s\n", source_name, basename(input_obj)))
      obj <- readRDS(input_obj)
    } else {
      cat(sprintf("  [LOAD] %s: 使用内存中的对象\n", source_name))
      obj <- input_obj
    }
    
    # 2. 智能筛选
    if (!is.null(idents)) {
      # 检查列名
      if (is.null(col)) stop(sprintf("错误: %s 设置了 subset_idents 但未设置 subset_col", source_name))
      if (!col %in% colnames(obj@meta.data)) stop(sprintf("错误: 列名 '%s' 在 %s 中不存在", col, source_name))
      
      cat(sprintf("  [SUBSET] %s: 按 '%s' 筛选 -> %s\n", source_name, col, paste(idents, collapse=",")))
      
      # 执行 subset
      Idents(obj) <- col
      obj <- subset(obj, idents = idents)
    } else {
      cat(sprintf("  [KEEP] %s: 未设置筛选，保留全部细胞\n", source_name))
    }
    
    # 3. 检查细胞并打标签
    if (ncol(obj) == 0) stop(sprintf("错误: %s 筛选后剩余 0 个细胞，请检查标签名是否正确！", source_name))
    cat(sprintf("     -> 保留细胞数: %d\n", ncol(obj)))
    
    obj[[output_col]] <- label_name
    return(obj)
  }
  
  # --- 处理两个对象 ---
  s1 <- process_single_obj(obj1, subset_col1, subset_idents1, label1, "Obj1")
  s2 <- process_single_obj(obj2, subset_col2, subset_idents2, label2, "Obj2")
  
  # --- 合并 (Merge) ---
  cat("  [MERGE] 正在合并数据集...\n")
  # add.cell.ids 防止条形码冲突
  combined <- merge(x = s1, y = s2, add.cell.ids = c(label1, label2), project = project_name)
  
  # --- Seurat V5 修复 (关键步骤) ---
  # V5 merge 后数据在不同 layer，必须合并，否则 CellChat/FindMarkers 会报错
  cat("  [FIX] 执行 JoinLayers (Seurat V5 兼容性修复)...\n")
  combined <- JoinLayers(combined)
  
  # --- 标准化 (可选) ---
  if (do_normalize) {
    cat("  [PROCESS] 执行 NormalizeData (LogNormalize)...\n")
    combined <- NormalizeData(combined, verbose = FALSE)
  }
  
  # --- 整理 Idents ---
  Idents(combined) <- combined@meta.data[[output_col]]
  
  # --- 输出结果 ---
  cat("-------------------------------------------------\n")
  cat("合并完成！新分组统计:\n")
  print(table(combined@meta.data[[output_col]]))
  
  # --- 保存 ---
  if (!is.null(output_file)) {
    saveRDS(combined, output_file)
    cat(sprintf("  [SAVE] 结果已保存至: %s\n", output_file))
  }
  
  cat("=================================================\n")
  return(combined)
}


# 3. 实际应用场景示例 ------------------------------------------------------

# 假设文件路径如下 (你可以直接传对象变量，也可以传路径字符串)
data1 <- "/home/lin/c_group/CD4+T_ident.rds"
data2 <- "/home/lin/c_group/HSCs_ident.rds"

# 为了方便演示，先读入内存 (也可以不读，直接传路径)
data1 <- readRDS(data1)
data2 <- readRDS(data2)

# 查看一下现有的簇，方便确认标签名
# levels(CD4T@meta.data$seurat_clusters)
# levels(CD8T@meta.data$seurat_clusters)


# ==========================================================================
# 场景 A: 你的核心需求 (CD4 只取 Treg + CD8 全集)
# ==========================================================================
cat("\n--- 运行场景 A ---\n")

seu_treg_cd8 <- merge_seurat_custom(
  # --- 对象1: CD4 (只取子集) ---
  obj1 = data1, 
  subset_col1 = "seurat_clusters",            # 依据这列筛选
  subset_idents1 = "CD4+ Regulatory T cells", # 只保留这一簇
  label1 = "CD4+ Treg",                       # 新名字
  
  # --- 对象2: CD8 (全都要) ---
  obj2 = data2,
  subset_col2 = NULL,     # 不筛选
  subset_idents2 = NULL,  # NULL = 保留所有
  label2 = "HSCs",      # 新名字
  
  # --- 结果设置 ---
  output_col = "celltype",
  output_file = "/home/lin/c_group/Treg_HSCs_combined.rds"
)


# ==========================================================================
# 场景 B: 两个数据都只取特定子集 (例如：只比较两个系的 Naive 细胞)
# ==========================================================================
cat("\n--- 运行场景 B ---\n")

# 假设 CD8 数据里也有 Naive 标签
seu_naive_compare <- merge_seurat_custom(
  # --- CD4 Naive ---
  obj1 = CD4T,
  subset_col1 = "seurat_clusters",
  subset_idents1 = "CD4+ Naive/Central Memory T cells",
  label1 = "CD4_Naive",
  
  # --- CD8 Naive (假设存在这个簇，演示用) ---
  # 注意：如果不知道 CD8 具体标签，这里会报错提醒
  obj2 = CD8T,
  subset_col2 = "seurat_clusters", 
  subset_idents2 = c("1", "2"), # 也可以写 Cluster ID
  label2 = "CD8_Naive",
  
  output_col = "celltype_naive",
  output_file = "/home/lin/c_group/Naive_Compare.rds"
)


# ==========================================================================
# 场景 C: 两个全集直接暴力合并 (All + All)
# ==========================================================================
cat("\n--- 运行场景 C ---\n")

seu_all_combined <- merge_seurat_custom(
  obj1 = path_cd4,   # 这里演示直接传路径，不占用额外内存
  label1 = "Total_CD4",
  
  obj2 = path_cd8,   # 这里演示直接传路径
  label2 = "Total_CD8",
  
  output_col = "major_lineage",
  output_file = "/home/lin/c_group/All_T_cells.rds"
)
