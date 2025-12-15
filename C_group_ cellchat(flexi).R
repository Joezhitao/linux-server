#############################################################
# 单细胞 CellChat 标准化分析脚本 (基于官方文档重构版)
# 功能：CellChat对象构建、通讯推断、可视化 (Circle/Bubble)
# 适配：Human, Mouse, Rat (Rat自动映射到Mouse库)
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(future)
})

# 设置并行处理 (官方推荐)
plan("multisession", workers = 20)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100 * 1024^3) # 增加内存限制

# 2. 核心分析函数定义 ------------------------------------------------------

run_cellchat_workflow <- function(
    seurat_obj, 
    group_by, 
    species = "mouse", 
    output_dir = "./cellchat_output",
    file_prefix = "Analysis",
    subset_mode = FALSE,         
    subset_col = NULL,           
    subset_idents = NULL
) {
  
  # --- 0. 初始化目录 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(paste0("\n========== 开始 CellChat 分析: ", file_prefix, " ==========\n"))
  
  # --- 1. 数据准备与亚群筛选 ---
  current_obj <- seurat_obj
  
  # 如果开启亚群模式
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("错误: 亚群模式下必须指定 subset_col 和 subset_idents")
    cat(sprintf("  [MODE] 正在提取亚群: %s in %s\n", paste(subset_idents, collapse=","), subset_col))
    Idents(current_obj) <- current_obj@meta.data[[subset_col]]
    current_obj <- subset(current_obj, idents = subset_idents)
  }
  
  # 设置后续分析用的分组 (例如 seurat_clusters)
  if (!group_by %in% colnames(current_obj@meta.data)) stop(paste("未找到分组列:", group_by))
  Idents(current_obj) <- current_obj@meta.data[[group_by]]
  
  # 提取数据 (兼容 Seurat V3/V4/V5)
  cat("  [DATA] 提取表达矩阵...\n")
  data.input <- tryCatch({
    GetAssayData(current_obj, assay = "RNA", layer = "data")
  }, error = function(e) {
    GetAssayData(current_obj, assay = "RNA", slot = "data")
  })
  meta.data <- current_obj@meta.data
  
  # --- 2. 创建 CellChat 对象 ---
  cat("  [INIT] 创建 CellChat 对象...\n")
  cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = group_by)
  
  # --- 3. 设置数据库 (Rat -> Mouse, Human -> Human) ---
  cat(sprintf("  [DB] 设置物种数据库: %s\n", species))
  species_lower <- tolower(species)
  
  if (species_lower == "human") {
    CellChatDB <- CellChatDB.human
  } else if (species_lower %in% c("mouse", "rat")) {
    if (species_lower == "rat") cat("    -> 检测到 Rat，映射至 Mouse 数据库\n")
    CellChatDB <- CellChatDB.mouse
  } else {
    stop("仅支持 human, mouse, rat")
  }
  
  # 使用 "Secreted Signaling" (官方推荐)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
  cellchat@DB <- CellChatDB.use
  
  # --- 4. 预处理 (识别高变配受体) ---
  cat("  [PREP] 预处理与高变基因识别...\n")
  cellchat <- subsetData(cellchat) # 即使全基因也推荐运行
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # --- 5. 推断通讯网络 (核心步骤) ---
  cat("  [RUN] 推断通讯概率 (type='triMean')...\n")
  # type = "triMean" producing fewer but stronger interactions (Robust)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  
  # 过滤少于10个细胞的通讯
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # 推断通路水平通讯
  cat("  [RUN] 推断通路水平通讯...\n")
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 聚合网络
  cellchat <- aggregateNet(cellchat)
  
  # --- 6. 保存对象 ---
  rds_path <- file.path(output_dir, paste0(file_prefix, "_cellchat.rds"))
  saveRDS(cellchat, rds_path)
  cat(sprintf("  [SAVE] CellChat对象已保存: %s\n", rds_path))
  
  # --- 7. 自动绘图 (图1 & 图2) ---
  plot_standard_charts(cellchat, output_dir, file_prefix)
  
  return(cellchat)
}

# 3. 标准绘图函数 (图1和图2) -----------------------------------------------

plot_standard_charts <- function(cellchat, output_dir, file_prefix) {
  cat("  [PLOT]正在生成标准图表...\n")
  groupSize <- as.numeric(table(cellchat@idents))
  
  # --- 图1: Circle Plot (互作数量与强度) ---
  pdf(file.path(output_dir, paste0(file_prefix, "_1_Circle_Net.pdf")), width = 10, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Interaction weights")
  dev.off()
  
  # --- 图2: Bubble Plot (全局所有显著通讯) ---
  # 注意：如果通讯过多，气泡图可能很大
  tryCatch({
    p <- netVisual_bubble(cellchat, remove.isolate = FALSE) + 
      ggtitle("All Significant Interactions")
    
    ggsave(file.path(output_dir, paste0(file_prefix, "_2_All_Bubble.pdf")), 
           plot = p, width = 12, height = 18, limitsize = FALSE)
  }, error = function(e) {
    cat("    [WARN] 全局气泡图绘制失败 (可能无显著结果): ", e$message, "\n")
  })
}

# 4. 自定义配受体查询绘图函数 (图3 - 升级版: 支持筛选细胞群) -----------------

plot_custom_bubble <- function(
    cellchat_obj,
    ligand = NULL,          # 筛选配体基因名 (例如 "Tgfb1")
    receptor = NULL,        # 筛选受体基因名 (例如 "Tgfbr1")
    source_groups = NULL,   # 筛选发送细胞群名称 (例如 c("Fibroblast", "Macrophage"))
    target_groups = NULL,   # 筛选接收细胞群名称 (例如 c("Epithelial"))
    output_dir,
    file_prefix = "Custom"
) {
  if (is.null(cellchat_obj)) stop("请提供 CellChat 对象！")
  
  cat("\n========== 开始自定义查询绘图 (Advanced) ==========\n")
  
  # 1. 检索所有推断出的互作
  # subsetCommunication 返回一个数据框，包含 source, target, ligand, receptor 等列
  df.net <- subsetCommunication(cellchat_obj)
  
  if (is.null(df.net) || nrow(df.net) == 0) {
    stop("该对象中没有显著的通讯结果。")
  }
  
  # --- 筛选逻辑开始 ---
  
  # 2.1 筛选配体和受体 (基因层面)
  gene_mask <- rep(TRUE, nrow(df.net))
  query_desc_gene <- ""
  
  if (!is.null(ligand)) {
    gene_mask <- gene_mask & (df.net$ligand == ligand)
    query_desc_gene <- paste0(query_desc_gene, "_L-", ligand)
  }
  if (!is.null(receptor)) {
    gene_mask <- gene_mask & (df.net$receptor == receptor)
    query_desc_gene <- paste0(query_desc_gene, "_R-", receptor)
  }
  
  if (is.null(ligand) && is.null(receptor)) {
    cat("  [INFO] 未指定配体或受体，将展示指定细胞群间的所有互作。\n")
    query_desc_gene <- "_AllGenes"
  }
  
  # 2.2 筛选发送者和接收者 (细胞群层面)
  cell_mask <- rep(TRUE, nrow(df.net))
  query_desc_cell <- ""
  
  if (!is.null(source_groups)) {
    # 检查输入的细胞群是否存在
    valid_sources <- intersect(source_groups, unique(df.net$source))
    if (length(valid_sources) == 0) {
      cat(sprintf("  [WARN] 输入的 source_groups 不存在于数据中。现有组: %s\n", paste(unique(df.net$source), collapse=",")))
      return(NULL)
    }
    cell_mask <- cell_mask & (df.net$source %in% valid_sources)
    query_desc_cell <- paste0(query_desc_cell, "_Src-", paste(valid_sources, collapse = "+"))
  }
  
  if (!is.null(target_groups)) {
    valid_targets <- intersect(target_groups, unique(df.net$target))
    if (length(valid_targets) == 0) {
      cat(sprintf("  [WARN] 输入的 target_groups 不存在于数据中。现有组: %s\n", paste(unique(df.net$target), collapse=",")))
      return(NULL)
    }
    cell_mask <- cell_mask & (df.net$target %in% valid_targets)
    query_desc_cell <- paste0(query_desc_cell, "_Tgt-", paste(valid_targets, collapse = "+"))
  }
  
  # 3. 应用筛选
  final_mask <- gene_mask & cell_mask
  df.subset <- df.net[final_mask, ]
  
  if (nrow(df.subset) == 0) {
    cat("  [WARN] 没有找到符合所有条件的互作结果。\n")
    return(NULL)
  }
  
  # 4. 绘图参数准备
  # netVisual_bubble 需要两个关键参数:
  # sources.use: 发送者索引/名称
  # targets.use: 接收者索引/名称
  # pairLR.use:  配受体对数据框
  
  # 如果用户指定了 group，直接用用户的；如果没指定，则使用筛选结果中存在的 group
  # 注意：netVisual_bubble 接受的是向量，用于限定气泡图的轴
  real_sources <- if(!is.null(source_groups)) source_groups else unique(df.subset$source)
  real_targets <- if(!is.null(target_groups)) target_groups else unique(df.subset$target)
  
  pairLR.use <- unique(data.frame(interaction_name = df.subset$interaction_name))
  
  cat(sprintf("  [PLOT] 筛选后剩余 %d 条互作关系，正在绘制...\n", nrow(df.subset)))
  
  # 5. 绘图
  # 构建完整的文件名描述，避免过长可截断，这里简化处理
  full_desc <- paste0(query_desc_gene, query_desc_cell)
  # 移除可能的开头下划线
  if(substr(full_desc, 1, 1) == "_") full_desc <- substr(full_desc, 2, nchar(full_desc)) 
  
  # 限制文件名长度，防止报错
  if(nchar(full_desc) > 50) full_desc <- paste0(substr(full_desc, 1, 45), "_etc")
  
  tryCatch({
    p <- netVisual_bubble(
      cellchat_obj, 
      sources.use = real_sources, 
      targets.use = real_targets, 
      pairLR.use = pairLR.use, 
      remove.isolate = FALSE
    ) + ggtitle(paste0("Custom Query"))
    
    out_file <- file.path(output_dir, paste0(file_prefix, "_3_Bubble_", full_desc, ".pdf"))
    
    # 动态调整高度：互作对越多，图片越高
    plot_height <- max(4, nrow(pairLR.use) * 0.4 + 2) 
    
    ggsave(out_file, plot = p, width = 8, height = plot_height, limitsize = FALSE)
    cat(sprintf("  [SUCCESS] 图片已保存: %s\n", out_file))
    
  }, error = function(e) {
    cat("  [ERROR] 绘图失败: ", e$message, "\n")
  })
}

# ==========================================================================
# 5. 参数设置与运行区 (用户只需修改这里) -----------------------------------
# ==========================================================================

# --- 基础参数 ---
RDS_FILE    <- "/home/lin/c_group/T_HSCs_combined.rds"
SPECIES     <- "rat"              # 选项: "human", "mouse", "rat"
GROUP_BY    <- "seurat_clusters"  # 分组列名
OUTPUT_ROOT <- "/home/lin/c_group/cellchat" # 输出目录

# 读取数据 (只读一次)
if (!exists("seurat_obj")) {
  cat("正在读取 Seurat RDS 文件...\n")
  seurat_obj <- readRDS(RDS_FILE)
  # 确保分组是因子 (Factor)
  if(GROUP_BY %in% colnames(seurat_obj@meta.data)){
    seurat_obj@meta.data[[GROUP_BY]] <- as.factor(seurat_obj@meta.data[[GROUP_BY]])
  }
}

# --------------------------------------------------------------------------
# 模式选择：请根据需要取消注释并运行某一块
# --------------------------------------------------------------------------

# 【场景 A: 全局分析】
# 运行后会生成对象 res_global 和两张基础图
res_global <- run_cellchat_workflow(
  seurat_obj = seurat_obj,
  group_by = GROUP_BY,
  species = SPECIES,
  output_dir = file.path(OUTPUT_ROOT, "Global"),
  file_prefix = "Global",
  subset_mode = FALSE
)


# 【场景 B: 亚群分析】 (示例：仅分析 orig.ident 为 "Tumor" 的细胞)
# res_sub <- run_cellchat_workflow(
#   seurat_obj = seurat_obj,
#   group_by = GROUP_BY,
#   species = SPECIES,
#   output_dir = file.path(OUTPUT_ROOT, "Subgroup"),
#   file_prefix = "Tumor_Only",
#   subset_mode = TRUE,
#   subset_col = "orig.ident",   # 筛选依据的列
#   subset_idents = c("Tumor")   # 保留的组名
# )


# --------------------------------------------------------------------------
# 【图3: 自定义查询绘图】 (独立运行)
# --------------------------------------------------------------------------
# 不需要重新跑分析，直接用上面生成的 res_global 或 res_sub 对象

# 例子 1: 查询所有以 Met 为受体的通讯 (Rat/Mouse 基因通常首字母大写)
plot_custom_bubble(
  cellchat_obj = res_global,
  ligand = "Tgfb1",         # 筛选基因
  source_groups = c("CD4+ Regulatory T cells", "CD8+ Effector Memory T cells"),   # 筛选发送细胞 (注意必须是字符串，如果是数字ID需加引号)
  target_groups = c("Fibrogenic HSC"),   # 筛选接收细胞
  output_dir = file.path(OUTPUT_ROOT, "Global"),
  file_prefix = "Tgfb1_C1_to_C2"
)
levels(as.factor(seurat_obj@meta.data$seurat_clusters))

# 例子 2: 查询特定配体 Hgf
# plot_custom_bubble(
#   cellchat_obj = res_global,
#   ligand = "Hgf",
#   receptor = NULL,
#   output_dir = file.path(OUTPUT_ROOT, "Global"),
#   file_prefix = "Query_Hgf"
# )
