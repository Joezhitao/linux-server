#############################################################
# 单细胞 irGSEA 参数化分析脚本 (v4.2: 增加亚群分析支持)
# 功能：多方法评分、差异通路识别、RRA聚合、可视化
# 更新：新增 subset_mode 参数，支持仅分析特定亚群 (如仅分析 Con 组)
# 日期：2023-12-07
#############################################################

rm(list = ls())
gc()

# 1. 加载包与环境设置 ------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(irGSEA)
  library(tidyverse)
  library(ggplot2)
  library(cowplot) # 辅助绘图
})

# 设置大内存支持 (irGSEA计算量大)
options(future.globals.maxSize = 1000 * 1024^3)

# 2. 辅助功能函数 ----------------------------------------------------------

# --- 函数 A: 稳健的 GMT 读取 ---
read_gmt_custom <- function(file_path) {
  if (!file.exists(file_path)) stop("GMT file not found: ", file_path)
  lines <- readLines(file_path)
  gene_sets <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) > 2) {
      gene_sets[[fields[1]]] <- fields[3:length(fields)]
    }
  }
  return(gene_sets)
}

# --- 函数 B: 基因名大写转换 (修复版) ---
convert_to_human_symbol <- function(seurat_obj, species) {
  species <- tolower(species)
  if (species %in% c("mouse", "rat")) {
    cat(sprintf("  [INFO] 检测到物种为 %s，正在将基因名转换为大写以匹配 Hallmark GMT...\n", species))
    
    # 1. 提取矩阵
    counts <- tryCatch({
      GetAssayData(seurat_obj, layer = "counts")
    }, error = function(e) {
      GetAssayData(seurat_obj, slot = "counts")
    })
    
    data <- tryCatch({
      GetAssayData(seurat_obj, layer = "data")
    }, error = function(e) {
      GetAssayData(seurat_obj, slot = "data")
    })
    
    # 2. 基因名转大写
    new_genes <- toupper(rownames(counts))
    
    # 3. 处理重复基因名
    if (any(duplicated(new_genes))) {
      cat("  [WARN] 转换后发现重复基因名，将保留第一个出现的探针...\n")
      keep_idx <- !duplicated(new_genes)
      counts <- counts[keep_idx, ]
      data <- data[keep_idx, ]
      new_genes <- new_genes[keep_idx]
    }
    
    rownames(counts) <- new_genes
    rownames(data) <- new_genes
    
    # 4. 创建新对象
    new_obj <- CreateSeuratObject(counts = counts, meta.data = seurat_obj@meta.data)
    new_obj <- SetAssayData(new_obj, slot = "data", new.data = data)
    Idents(new_obj) <- Idents(seurat_obj)
    
    # 5. 复制降维信息
    if (length(seurat_obj@reductions) > 0) {
      for (red in names(seurat_obj@reductions)) {
        new_obj@reductions[[red]] <- seurat_obj@reductions[[red]]
      }
    }
    return(new_obj)
  }
  return(seurat_obj)
}

# 3. 主分析函数 (已修改支持亚群) --------------------------------------------

run_irgsea_workflow <- function(
    rds_path,
    gmt_path,
    species = "human",
    group_by = "seurat_clusters",
    output_dir = "./irGSEA_Result",
    top_n_heatmap = 50,
    n_cores = 10,
    # --- 新增亚群参数 ---
    subset_mode = FALSE,   # 是否开启亚群过滤
    subset_col = NULL,     # 过滤依据的列 (如 "orig.ident")
    subset_idents = NULL   # 保留的ID (如 "Con")
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始 irGSEA Hallmark 分析流程\n")
  
  # --- 1. 读取数据 ---
  cat("  [DATA] 读取 Seurat 对象...\n")
  pbmc <- readRDS(rds_path)
  
  # --- 1.1 亚群过滤 (新增逻辑) ---
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("亚群模式下必须指定 subset_col 和 subset_idents")
    cat(sprintf("  [FILTER] 正在筛选亚群: %s 属于 {%s} ...\n", subset_col, paste(subset_idents, collapse=",")))
    
    # 临时切换 Ident 进行筛选
    Idents(pbmc) <- pbmc@meta.data[[subset_col]]
    pbmc <- subset(pbmc, idents = subset_idents)
    
    cat(sprintf("  [INFO] 筛选后剩余细胞数: %d\n", ncol(pbmc)))
  }
  
  # --- 1.2 设置主分析分组 ---
  # 检查分组
  if (!group_by %in% colnames(pbmc@meta.data)) stop(paste("未找到分组列:", group_by))
  Idents(pbmc) <- pbmc@meta.data[[group_by]]
  cat(sprintf("  [INFO] 分析分组: %s (共 %d 个亚群)\n", group_by, length(levels(Idents(pbmc)))))
  
  # --- 2. 基因名适配 ---
  pbmc <- convert_to_human_symbol(pbmc, species)
  
  # --- 3. 读取 GMT ---
  cat("  [GMT] 读取基因集...\n")
  gene_sets <- read_gmt_custom(gmt_path)
  
  # 检查基因重叠情况
  all_genes <- rownames(pbmc)
  gmt_genes <- unique(unlist(gene_sets))
  overlap <- intersect(all_genes, gmt_genes)
  cat(sprintf("  [CHECK] 数据集与GMT的重叠基因数: %d / %d (GMT总数)\n", length(overlap), length(gmt_genes)))
  
  if (length(overlap) < 100) stop("错误：基因重叠太少！请检查物种设置或基因名格式。")
  
  # --- 4. irGSEA 评分 (Scoring) ---
  cat("  [RUN] 开始 irGSEA 评分 (Method: AUCell, UCell, singscore, ssgsea)...\n")
  pbmc_irGSEA <- irGSEA.score(
    object      = pbmc,
    assay       = "RNA",
    slot        = "data",
    geneset     = gene_sets, 
    custom      = TRUE,       
    msigdb      = FALSE,      
    method      = c("AUCell", "UCell", "singscore", "ssgsea"),
    species     = "Homo sapiens", 
    kcdf        = 'Gaussian',
    ncores      = n_cores,
    seeds       = 123
  )
  
  saveRDS(pbmc_irGSEA, file.path(output_dir, "pbmc_irGSEA_scored.rds"))
  
  # --- 5. 差异整合分析 (Integration) ---
  cat("  [RUN] 整合差异结果 (RRA)...\n")
  # 注意：这里 group.by 使用的是函数传入的 group_by 参数
  result_dge <- irGSEA.integrate(
    object   = pbmc_irGSEA,
    group.by = group_by,
    method   = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  
  saveRDS(result_dge, file.path(output_dir, "irGSEA_DGE_result.rds"))
  
  # --- 6. 可视化模块 ---
  cat("  [PLOT] 生成可视化图表...\n")
  plot_dir <- file.path(output_dir, "Plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
  save_p <- function(p, name, w=10, h=8) {
    pdf(file.path(plot_dir, name), width=w, height=h)
    print(p)
    dev.off()
  }
  
  try({
    p_heat <- irGSEA.heatmap(object = result_dge, method = "RRA", top = top_n_heatmap, show.geneset = NULL)
    save_p(p_heat, "01_RRA_Global_Heatmap.pdf", w = 12, h = 10)
  })
  
  try({
    p_bubble <- irGSEA.bubble(object = result_dge, method = "RRA", top = 10)
    save_p(p_bubble, "02_RRA_Group_Bubble.pdf", w = 14, h = 8)
  })
  
  try({
    p_bar <- irGSEA.barplot(object = result_dge, method = c("AUCell", "UCell", "singscore", "ssgsea"))
    save_p(p_bar, "03_Method_Barplot.pdf", w = 12, h = 8)
  })
  
  if (!is.null(result_dge$RRA)) {
    top_pathway <- result_dge$RRA$Name[1]
    cat(sprintf("  [PLOT] 正在绘制 Top 1 通路: %s 的密度图...\n", top_pathway))
    
    try({
      p_density <- irGSEA.density.scatterplot(
        object = pbmc_irGSEA,
        method = "UCell",
        show.geneset = top_pathway,
        reduction = "umap"
      )
      save_p(p_density, paste0("04_Density_UMAP_", top_pathway, ".pdf"), w = 8, h = 7)
    })
    
    try({
      p_vln <- irGSEA.halfvlnplot(object = pbmc_irGSEA, method = "AUCell", show.geneset = top_pathway)
      save_p(p_vln, paste0("05_HalfVln_", top_pathway, ".pdf"), w = 8, h = 6)
    })
  }
  
  cat("=================================================\n")
  cat(sprintf("分析完成！结果目录: %s\n", output_dir))
  return(result_dge)
}


# 4. 参数设置与运行 (场景选择) ----------------------------------------------

RDS_PATH  <- "/home/lin/c_group/hep.rds"
GMT_PATH  <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"
SPECIES   <- "rat"
GROUP_COL <- "cellcluster"  # 无论全局还是亚群，最终都按这个聚类画图
CORES     <- 20

# =======================================================
# 场景 A: 全局分析 (分析所有细胞)
# =======================================================
cat("\n>>> 运行场景 A: 全局分析...\n")
res_global <- run_irgsea_workflow(
  rds_path      = RDS_PATH,
  gmt_path      = GMT_PATH,
  species       = SPECIES,
  group_by      = GROUP_COL,
  output_dir    = "/home/lin/c_group/subpopulation_analysis/gsea_global", # 输出到 global 文件夹
  top_n_heatmap = 50,
  n_cores       = CORES,
  subset_mode   = FALSE   # 关闭亚群过滤
)

# =======================================================
# 场景 B: 亚群分析 (仅分析 Con 组)
# =======================================================
cat("\n>>> 运行场景 B: 亚群分析 (Con 组)...\n")
res_sub <- run_irgsea_workflow(
  rds_path      = RDS_PATH,
  gmt_path      = GMT_PATH,
  species       = SPECIES,
  group_by      = GROUP_COL,
  output_dir    = "/home/lin/c_group/subpopulation_analysis/gsea_sub_Con", # 输出到 sub_Con 文件夹
  top_n_heatmap = 20,
  n_cores       = CORES,
  
  # --- 亚群核心设置 ---
  subset_mode   = TRUE,             # 开启亚群过滤
  subset_col    = "orig.ident",     # 过滤依据的列 (你的分组列)
  subset_idents = "Con"             # 仅保留 Con 组
)
