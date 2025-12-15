#############################################################
# 单细胞 irGSEA 参数化分析脚本 (v4.1: 修复版)
# 功能：多方法评分、差异通路识别、RRA聚合、可视化
# 修复：修正了 SetAssayData 参数拼写错误 (new_data -> new.data)
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
    
    # 1. 提取矩阵 (兼容 Seurat v4/v5 写法)
    # 尝试用 layer，失败则用 slot
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
    
    # 赋给新矩阵
    rownames(counts) <- new_genes
    rownames(data) <- new_genes
    
    # 4. 创建新对象 (重建以确保索引正确)
    new_obj <- CreateSeuratObject(counts = counts, meta.data = seurat_obj@meta.data)
    
    # 【修复点】这里参数名必须是 new.data (点号)，不能是下划线
    new_obj <- SetAssayData(new_obj, slot = "data", new.data = data)
    
    # 恢复 Idents
    Idents(new_obj) <- Idents(seurat_obj)
    
    # 5. 复制降维信息 (UMAP/TSNE 坐标)
    # 虽然基因名变了导致 Loadings 可能失效，但 Cell Embeddings (坐标) 是基于细胞的，画图需要用
    if (length(seurat_obj@reductions) > 0) {
      for (red in names(seurat_obj@reductions)) {
        new_obj@reductions[[red]] <- seurat_obj@reductions[[red]]
      }
    }
    
    return(new_obj)
  }
  return(seurat_obj)
}

# 3. 主分析函数 ------------------------------------------------------------

run_irgsea_workflow <- function(
    rds_path,
    gmt_path,
    species = "human",
    group_by = "seurat_clusters",
    output_dir = "./irGSEA_Result",
    top_n_heatmap = 50,  # 热图展示多少条通路
    n_cores = 10         # 线程数
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始 irGSEA Hallmark 分析流程\n")
  
  # --- 1. 读取数据 ---
  cat("  [DATA] 读取 Seurat 对象...\n")
  pbmc <- readRDS(rds_path)
  
  # 检查分组
  if (!group_by %in% colnames(pbmc@meta.data)) stop(paste("未找到分组列:", group_by))
  Idents(pbmc) <- pbmc@meta.data[[group_by]]
  cat(sprintf("  [INFO] 分析分组: %s (共 %d 个亚群)\n", group_by, length(levels(Idents(pbmc)))))
  
  # --- 2. 基因名适配 ---
  pbmc <- convert_to_human_symbol(pbmc, species)
  
  # --- 3. 读取 GMT ---
  cat("  [GMT] 读取基因集...\n")
  gene_sets <- read_gmt_custom(gmt_path)
  cat(sprintf("  [INFO] 成功加载 %d 条通路\n", length(gene_sets)))
  
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
    species     = "Homo sapiens", # 这里填Homo sapiens是因为我们已经手动处理了基因名
    kcdf        = 'Gaussian',
    ncores      = n_cores,
    seeds       = 123
  )
  
  saveRDS(pbmc_irGSEA, file.path(output_dir, "pbmc_irGSEA_scored.rds"))
  
  # --- 5. 差异整合分析 (Integration) ---
  cat("  [RUN] 整合差异结果 (RRA)...\n")
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
  
  # 辅助保存函数
  save_p <- function(p, name, w=10, h=8) {
    pdf(file.path(plot_dir, name), width=w, height=h)
    print(p)
    dev.off()
  }
  
  # 6.1 全局热图 (您的核心需求)
  # RRA 热图：展示综合排名靠前的通路
  try({
    p_heat <- irGSEA.heatmap(object = result_dge, method = "RRA", top = top_n_heatmap, show.geneset = NULL)
    save_p(p_heat, "01_RRA_Global_Heatmap.pdf", w = 12, h = 10)
  })
  
  # 6.2 气泡图 (展示每组显著富集的通路)
  try({
    p_bubble <- irGSEA.bubble(object = result_dge, method = "RRA", top = 10) # 每组展示 Top 10
    save_p(p_bubble, "02_RRA_Group_Bubble.pdf", w = 14, h = 8)
  })
  
  # 6.3 棒棒糖图
  try({
    p_bar <- irGSEA.barplot(object = result_dge, method = c("AUCell", "UCell", "singscore", "ssgsea"))
    save_p(p_bar, "03_Method_Barplot.pdf", w = 12, h = 8)
  })
  
  # 6.4 密度散点图 (Top Pathway Density Plot)
  # 自动提取综合排名第一的通路画图
  if (!is.null(result_dge$RRA)) {
    top_pathway <- result_dge$RRA$Name[1]
    cat(sprintf("  [PLOT] 正在绘制 Top 1 通路: %s 的密度图...\n", top_pathway))
    
    try({
      p_density <- irGSEA.density.scatterplot(
        object = pbmc_irGSEA,
        method = "UCell", # 使用 UCell 分数展示密度
        show.geneset = top_pathway,
        reduction = "umap" # 默认使用 umap
      )
      save_p(p_density, paste0("04_Density_UMAP_", top_pathway, ".pdf"), w = 8, h = 7)
    })
    
    # 6.5 半小提琴图
    try({
      p_vln <- irGSEA.halfvlnplot(object = pbmc_irGSEA, method = "AUCell", show.geneset = top_pathway)
      save_p(p_vln, paste0("05_HalfVln_", top_pathway, ".pdf"), w = 8, h = 6)
    })
  }
  
  cat("=================================================\n")
  cat(sprintf("分析完成！结果目录: %s\n", output_dir))
  return(result_dge)
}


# 4. 参数设置与运行 --------------------------------------------------------

# --- 核心参数设置 ---
RDS_PATH      <- "/home/lin/c_group/hep.rds"
GMT_PATH      <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"
SPECIES       <- "rat"              # 选项: "human", "mouse", "rat" (如果是rat会自动转大写匹配hallmark)
GROUP_COL     <- "cellcluster"      # 分组列
OUTPUT_DIR    <- "/home/lin/c_group/subpopulation_analysis/gsea"
TOP_N_HEATMAP <- 20                 # 热图展示多少条通路

# --- 执行分析 ---
res <- run_irgsea_workflow(
  rds_path      = RDS_PATH,
  gmt_path      = GMT_PATH,
  species       = SPECIES,
  group_by      = GROUP_COL,
  output_dir    = OUTPUT_DIR,
  top_n_heatmap = TOP_N_HEATMAP,
  n_cores       = 20
)
